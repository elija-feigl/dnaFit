import sys
import attr
from pathlib import Path
from typing import Any
import subprocess
import numpy as np
import mrcfile
import MDAnalysis as mda

from ..core.utils import _get_executable
from . import templates

try:
    import importlib.resources as resources
except ImportError:
    import importlib_resources as resources

""" DESCR:
    mrDNA driven cascade fitting simulation class.
"""


def external_docking_loop(prefix):
    ##########################################################################
    # TODO -low-: implement automatic alignmant or python UI based alignment
    # BREAK: reorient helix and fit still external
    conf = Path("./{}.pdb".format(prefix))
    while True:
        input(
            "You need to manualy align mrc and pdb using VMD. save as {} and press ENTER".format(conf))
        if conf.is_file():
            return conf


@attr.s
class Cascade(object):
    """ Cascade simulation class
    """
    conf: Path = attr.ib()
    top: Path = attr.ib()
    exb: Path = attr.ib()
    mrc: Path = attr.ib()
    recenter: bool = attr.ib()
    is_docked: bool = attr.ib()

    def __attrs_post_init__(self) -> None:
        self.prefix: str = self.top.stem

        _ = self._recenter_mrc()
        self._recenter_conf(conf=self.conf)

        if not self.is_docked:
            self.conf = external_docking_loop(self.prefix)

        # NOTE: splitting here allows starting from an old mrDNA or enrgMD setup.
        #   as long as .exb is orderer (extra script)
        self._split_exb_file()

        self.charmrun = _get_executable("charmrun")
        self.namd2 = _get_executable("namd2")
        self.vmd = _get_executable("vmd")

    def _recenter_mrc(self, reverse=False) -> Any:
        with mrcfile.open(self.mrc, mode='r+') as mrc:
            c = np.array(mrc.header["cella"])
            cell = np.array([c["x"], c["y"], c["z"]])
        if reverse:
            shift = cell * 0.5 if reverse else cell * - 0.5
        mrc.header["origin"] = tuple(shift)
        return shift

    def _recenter_conf(self, conf: Path, to_position=np.array([0.0, 0.0, 0.0])) -> None:
        u = mda.Universe(str(self.top), str(conf))
        translation = to_position - u.atoms.center_of_geometry
        u.atoms.translate(translation)
        u.atoms.write(str(conf))

    def _split_exb_file(self) -> None:
        # NOTE: assuming modifed mrDNA > march 2021, (annotated, sorted .exb files)
        with self.exb.open(mode='r') as f:
            exb_data = f.readlines()
        exb_split = list()
        bond_list = list()
        for line in exb_data:
            if line.startswith("#"):
                exb_split.append(bond_list)
                bond_list = list()
            bond_list.append(line)

        with self.exb.with_name("{}-HO.exb".format(self.prefix)).open(mode='w') as f:
            for bond_list in exb_split:
                if bond_list[0].startswith("# PAIR"):
                    f.writelines(bond_list)

        with self.exb.with_name("{}-SR.exb".format(self.prefix)).open(mode='w') as f:
            for bond_list in exb_split:
                if not bond_list[0].startswith("# PUSHBONDS"):
                    f.writelines(bond_list)

    def run_cascaded_fitting(self, base_time_steps: int, resolution: float):
        # NOTE: assuming we can put all files in the current folder!

        def create_sh_file(sh_file):
            with sh_file.open(mode='w') as f:
                namd_path = self.namd2.resolve().parent
                vmd_path = self.vmd.resolve().parent
                sh_base = resources.read_text(
                    templates, "c-mrDNA-MDff-cascade-sh.txt")
                sh_parameters = "readonly UBIN = {namd_path}\n\
                    readonly VBIN = {vmd_path}\n\n\
                    # general\n\
                    readonly DESIGNNAME = {prefix}\n\
                    declare - ri TSTEPS = {time_steps}\n\n\
                    # cascade\n\
                    declare - ri NCASCADE = {n_cascade}\n\
                    readonly GFMAX = {resolution_max}\n\
                    readonly MAPRESOLUTION = {resolution}\n".format(**locals())
                f.write("\n".join([sh_parameters, sh_base]))

        def create_namd_file(namd_file):
            with namd_file.open(mode='w') as f:
                namd_base = resources.read_text(templates, "namd.txt")
                namd_header = resources.read_text(
                    templates, "namd_header.txt")
                namd_parameters = "set PREFIX {prefix}\n\
                    set TS {time_steps}\n\
                    set MS {minimisation_steps}\n\
                    set GRIDON {grid_on}\n\
                    set DIEL{dielectr_constant}\n\
                    set GSCALE {gscale}\n\
                    set GRIDFILE {grid_file}\n\
                    set GRIDPDB {grid_pdb}\n\n\
                    set ENRGMDON {enrgmd_on}\n\
                    set ENRGMDBONDS {enrgmd_file}\n\n\
                    set OUTPUTNAME{output_name}\n".format(**locals())

                f.write("\n".join([namd_header, namd_parameters, namd_base]))
            return

        # TODO: allow additional parameter changes
        #       TODO: change from energymin to fixed.pdb protocol
        prefix = self.prefix
        time_steps = base_time_steps
        minimisation_steps = base_time_steps
        n_cascade = 8
        resolution_max = resolution + 14
        gscale = 0.3
        dielectr_constant = 1
        grid_on = "1"
        grid_file = "1.dx"
        grid_pdb = "grid.pdb"
        enrgmd_on = "on"
        enrgmd_file = "$PREFIX.exb"
        output_name = "$N/$PREFIX"

        namd_file = self.top.with_name(
            "{}_c-mrDNA-MDff.namd".format(self.prefix))
        create_namd_file(namd_file)

        ######################################################################
        # VERSION 1: execute with sh file
        # NOTE: sh execute changes namd-file
        sh_file = self.top.with_name(
            "{}_c-mrDNA-MDff.sh".format(self.prefix))
        create_sh_file(sh_file)

        cmd = ("sh", sh_file)
        # TODO: logger: "Starting cascade with sh script {}".format(cmd)

        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, universal_newlines=True)
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()

        ######################################################################
        # VERSION 1: pytohnised fitting script
        # TODO: ...

        final_conf = self.conf.with_name(
            "{}-final.pdb".format(self.prefix))
        if final_conf.is_file():
            raise Exception("ERROR: cascaded fit incomplete")

        mrc_shift = self._recenter_mrc(reverse=True)
        self._recenter_conf(conf=final_conf, to_position=mrc_shift)
        return AtomicModelFit(
            conf=final_conf, top=self.top, mrc=self.mrc)


@attr.s
class AtomicModelFit(object):
    conf: Path = attr.ib()
    top: Path = attr.ib()
    mrc: Path = attr.ib()

    def write_linkage(self, cad_file: Path, seq_file: Path):
        # TODO: write persistent and human readable linkage (Fbp, FidDid)
        raise NotImplementedError

    def write_output(self, write_mmCif=True, crop_mrc=True):
        # TODO: save final pdb/mmcif together with mrc and masked mrc
        raise NotImplementedError
