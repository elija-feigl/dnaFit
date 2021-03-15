import sys
import attr
import inspect
import subprocess
import logging

from pathlib import Path
from typing import List

import numpy as np
import MDAnalysis as mda

from .atomic_model_fit import AtomicModelFit
from .. import get_resource
from ..core.utils import _get_executable
from ..data.mrc import recenter_mrc

import warnings
warnings.filterwarnings('ignore')

""" DESCR:
    mrDNA driven cascade fitting simulation class.
"""


def external_docking_loop(prefix):
    ##########################################################################
    # TODO -low-: implement automatic alignmant or python UI based alignment
    # BREAK: reorient helix and fit still external
    conf = Path(f"./{prefix}.pdb")
    while True:
        input(
            f"You need to manualy align mrc and pdb using VMD. save as {conf} and press ENTER")
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
        self.logger = logging.getLogger(__name__)

        # intenally moving to center at origin simplifies rotation in vmd

        self.logger.info("internal recenter at origin for rotation in vmd")
        self.mrc_shift = recenter_mrc(self.mrc)
        self.logger.debug(f"shifted mrc file by: {self.mrc_shift}")
        self._recenter_conf(conf=self.conf)

        if not self.is_docked:
            self.conf = external_docking_loop(self.prefix)

        self._split_exb_file()

        self.charmrun = _get_executable("charmrun")
        self.namd2 = _get_executable("namd2")
        self.vmd = _get_executable("vmd")

    def _recenter_conf(self, conf: Path, to_position=np.array([0.0, 0.0, 0.0])) -> None:
        u = mda.Universe(str(self.top), str(conf))
        translation = to_position - u.atoms.center_of_geometry()
        u.atoms.translate(translation)
        u.atoms.write(str(conf))

    def _split_exb_file(self) -> None:
        self.logger.info(
            "assuming annotated, sorted .exb files (mrDNA > march 2021)")
        with self.exb.open(mode='r') as f:
            exb_data = f.readlines()
        exb_split: List[List[str]] = list()
        bond_list: List[str] = list()
        for line in exb_data:
            if line.startswith("#"):
                exb_split.append(bond_list)
                bond_list = list()
            bond_list.append(line)

        with self.exb.with_name(f"{self.prefix}-HO.exb").open(mode='w') as f:
            for bond_list in exb_split:
                if bond_list[0].startswith("# PAIR"):
                    f.writelines(bond_list)

        with self.exb.with_name(f"{self.prefix}-SR.exb").open(mode='w') as f:
            for bond_list in exb_split:
                if not bond_list[0].startswith("# PUSHBONDS"):
                    f.writelines(bond_list)

    def run_cascaded_fitting(self, base_time_steps: int, resolution: float):
        """ creating files and externally executing sh-script for cascaded flexible fitting
        """
        def create_sh_file(sh_file):
            """ create enrgMD-driven cascaded flexible fitting script from template
            """
            with sh_file.open(mode='w') as f:
                namd_path = self.namd2.resolve().parent
                vmd_path = self.vmd.resolve().parent
                sh_base = get_resource(
                    "c-mrDNA-MDff-cascade-sh.txt").read_text()
                sh_parameters = inspect.cleandoc(f"""
                    readonly UBIN = {namd_path}
                    readonly VBIN = {vmd_path}

                    # general
                    readonly DESIGNNAME = {prefix}
                    declare - ri TSTEPS = {time_steps}

                    # cascade
                    declare - ri NCASCADE = {n_cascade}
                    readonly GFMAX = {resolution_max}
                    readonly MAPRESOLUTION = {resolution}
                    """)
                f.write("\n".join([sh_parameters, sh_base]))

        def create_namd_file(namd_file):
            with namd_file.open(mode='w') as f:
                namd_base = get_resource("namd.txt").read_text()
                namd_header = get_resource("namd_header.txt").read_text()
                namd_parameters = inspect.cleandoc(f"""
                    set PREFIX {prefix}
                    set TS {time_steps}
                    set MS {minimisation_steps}
                    set GRIDON {grid_on}
                    set DIEL{dielectr_constant}
                    set GSCALE {gscale}
                    set GRIDFILE {grid_file}
                    set GRIDPDB {grid_pdb}

                    set ENRGMDON {enrgmd_on}
                    set ENRGMDBONDS {enrgmd_file}

                    set OUTPUTNAME{output_name}
                    """)
                f.write("\n".join([namd_header, namd_parameters, namd_base]))
            return

        # TODO: allow additional parameter changes
        #       NOTE: change from energymin to fixed.pdb protocol
        #       NOTE: actually mrDNA includes enrgMD
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

        namd_file = self.top.with_name(f"{self.prefix}_c-mrDNA-MDff.namd")
        create_namd_file(namd_file)
        self.logger.debug(f"initializing namd file {namd_file}")

        ######################################################################
        # VERSION 1: execute with sh file
        # NOTE: sh execute changes namd-file
        sh_file = self.top.with_name(f"{self.prefix}_c-mrDNA-MDff.sh")
        self.logger.debug(f"changing sh file {sh_file}")
        create_sh_file(sh_file)

        cmd = ("sh", sh_file)
        self.logger.info(f"cascade:  with {cmd}")
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, universal_newlines=True)
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()

        ######################################################################
        # VERSION 2: pythonised fitting script
        # TODO: ...

        final_conf = self.conf.with_name(f"{self.prefix}-last.pdb")
        if not final_conf.is_file():
            self.logger.warning(
                f"cascaded fit incomplete. {final_conf} not found ")
            raise Exception("cascade incomplete")

        # revert internal recentering to align with original mrc data
        self.logger.debug(
            f"shifted mrc & pdb file back using previous shift: {self.mrc_shift}")
        _ = recenter_mrc(self.mrc, to_position=self.mrc_shift)
        self._recenter_conf(conf=final_conf, to_position=self.mrc_shift)
        return AtomicModelFit(
            conf=final_conf, top=self.top, mrc=self.mrc)
