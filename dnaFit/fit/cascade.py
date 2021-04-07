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

""" mrDNA driven cascade fitting simulation class.
"""
###############################################################################
# PRESET PARAMETERS
N_CASCADE = 8  # number of steps in the cascade
N_REPEAT = 2  # number of consecutive cascades
FIRST_STOP = 2  # abort first cascade after as many steps
LR_STOP = 3  # turn off long range enrgMD bonds after as many steps
GSCALE = 0.3  # scaling factor for map potential
DIEL_CONST = 1  # dielectric constant (enrgMD==1)
RES_SPAN = 14  # resolution range covered by low pass filtering
MAPTHRES = 0.0  # threshold for cropping mrc data (voxel smaller than)
TS_ENRG = 18000  # steps for initial enrgMD
TS_RELAX = 18000  # steps for relax enrgMD
###############################################################################


def external_docking_loop(prefix):
    ##########################################################################
    # TODO -low-: implement automatic alignment or python UI based alignment
    # BREAK: reorient helix and fit still external
    conf = Path(f"./{prefix}.pdb")
    while True:
        input(
            f"You need to manually align mrc and pdb using VMD. save as {conf} and press ENTER")
        if conf.is_file():
            return conf


def _exec(cmd):
    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, universal_newlines=True)
    for line in process.stdout:
        sys.stdout.write(line)
        sys.stdout.flush()


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

        # internally moving to center at origin simplifies rotation in vmd
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

        with self.exb.with_name(f"{self.prefix}-BP.exb").open(mode='w') as f:
            for bond_list in exb_split:
                if bond_list[0].startswith("# PAIR"):
                    f.writelines(bond_list)

        with self.exb.with_name(f"{self.prefix}-SR.exb").open(mode='w') as f:
            for bond_list in exb_split:
                if not bond_list[0].startswith("# PUSHBONDS"):
                    f.writelines(bond_list)

    def run_cascaded_fitting(self, base_time_steps: int, resolution: float, is_SR=False, is_film=False):
        def run_namd():
            # TODO: skipping folders that already are complete
            cmd = f"{self.charmrun} + p32 {self.namd2} + netpoll $1 2 > &1 | tee {output_name}.log"
            self.logger.info(f"cascade:  with {cmd}")
            _exec(cmd)
            return output_name

        def create_namd_file(namd_file, ts, ms=0, mdff="1"):
            with namd_file.open(mode='w') as f:
                if not is_film:
                    namd_base = get_resource("namd.txt").read_text()
                else:
                    namd_base = get_resource("namd_film.txt").read_text()
                namd_header = get_resource("namd_header.txt").read_text()
                namd_parameters = inspect.cleandoc(f"""
                    set PREFIX {prefix}
                    set TS {ts}
                    set MS {ms}
                    set GRIDON {mdff}
                    set DIEL{dielectr_constant}
                    set GSCALE {gscale}
                    set GRIDFILE {grid_file}
                    set GRIDPDB {grid_pdb}

                    set ENRGMDON on
                    set ENRGMDBONDS {enrgmd_file}

                    set OUTPUTNAME{output_name}
                    set TSLAST {time_steps_last}
                    set N {step}
                    set PREVIOUS {folder_last}
                    """)
                f.write("\n".join([namd_header, namd_parameters, namd_base]))
            return (time_steps_last + ts)

        def vmd_prep():
            vmd_prep = Path("mdff-prep.vmd")
            grid_pdb = Path("grid.pdb")
            lines = ["package require volutil", "package require mdff"]

            lines.append(
                f"volutil -clamp {MAPTHRES}:1.0 {self.mrc} -o base.dx")
            n_layers = n_cascade - 1
            for n in range(n_cascade):
                gfilter = RES_SPAN/n_layers*(n_layers-n) + resolution
                lines.append(
                    f"volutil -smooth {gfilter} base.dx -o {n}.dx")
                lines.append(f"mdff griddx -i {n}.dx -o grid-{n}.dx")
            lines.append(
                f"mdff gridpdb -psf {self.top} -pdb {self.conf} -o {grid_pdb}")
            lines.append("exit")

            with Path("mdff-prep.vmd").open() as f:
                f.writelines(lines)
            cmd = f"{self.vmd} -dispdev text -eofexit -e {vmd_prep}"
            self.logger.info(f"vmd prep:  with {cmd}")
            _exec(cmd)
            return grid_pdb

        def vmd_post():
            vmd_post = Path("mdff-post.vmd")
            lines = ["package require volutil", "package require mdff"]
            lines.append(f"mol new {self.top}")
            name = self.top.stem
            # full dcd
            lines.append(
                f"mol addfile ./enrgMD/{name}.dcd start 0 step 1 waitfor all")
            for n in range(n_cascade):
                lines.append(
                    f"mol addfile ./{n}/{name}.dcd start 0 step 1 waitfor all")
            lines.append(
                f"mol addfile ./final/{name}.dcd start 0 step 1 waitfor all")

            lines.append(f"animate write dcd {name}.dcd")
            # last frame pdb
            lines.append("set sel [atomselect top all]")
            lines.append(f"$sel writepdb {name}-last.pdb")
            lines.append("exit")
            with vmd_post.open() as f:
                f.writelines(lines)

            cmd = f"{self.vmd} -dispdev text -eofexit -e {vmd_post}"
            self.logger.info(f"vmd postprocess:  with {cmd}")
            _exec(cmd)

        prefix = self.prefix
        n_cascade = N_CASCADE
        n_repeat = N_REPEAT
        first_stop = FIRST_STOP
        longrange_stop = LR_STOP
        gscale = GSCALE
        dielectr_constant = DIEL_CONST

        if is_film:
            ts_enrg = 2400
            ms_enrg = 1200
            ts_relax = 9600
            base_time_steps = 2400
        else:
            ts_enrg = TS_ENRG
            ms_enrg = TS_ENRG
            ts_relax = TS_RELAX

        self.logger.debug("VMD based cascade MDff prep.")
        grid_pdb = vmd_prep()
        namd_file = self.top.with_name(f"{prefix}_c-mrDNA-MDff.namd")

        # pure enrgMD run without map
        step = -1
        time_steps_last = 0
        folder_last = ""
        output_name = "enrgMD"
        time_steps_last = create_namd_file(
            namd_file, ts=ts_enrg, ms=ms_enrg, mdff="0")
        self.logger.debug(f"{output_name}: ts={ts_enrg}, ms={ms_enrg}, mdff=0")
        folder_last = run_namd()

        # cascades
        for cascade in range(n_repeat):
            for n in range(n_cascade):
                # for bad docking the system will be relaxed with pure enrgMD after first iteration
                if cascade == 1 and n == 0:
                    step += 1
                    enrgmd_file = f"{prefix}.exb"
                    output_name = f"{step}/{prefix}"
                    time_steps_last = create_namd_file(
                        namd_file, ts=ts_relax, mdff="0")
                    self.logger.debug(
                        f"{step}-{cascade}-{n}: ts={TS_RELAX}, mdff=0, enrgMD relax")
                    folder_last = run_namd()

                if cascade == 1 and n > first_stop:
                    break

                if n > longrange_stop:
                    enrgmd_file = f"{prefix}-SR.exb"
                else:
                    enrgmd_file = f"{prefix}.exb"
                grid_file = f"grid-{n}.dx"
                step += 1
                output_name = f"{step}/{prefix}"
                time_steps_last = create_namd_file(
                    namd_file, ts=base_time_steps)
                self.logger.debug(
                    f"{step}-{cascade}-{n}: ts={base_time_steps}, mdff=1, {enrgmd_file}, {grid_file}")
                folder_last = run_namd()

        # refine by removing intrahelical bonds
        if not is_SR:
            step += 1
            enrgmd_file = f"{prefix}-BP.exb"
            time_steps_last = create_namd_file(namd_file, ts=base_time_steps)
            self.logger.debug(
                f"{step}: ts={base_time_steps}, mdff=1, {enrgmd_file}, {grid_file}")
            folder_last = run_namd()

        # energy minimization with increased gscale to relax bonds
        gscale = 1.0
        output_name = "final"
        time_steps_last = create_namd_file(
            namd_file, ts=0, ms=base_time_steps)
        self.logger.debug(
            f"{output_name}: ts={base_time_steps}, mdff=1, {enrgmd_file}, {grid_file}")
        _ = run_namd()

        self.logger.debug("VMD based cascade MDff prep.")
        vmd_post()
        # TODO: cleanup??

        final_conf = self.conf.with_name(f"{self.prefix}-last.pdb")
        if not final_conf.is_file():
            self.logger.error(
                f"cascaded fit incomplete. {final_conf} not found ")
            sys.exit(1)
        # revert internal recentering to align with original mrc data
        self.logger.debug(
            f"shifted mrc & pdb file back using previous shift: {self.mrc_shift}")
        _ = recenter_mrc(self.mrc, to_position=self.mrc_shift)
        self._recenter_conf(conf=final_conf, to_position=self.mrc_shift)
        return AtomicModelFit(
            conf=final_conf, top=self.top, mrc=self.mrc)
