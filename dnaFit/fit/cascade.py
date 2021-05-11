import inspect
import logging
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import List

from .. import get_resource
from ..core.utils import _get_executable, _exec
from .atomic_model_fit import AtomicModelFit

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
RES_SPAN = 24  # resolution range covered by low pass filtering
MAPTHRES = 0.0  # threshold for cropping mrc data (voxel smaller than)
TS_ENRG = 18000  # steps for initial enrgMD
TS_RELAX = 18000  # steps for relax enrgMD
###############################################################################


@dataclass
class Cascade(object):
    """ Cascade simulation class
    """
    conf: Path
    top: Path
    exb: Path
    mrc: Path

    def __post_init__(self) -> None:
        self.prefix: str = self.top.stem
        self.logger = logging.getLogger(__name__)
        self._split_exb_file()

        self.charmrun = _get_executable("charmrun")
        self.namd2 = _get_executable("namd2")
        self.vmd = _get_executable("vmd")

    def _split_exb_file(self) -> None:
        self.logger.info(
            "assuming annotated, sorted .exb files (mrDNA > march 2021)")
        # TODO: alternative splittings
        with self.exb.open(mode='r') as f:
            exb_data = f.readlines()
        exb_split: List[List[str]] = list()
        bond_list: List[str] = list()
        for line in exb_data:
            if line.startswith("#"):
                exb_split.append(bond_list)
                bond_list = list()
            bond_list.append(line)
        exb_split.pop(0)

        with self.exb.with_name(f"{self.prefix}-BP.exb").open(mode='w') as f:
            for bond_list in exb_split:
                if bond_list[0].lower().startswith("# pair"):
                    f.writelines(bond_list)

        with self.exb.with_name(f"{self.prefix}-SR.exb").open(mode='w') as f:
            for bond_list in exb_split:
                if not bond_list[0].lower().startswith("# pushbonds"):
                    f.writelines(bond_list)

    def run_cascaded_fitting(self, base_time_steps: int, resolution: float,
                             is_SR=False, is_film=False, is_enrgMD=False):
        def create_namd_file(namd_file, ts, ms=0, mdff="1"):
            with namd_file.open(mode='w') as f:
                if not is_film:
                    namd_base = get_resource("namd.txt").read_text()
                else:
                    namd_base = get_resource("namd_film.txt").read_text()
                namd_header = get_resource("namd_header.txt").read_text()
                namd_parameters = inspect.cleandoc(f"""
                    set INIT {init}
                    set PREFIX {prefix}
                    set TS {ts}
                    set MS {ms}
                    set GRIDON {mdff}
                    set DIEL {dielectr_constant}
                    set GSCALE {gscale}
                    set GRIDFILE {grid_file}
                    set GRIDPDB {grid_pdb}

                    set ENRGMDON on
                    set ENRGMDBONDS {enrgmd_file}

                    set TSLAST {time_steps_last}
                    set CURRENT {folder}
                    set PREVIOUS {previous_folder}
                    """)
                f.write("\n".join([namd_header, namd_parameters, namd_base]))
            return (time_steps_last + ts)

        # TODO: create individual dataclass for each step to allow easy setting passing
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

        grid_pdb = Path("grid.pdb")
        if not grid_pdb.exists():
            self._vmd_prep(resolution=resolution,
                           n_cascade=n_cascade,
                           grid_pdb=grid_pdb)
        init = 1  # init on first run
        step = 0
        time_steps_last = 0
        previous_folder = "none"

        if is_enrgMD:  # pure enrgMD run without map
            folder = "enrgMD"
            grid_file = "none"
            enrgmd_file = "none"
            namd_file = self.top.with_name(
                f"{prefix}_c-mrDNA-MDff-{folder}.namd")
            time_steps_last = create_namd_file(
                namd_file, ts=ts_enrg, ms=ms_enrg, mdff="0")
            self.logger.debug(
                f"{folder}: ts={ts_enrg}, ms={ms_enrg}, mdff=0")
            previous_folder = self._run_namd(
                folder=folder, namd_file=namd_file)
            init = 0

        # cascades
        for cascade in range(n_repeat):
            for n in range(n_cascade):
                # for bad docking the system will be relaxed with pure enrgMD after first iteration
                # TODO: better strong annealing instead?
                if cascade == 1 and n == 0:
                    enrgmd_file = f"{prefix}.exb"
                    folder = f"{step}"
                    namd_file = self.top.with_name(
                        f"{prefix}_c-mrDNA-MDff-{folder}.namd")

                    time_steps_last = create_namd_file(
                        namd_file, ts=ts_relax, mdff="0")
                    self.logger.debug(
                        f"{step}-{cascade}-{n}: ts={TS_RELAX}, mdff=0, enrgMD relax")
                    previous_folder = self._run_namd(
                        folder=folder, namd_file=namd_file)
                    step += 1

                if cascade == 0 and n > first_stop:
                    break

                if n > longrange_stop:
                    enrgmd_file = f"{prefix}-SR.exb"
                else:
                    enrgmd_file = f"{prefix}.exb"
                grid_file = f"grid-{n}.dx"
                folder = f"{step}"
                namd_file = self.top.with_name(
                    f"{prefix}_c-mrDNA-MDff-{folder}.namd")

                time_steps_last = create_namd_file(
                    namd_file, ts=base_time_steps)
                self.logger.debug(
                    f"{step}-{cascade}-{n}: ts={base_time_steps}, mdff=1, {enrgmd_file}, {grid_file}")
                previous_folder = self._run_namd(
                    folder=folder, namd_file=namd_file)
                step += 1
                init = 0

        # TODO: add annealing step?

        # refine by removing intrahelical bonds
        if not is_SR:
            enrgmd_file = f"{prefix}-BP.exb"
            folder = f"{step}"
            namd_file = self.top.with_name(
                f"{prefix}_c-mrDNA-MDff-{folder}.namd")

            time_steps_last = create_namd_file(namd_file, ts=base_time_steps)
            self.logger.debug(
                f"{step}: ts={base_time_steps}, mdff=1, {enrgmd_file}, {grid_file}")
            previous_folder = self._run_namd(
                folder=folder, namd_file=namd_file)
            step += 1

        # energy minimization with increased gscale to relax bonds
        gscale = 1.0
        folder = "final"
        namd_file = self.top.with_name(
            f"{prefix}_c-mrDNA-MDff-{folder}.namd")

        time_steps_last = create_namd_file(
            namd_file, ts=0, ms=base_time_steps)
        self.logger.debug(
            f"{folder}: ts={base_time_steps}, mdff=1, {enrgmd_file}, {grid_file}")
        _ = self._run_namd(folder=folder, namd_file=namd_file)

        self.logger.debug("VMD based cascade MDff prep.")
        self._vmd_post(step)

        final_conf = self.conf.with_name(f"{self.prefix}-last.pdb")
        if not final_conf.is_file():
            self.logger.error(
                f"cascaded fit incomplete. {final_conf} not found ")
            sys.exit(1)

        return AtomicModelFit(
            conf=final_conf, top=self.top, mrc=self.mrc)

    def _run_namd(self, namd_file, folder):
        if not Path(folder).is_dir():
            Path(folder).mkdir()
            cmd = f"{self.charmrun} +p32 {self.namd2} +netpoll {namd_file}".split()
            self.logger.info(f"cascade: step {folder} with {cmd}")
            _exec(cmd)
            if not any(Path(folder).iterdir()):
                Path(folder).rmdir()
                self.logger.error("No files written. Abort cascade")
                sys.exit(1)
        else:
            self.logger.info(f"skipping existing cascade: {folder}")
        return folder

    def _vmd_prep(self, n_cascade, resolution, grid_pdb):
        vmd_prep = Path("mdff-prep.vmd")
        lines = ["package require volutil", "package require mdff"]

        lines.append(
            f"volutil -clamp {MAPTHRES}:1.0 {self.mrc} -o base.dx")
        lines.append("mdff griddx -i base.dx -o grid-base.dx")

        n_layers = n_cascade - 1
        for n in range(n_cascade):
            gfilter = (n * (resolution - RES_SPAN))/n_layers + RES_SPAN
            lines.append(
                f"volutil -smooth {gfilter} base.dx -o {n}.dx")
            lines.append(f"mdff griddx -i {n}.dx -o grid-{n}.dx")
        lines.append(
            f"mdff gridpdb -psf {self.top} -pdb {self.conf} -o {grid_pdb}")
        # lines.append("exit")  # TODO: check if cause of file error

        with vmd_prep.open(mode='w') as f:
            f.write('\n'.join(lines))
        cmd = f"{self.vmd} -dispdev text -eofexit -e {vmd_prep}".split()
        self.logger.info(f"vmd prep:  with {cmd}")
        _exec(cmd)

    def _vmd_post(self, step):
        vmd_post = Path("mdff-post.vmd")
        lines = ["package require volutil", "package require mdff"]
        lines.append(f"mol new {self.top}")
        name = self.top.stem
        # full dcd
        if Path("./enrgMD").is_dir():
            lines.append(
                f"mol addfile ./enrgMD/{name}.dcd start 0 step 1 waitfor all")
        for n in range(step):
            lines.append(
                f"mol addfile ./{n}/{name}.dcd start 0 step 1 waitfor all")
        if Path("./final").is_dir():
            lines.append(
                f"mol addfile ./final/{name}.dcd start 0 step 1 waitfor all")

        lines.append(f"animate write dcd {name}.dcd")
        # last frame pdb
        lines.append("set sel [atomselect top all]")
        lines.append(f"$sel writepdb {name}-last.pdb")
        # lines.append("exit")  # TODO: check if cause of file error
        with vmd_post.open(mode='w') as f:
            f.write('\n'.join(lines))

        cmd = f"{self.vmd} -dispdev text -eofexit -e {vmd_post}".split()
        self.logger.info(f"vmd postprocess:  with {cmd}")
        _exec(cmd)
