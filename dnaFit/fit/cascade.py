#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
import logging
import multiprocessing
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from shutil import copyfile
from typing import List
from typing import Optional
from typing import Tuple

from ..core.utils import _exec
from ..core.utils import _get_executable
from ..link.linker import Linker
from .atomic_model_fit import AtomicModelFit
from .cascade_data import CascadeData
from .cascade_data import MAPTHRES
from .cascade_data import RES_SPAN

warnings.filterwarnings("ignore")


@dataclass
class Cascade:
    """Cascade simulation class"""

    conf: Path
    top: Path
    exb: Path
    json: Path
    seq: Path
    mrc: Path

    generated_with_mrdna: bool = True
    custom_grid_pdb: Optional[Path] = None

    def __post_init__(self) -> None:
        self.logger = logging.getLogger(__name__)

        self.prefix: str = self.top.stem
        self.data = CascadeData(self.prefix)

        run_folder = Path("run")
        if run_folder.exists():
            self.logger.info("Re-using existing folder 'run' for cascade data.")
        else:
            run_folder.mkdir(exist_ok=True)
        if self.custom_grid_pdb is not None:
            self.logger.info("Using existing grid.pdb file from 'run' folder.")
            copyfile(self.custom_grid_pdb, "run/grid.pdb")

        self._split_exb_file()
        try:
            self.namd2 = _get_executable("namd2")
            self.vmd = _get_executable("vmd")
        except OSError as exc:
            self.logger.exception("Abort. %s", exc)
            sys.exit(1)

    def _split_exb_file(self) -> None:
        self.logger.info("assuming annotated, sorted .exb files (mrdna > march 2021)")
        # TODO-IMPROVEMENT: include custom extrabonds

        exb_bp = Path(f"{self.prefix}-BP.exb")
        exb_sr = Path(f"{self.prefix}-SR.exb")
        if exb_bp.exists() and exb_sr.exists():
            self.logger.info("Split .exb file exist. Skip recalculation.")
            return

        with self.exb.open(mode="r") as exb_file:
            exb_data = exb_file.readlines()
        exb_split: List[List[str]] = list()
        bond_list: List[str] = list()
        for line in exb_data:
            if line.startswith("#"):
                exb_split.append(bond_list)
                bond_list = list()
            bond_list.append(line)
        exb_split.pop(0)

        with exb_bp.open(mode="w") as exb_file_bp:
            for bond_list in exb_split:
                if bond_list[0].lower().startswith("# pair"):
                    exb_file_bp.writelines(bond_list)

        with exb_sr.open(mode="w") as exb_file_sr:
            for bond_list in exb_split:
                if not bond_list[0].lower().startswith("# pushbonds"):
                    exb_file_sr.writelines(bond_list)

    def _prep_gridpdb(self, resolution, exclude_ss) -> str:
        grid_pdb = "grid.pdb"
        if not Path("run/grid-base.dx").exists():
            self.logger.debug("Running VMD prep for generating cascade maps.")
            self._vmd_prep(resolution=resolution, n_cascade=self.data.n_cascade)
        if not Path(f"run/{grid_pdb}").exists():
            if exclude_ss:
                self.logger.info("Generate custom gridPdb that excludes ssDNA.")
            else:
                self.logger.info("Generate custom gridPdb that includes ssDNA.")
            linker = Linker(
                conf=self.conf,
                top=self.top,
                json=self.json,
                seq=self.seq,
                generated_with_mrdna=self.generated_with_mrdna,
            )
            linker.write_custom_gridpdb(dest=Path(f"run/{grid_pdb}"), exclude_ss=exclude_ss)
        return grid_pdb

    def _run_namd(self, namd_file, num_procs=None):
        exec_folder = Path(f"run/{self.data.folder}")
        if not exec_folder.is_dir():
            exec_folder.mkdir()
            if num_procs is None:
                num_procs = max(1, multiprocessing.cpu_count() - 1)

            cmd = f"{self.namd2} +p{num_procs} {namd_file}".split()
            self.logger.info("cascade: folder %s with %s", self.data.folder, cmd)
            logfile = Path(f"run/{self.data.folder}/namd-out.log")
            _exec(cmd, logfile)
            if not any(exec_folder.glob("*.xst")):
                self.logger.warning("No .xst data written -> NAMD2 critical error. Abort cascade")
                sys.exit(1)
        else:
            self.logger.info("Skipping existing cascade: %s", self.data.folder)

    def _relax(self, exb_file: str):
        """pure enrgMD run without map"""
        namd_file = Path(f"run/{self.prefix}_c-mrdna-MDff-enrgMD.namd")
        self.data.create_namd_file(
            namd_path=namd_file,
            ts=self.data.ts_relax,
            ms=self.data.ms_relax,
            exb_file=exb_file,
        )
        self.logger.debug(
            "%s: gscale=%s, ts=%s, ms=%s, exb=%s, grid=%s",
            self.data.folder,
            0,
            self.data.ts_relax,
            self.data.ms_relax,
            exb_file,
            "off",
        )
        self._run_namd(namd_file=namd_file)
        # NOTE: no self.data.step_counter increase (only count fitting)

    def _step(self, grid_file: str, exb_file: str):
        """regular cascade step"""
        namd_filepath = Path(f"run/{self.prefix}_c-mrdna-MDff-{self.data.step_counter}.namd")
        self.data.create_namd_file(
            namd_path=namd_filepath,
            ts=self.data.ts,
            gscale=self.data.gscale,
            grid_file=grid_file,
            exb_file=exb_file,
        )
        self.logger.debug(
            "%s: gscale=%s, ts=%s, ms=%s, exb=%s, grid=%s",
            self.data.folder,
            self.data.gscale,
            self.data.ts,
            0,
            exb_file,
            grid_file,
        )
        self._run_namd(namd_file=namd_filepath)
        self.data.step_counter += 1

    def _annealing(self, grid_file: str, exb_file: str, temperatures: Tuple[int, int]):
        """namd annealing step"""
        namd_filepath = Path(f"run/{self.prefix}_c-mrdna-MDff-{self.data.step_counter}.namd")
        self.data.create_namd_file(
            namd_path=namd_filepath,
            ts=self.data.ts_relax,
            grid_file=grid_file,
            exb_file=exb_file,
            annealing=temperatures,
            gscale=self.data.gscale_anneal,
        )
        self.logger.debug(
            "%s: gscale=%s, ts=%s, ms=%s, exb=%s, grid=%s, anneal=>%s",
            self.data.folder,
            self.data.gscale_anneal,
            self.data.ts_relax,
            0,
            exb_file,
            grid_file,
            temperatures[-1],
        )
        self._run_namd(namd_file=namd_filepath)
        self.data.step_counter += 1

    def _emin(self, grid_file: str, exb_file: str):
        """energy minimization step"""
        namd_filepath = Path(f"run/{self.prefix}_c-mrdna-MDff-final.namd")
        self.data.create_namd_file(
            namd_path=namd_filepath,
            ts=0,
            ms=self.data.ts,
            grid_file=grid_file,
            exb_file=exb_file,
            gscale=self.data.gscale_min,
        )
        self.logger.debug(
            "%s: gscale=%s, ts=%s, ms=%s, exb=%s, grid=%s",
            self.data.folder,
            self.data.gscale_min,
            0,
            self.data.ts,
            exb_file,
            grid_file,
        )
        self._run_namd(namd_file=namd_filepath)
        # NOTE: no self.data.step_counter increase

    def _fitting(self, is_sr: bool) -> None:

        self.logger.debug("Start VMD based cascade MDff prep.")
        exb = f"../{self.prefix}"

        # PART 1: regular relaxation run
        if not self.generated_with_mrdna:
            self.data.set_folder("relax")
            self._relax(exb_file=f"{exb}.exb")

        # PART 2: fitting cascade
        for cascade in range(self.data.n_repeat):
            for step in range(self.data.n_cascade):
                # PART 2.1: short & coarse 1st cascade for rough alignment
                if cascade == 0 and step > self.data.first_stop:
                    break

                # PART 2.2: annealing with decreased gscale (2nd cascade)
                if cascade == 1 and step == 0:
                    self.data.set_folder(f"{self.data.step_counter}")
                    self._annealing(
                        grid_file=f"grid-{step}.dx",
                        exb_file=f"{exb}.exb",
                        temperatures=(self.data.temperature, 400),
                    )
                    self.data.set_folder(f"{self.data.step_counter}")
                    self._annealing(
                        grid_file=f"grid-{step}.dx",
                        exb_file=f"{exb}.exb",
                        temperatures=(400, self.data.temperature),
                    )

                # PART 2.3: elastic network reduction
                exb_type = "-SR" if step > self.data.longrange_stop else ""
                self.data.set_folder(f"{self.data.step_counter}")
                self._step(grid_file=f"grid-{step}.dx", exb_file=f"{exb}{exb_type}.exb")

        # (PART 2.4): elastic network reduction
        # TODO-IMPROVEMENT: add additional annealing step?

        # (PART 2.5): remove all but base pair elastic network
        if not is_sr:
            exb_type = "-BP"
            self.data.set_folder(f"{self.data.step_counter}")
            self._step(grid_file="grid-base.dx", exb_file=f"{exb}{exb_type}.exb")

        # PART 3: energy minimization with increased gscale to relax bonds
        self.data.set_folder("final")
        self._emin(grid_file="grid-base.dx", exb_file=f"{exb}{exb_type}.exb")

    def run_cascaded_fitting(
        self, time_steps: int, resolution: float, is_sr=False, is_film=False, exclude_ss=True
    ):
        if is_film:
            self.data.movie_setting()
        self.data.grid_pdb = self._prep_gridpdb(resolution=resolution, exclude_ss=exclude_ss)
        self.data.ts = time_steps

        self._fitting(is_sr=is_sr)

        self.logger.debug("VMD based cascade MDff prep.")
        self._vmd_post(self.data.step_counter)

        final_conf = self.conf.with_name(f"{self.prefix}-last.pdb")
        if not final_conf.is_file():
            self.logger.critical("Cascaded fit incomplete. %s not found. Abort.", final_conf)
            sys.exit(1)
        return AtomicModelFit(conf=final_conf, top=self.top, mrc=self.mrc)

    def _vmd_prep(self, n_cascade: int, resolution: float) -> None:
        vmd_prep_path = Path("run/mdff-prep.vmd")
        lines = ["package require volutil", "package require mdff"]

        # convert mrc file to dx format
        lines.append(f"volutil -clamp {MAPTHRES}:1.0 {self.mrc} -o run/base.dx")

        # create grid_file for cascade
        n_layers = n_cascade - 1
        for step in range(n_cascade):
            gfilter = (step * (resolution - RES_SPAN)) / n_layers + RES_SPAN - resolution
            lines.append(f"volutil -smooth {gfilter} run/base.dx -o run/{step}.dx")
            # TODO-IMPROVEMENT: binary maps for resolution > 15A ?
            lines.append(f"mdff griddx -i run/{step}.dx -o run/grid-{step}.dx")
        lines.append("mdff griddx -i run/base.dx -o run/grid-base.dx")
        lines.append("exit")
        with vmd_prep_path.open(mode="w") as vmd_file:
            vmd_file.write("\n".join(lines))

        cmd = f"{self.vmd} -dispdev text -e {vmd_prep_path}".split()
        self.logger.info("vmd prep: with %s", cmd)
        logfile = Path("run/vmd_prep-out.log")
        _exec(cmd, logfile)

    def _vmd_post(self, n_steps: int) -> None:
        vmd_post_path = Path("run/mdff-post.vmd")
        lines = ["package require volutil", "package require mdff"]
        lines.append(f"mol new {self.top}")
        name = self.top.stem

        # full concatenated dcd
        if Path("./run/enrgMD").is_dir():
            lines.append(f"mol addfile ./run/relax/{name}.dcd start 0 step 1 waitfor all")
        for step in range(n_steps):
            lines.append(f"mol addfile ./run/{step}/{name}.dcd start 0 step 1 waitfor all")
        if Path("./run/final").is_dir():
            lines.append(f"mol addfile ./run/final/{name}.dcd start 0 step 1 waitfor all")
        lines.append(f"animate write dcd {name}.dcd")

        # last frame pdb
        lines.append("set sel [atomselect top all]")
        lines.append(f"$sel writepdb {name}-last.pdb")
        lines.append("exit")
        with vmd_post_path.open(mode="w") as vmd_file:
            vmd_file.write("\n".join(lines))

        cmd = f"{self.vmd} -dispdev text -e {vmd_post_path}".split()
        self.logger.info("vmd postprocess:  with %s", cmd)
        logfile = Path("run/vmd_post-out.log")
        _exec(cmd, logfile)
