#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2021  Elija Feigl
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see https://www.gnu.org/licenses/gpl-3.0.html.

import inspect
import logging
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from shutil import copyfile
from typing import List, Optional

from .. import get_resource
from ..core.utils import _exec, _get_executable
from ..link.linker import Linker
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
GSCALE_PRE = 0.3  # gscale for initial docking
DIEL_CONST = 1  # dielectric constant (enrgMD==1)
RES_SPAN = 24  # resolution range covered by low pass filtering
MAPTHRES = 0.0  # threshold for cropping mrc data (voxel smaller than)
TS_ENRG = 18000  # steps for initial enrgMD
TS_RELAX = 12000  # steps for relax enrgMD
###############################################################################


@dataclass
class Cascade(object):
    """ Cascade simulation class
    """
    conf: Path
    top: Path
    exb: Path
    json: Path
    seq: Path
    mrc: Path
    generated_with_mrdna: bool = True
    grid_pdb: Optional[Path] = None

    def __post_init__(self) -> None:
        self.prefix: str = self.top.stem
        self.logger = logging.getLogger(__name__)
        Path("run").mkdir(exist_ok=True)
        if self.grid_pdb is not None:
            copyfile(self.grid_pdb, "run/grid.pdb")

        self._split_exb_file()
        self.charmrun = _get_executable("charmrun")
        self.namd2 = _get_executable("namd2")
        self.vmd = _get_executable("vmd")

    def _split_exb_file(self) -> None:
        self.logger.info(
            "assuming annotated, sorted .exb files (mrDNA > march 2021)")
        # TODO: alternative splittings
        # TODO: include custom extrabonds
        # TODO: skip if exists
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

        with Path(f"run/{self.prefix}-BP.exb").open(mode='w') as f:
            for bond_list in exb_split:
                if bond_list[0].lower().startswith("# pair"):
                    f.writelines(bond_list)

        with Path(f"run/{self.prefix}-SR.exb").open(mode='w') as f:
            for bond_list in exb_split:
                if not bond_list[0].lower().startswith("# pushbonds"):
                    f.writelines(bond_list)

    def run_cascaded_fitting(self, base_time_steps: int, resolution: float,
                             is_SR=False, is_film=False, include_ss=False):
        def create_namd_file(namd_file, ts, ms=0, mdff="1", i_temp=300, f_temp=300):
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
                    set ITEMP {i_temp}
                    set FTEMP {f_temp}
                    """)
                f.write("\n".join([namd_header, namd_parameters, namd_base]))
            return (time_steps_last + ts)

        def _regular(step, folder):
            namd_file = Path(f"run/{prefix}_c-mrDNA-MDff-{step}.namd")
            time_steps_last = create_namd_file(namd_file, ts=base_time_steps)
            self.logger.debug(
                f"{step}: ts={base_time_steps}, mdff={gscale}, exb={enrgmd_file}, grid={grid_file}")
            previous_folder = self._run_namd(
                folder=folder, namd_file=namd_file)
            step += 1
            return time_steps_last, previous_folder, step

        def _annealing_incr(step, folder, gscale=0.01):
            # increase T with smaller gscale
            namd_file = Path(f"run/{prefix}_c-mrDNA-MDff-{step}.namd")

            time_steps_last = create_namd_file(
                namd_file, ts=ts_relax, i_temp=300, f_temp=400)
            self.logger.debug(
                f"{step}-{cascade}-{n}: annealing ts={TS_RELAX}, mdff={gscale}, annealing heating")
            previous_folder = self._run_namd(
                folder=folder, namd_file=namd_file)
            step += 1
            return time_steps_last, previous_folder, step

        def _annealing_decr(step, folder, gscale=0.01):
            # decrease T with smaller gscale
            namd_file = Path(f"run/{prefix}_c-mrDNA-MDff-{step}.namd")
            time_steps_last = create_namd_file(
                namd_file, ts=ts_relax, i_temp=400, f_temp=300)
            self.logger.debug(
                f"{step}-{cascade}-{n}: ts={TS_RELAX},  mdff={gscale}, annealing cooling")
            previous_folder = self._run_namd(
                folder=folder, namd_file=namd_file)
            step += 1
            return time_steps_last, previous_folder, step

        # TODO: create individual dataclass for each step to allow easy setting passing
        prefix = self.prefix
        n_cascade = N_CASCADE
        n_repeat = N_REPEAT
        first_stop = FIRST_STOP
        longrange_stop = LR_STOP
        gscale = GSCALE_PRE
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

        grid_pdb = "grid.pdb"
        if not Path("run/grid-base.dx").exists():
            self.logger.debug("Running VMD prep for generating cascade maps.")
            self._vmd_prep(resolution=resolution, n_cascade=n_cascade)
        if not Path(f"run/{grid_pdb}").exists():
            self.logger.debug("Generate custom gridPdb that excludes ssDNA.")
            linker = Linker(conf=self.conf, top=self.top, json=self.json, seq=self.seq,
                            generated_with_mrdna=self.generated_with_mrdna)
            linker.write_custom_gridpdb(
                dest=Path(f"run/{grid_pdb}"), exclude_ss=(not include_ss))
        init = 1  # init on first run
        step = 0
        time_steps_last = 0
        previous_folder = "none"

        if not self.generated_with_mrdna:  # pure enrgMD run without map
            folder = "enrgMD"
            grid_file = "none"
            enrgmd_file = "none"
            namd_file = Path(f"run/{prefix}_c-mrDNA-MDff-enrgMD.namd")
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

                # for bad docking the system will be relaxed reduced gscale and increased temperatur
                if cascade == 1 and n == 0:
                    enrgmd_file = f"../{prefix}.exb"
                    grid_file = "grid-0.dx"
                    folder = f"{step}"
                    time_steps_last, previous_folder, step = _annealing_incr(
                        step, folder)
                    folder = f"{step}"
                    time_steps_last, previous_folder, step = _annealing_decr(
                        step, folder)

                if cascade == 0 and n > first_stop:
                    break

                if n > longrange_stop:
                    enrgmd_file = f"{prefix}-SR.exb"
                    gscale = GSCALE
                else:
                    enrgmd_file = f"../{prefix}.exb"
                grid_file = f"grid-{n}.dx"
                folder = f"{step}"
                time_steps_last, previous_folder, step = _regular(step, folder)
                init = 0

        # TODO: add additional annealing step?

        # refine by removing intrahelical bonds
        grid_file = "grid-base.dx"
        if not is_SR:
            enrgmd_file = f"{prefix}-BP.exb"
            folder = f"{step}"
            time_steps_last, previous_folder, step = _regular(step, folder)

        # energy minimization with increased gscale to relax bonds
        gscale = 1.0
        folder = "final"
        namd_file = Path(f"run/{prefix}_c-mrDNA-MDff-final.namd")

        time_steps_last = create_namd_file(
            namd_file, ts=0, ms=base_time_steps)
        self.logger.debug(
            f"{folder}: ts={base_time_steps}, mdff={gscale}, exb={enrgmd_file}, {grid_file}")
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
        exec_folder = Path(f"run/{folder}")
        if not exec_folder.is_dir():
            exec_folder.mkdir()
            cmd = f"{self.charmrun} +p32 {self.namd2} +netpoll {namd_file}".split()
            self.logger.info(f"cascade: step {folder} with {cmd}")
            logfile = Path(f"run/{folder}/namd-out.log")
            _exec(cmd, logfile)
            if not any(exec_folder.glob('*.xst')):
                self.logger.warn(
                    "No .xst data written -> NAMD2 critical error. Abort cascade")
                sys.exit(1)
        else:
            self.logger.info(f"Skipping existing cascade: {folder}")
        return folder

    def _vmd_prep(self, n_cascade, resolution):
        vmd_prep = Path("run/mdff-prep.vmd")
        lines = ["package require volutil", "package require mdff"]

        lines.append(
            f"volutil -clamp {MAPTHRES}:1.0 {self.mrc} -o run/base.dx")

        n_layers = n_cascade - 1
        for n in range(n_cascade):
            gfilter = (n * (resolution - RES_SPAN)) / \
                n_layers + RES_SPAN - resolution
            lines.append(
                f"volutil -smooth {gfilter} run/base.dx -o run/{n}.dx")
            lines.append(f"mdff griddx -i run/{n}.dx -o run/grid-{n}.dx")
        lines.append("mdff griddx -i run/base.dx -o run/grid-base.dx")
        lines.append("exit")

        with vmd_prep.open(mode='w') as f:
            f.write('\n'.join(lines))
        cmd = f"{self.vmd} -dispdev text -e {vmd_prep}".split()
        self.logger.info(f"vmd prep:  with {cmd}")
        logfile = Path("run/vmd_prep-out.log")
        _exec(cmd, logfile)

    def _vmd_post(self, step):
        vmd_post = Path("run/mdff-post.vmd")
        lines = ["package require volutil", "package require mdff"]
        lines.append(f"mol new {self.top}")
        name = self.top.stem
        # full dcd
        if Path("./run/enrgMD").is_dir():
            lines.append(
                f"mol addfile ./run/enrgMD/{name}.dcd start 0 step 1 waitfor all")
        for n in range(step):
            lines.append(
                f"mol addfile ./run/{n}/{name}.dcd start 0 step 1 waitfor all")
        if Path("./run/final").is_dir():
            lines.append(
                f"mol addfile ./run/final/{name}.dcd start 0 step 1 waitfor all")

        lines.append(f"animate write dcd {name}.dcd")
        # last frame pdb
        lines.append("set sel [atomselect top all]")
        lines.append(f"$sel writepdb {name}-last.pdb")
        lines.append("exit")
        with vmd_post.open(mode='w') as f:
            f.write('\n'.join(lines))

        cmd = f"{self.vmd} -dispdev text -e {vmd_post}".split()
        self.logger.info(f"vmd postprocess:  with {cmd}")
        logfile = Path("run/vmd_post-out.log")
        _exec(cmd, logfile)
