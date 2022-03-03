#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
""" mrdna driven cascade fitting simulation class.
    fitting pipeline 12.2021:
        * prep files using VMD:
            - apply N_CASCADE gaussian filters to map of range (0, RES_SPAN-map_resolution)
            - [TODO-IMPROVEMENT: for >15A maps use binary maps instead]
            - create grid.pdb [or read supplied grid.pdb]
        * [TS_ENRG pure enrgMD if not mrDNA]
        * N_REPEAT cascades
            * N_CASCADE steps (n)
                - 1. cascade and n > FIRST_STOP ends: cascade early
                - 2. cascade: simulated annealing
                - n > LR_STOP: only short range extrabonds
        * [ single step switching to base pair extrabonds:]
        * energy minimization with increased gscale
        * processing files using VMD:
            - merge step-trajectories
            - create pdb file for last frame

    film pipeline 01.2021:
        - shorter timesteps, increased dcd output
"""
import inspect
from dataclasses import dataclass

from .. import get_resource


###############################################################################
# PRESET PARAMETERS
N_CASCADE = 8  # number of steps in the cascade
N_REPEAT = 2  # number of consecutive cascades
FIRST_STOP = 2  # abort first cascade after as many steps
LR_STOP = 3  # turn off long range enrgMD bonds after as many steps

RES_SPAN = 24  # resolution range covered by low pass filtering
MAPTHRES = 0.0  # threshold for cropping mrc data (voxel smaller than)
###############################################################################


@dataclass
class CascadeData:
    prefix: str

    n_cascade = N_CASCADE
    n_repeat = N_REPEAT
    first_stop = FIRST_STOP
    longrange_stop = LR_STOP

    gscale = 0.3  # scaling factor for map potential
    gscale_min = 1.0  # scaling factor for map potential for energy minimization
    gscale_anneal = 0.01
    temperature = 300
    dielectric_constant = 1  # dielectric constant (enrgMD==1)
    ts = 12000  # number of time steps (default)
    ts_relax = 12000  # steps for relax enrgMD and anneling
    ms_relax = 12000  # minimization steps for relax enrgMD and anneling
    namd_resource = "namd.txt"

    is_initialized: bool = False
    step_counter: int = 0
    time_steps_counter: int = 0
    previous_time_steps_counter: int = 0
    folder: str = "none"
    previous_folder: str = "none"
    grid_pdb: str = "grid.pdb"

    def set_folder(self, new_folder: str) -> None:
        self.previous_folder = self.folder
        self.folder = new_folder

    def update_time_steps_counter(self, incr: int) -> None:
        self.previous_time_steps_counter = self.time_steps_counter
        self.time_steps_counter += incr

    def movie_setting(self) -> None:
        self.ts_relax = 2400
        self.ms_relax = 1200
        self.namd_resource = "namd_film.txt"

    def create_namd_file(
        self, namd_path, ts, ms=0, gscale=None, grid_file=None, exb_file=None, annealing=None
    ):

        is_fitting = False if grid_file is None else True
        gscale = self.gscale if gscale is None else gscale
        annealing = (self.temperature, self.temperature) if annealing is None else annealing
        init = 0 if self.is_initialized else 1
        i_temp, f_temp = annealing
        mdff = "1" if is_fitting else "0"

        with namd_path.open(mode="w") as file:
            namd_base = get_resource(self.namd_resource).read_text()
            namd_header = get_resource("namd_header.txt").read_text()
            namd_parameters = inspect.cleandoc(
                f"""
                set INIT {init}
                set PREFIX {self.prefix}
                set TS {ts}
                set MS {ms}
                set GRIDON {mdff}
                set DIEL {self.dielectric_constant}
                set GSCALE {gscale}
                set GRIDFILE {grid_file}
                set GRIDPDB {self.grid_pdb}

                set ENRGMDON on
                set ENRGMDBONDS {exb_file}

                set TSLAST {self.previous_time_steps_counter}
                set CURRENT {self.folder}
                set PREVIOUS {self.previous_folder}
                set ITEMP {i_temp}
                set FTEMP {f_temp}
                """
            )
            file.write("\n".join([namd_header, namd_parameters, namd_base]))

        self.update_time_steps_counter(incr=ts)
        self.is_initialized = True
