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

""" custom mrDNA setup
"""

import logging
import os
import sys
from pathlib import Path
from shutil import copyfile

import MDAnalysis as mda
import numpy as np

from ..data.mrc import recenter_mrc
from .utils import _exec, _get_executable

logger = logging.getLogger(__name__)


def run_mrdna(cad_file: Path, seq_file: Path, prefix: str, directory="mrdna",
              gpu: int = 0, multidomain=False, bond_cutoff=300):
    """ running a mrdna simulation by executing mrdna externally
            all files are automatically written into a folder "mrdna"
            checks of completion of mrdna run
    """
    mrdna = _get_executable("mrdna")

    cmd = [str(mrdna), "-o", str(prefix), "-d",
           str(directory), "-g", str(gpu), "--run-enrg-md", "--enrg-md-steps", "1e5"]
    if multidomain:
        cmd += ["--coarse-steps", "5e7", "--crossover-to-intrahelical-cutoff",
                "25", "--coarse-bond-cutoff", str(bond_cutoff)]
    export_idx = 4 if multidomain else 3
    cmd += ["--sequence-file", str(seq_file), str(cad_file)]
    logger.info("starting mrdna: creates folder ./%s.", directory)
    logfile = Path("mrdna-interal.log")
    _exec(cmd, logfile)

    logger.info("mrdna: finished. Checking for final files.")
    is_ok = (os.path.isfile(f"./{directory}/{prefix}-{export_idx}.psf")
             and os.path.isfile(f"./{directory}/{prefix}-{export_idx}.pdb")
             and os.path.isfile(f"./{directory}/{prefix}-{export_idx}.exb"))
    if not is_ok:
        logger.error(
            "At least one of mrdna's *-%s. files was not created", export_idx)
        sys.exit(1)


def prep_cascaded_fitting(prefix: str, cad_file: Path, seq_file: Path,
                          mrc_file: Path, directory="mrdna", multidomain=False):
    """ prep cascaded fitting:
            prepares a new folder "prep" with files from a finished mrdna run
            in subfolder "mrdna" and copies necessary files.
    """
    home_directory = os.getcwd()
    try:
        Path("prep").mkdir(parents=True, exist_ok=True)
        os.chdir("prep")
        logger.debug("changing directory to: %s", os.getcwd())

        copyfile(cad_file, f"./{prefix}.json")
        copyfile(seq_file, f"./{prefix}.seq")
        copyfile(mrc_file, f"./{prefix}.mrc")
        export_idx = 4 if multidomain else 3
        pre_mrdna = f"../{directory}/{prefix}-{export_idx}"
        copyfile(f"{pre_mrdna}.psf", f"./{prefix}.psf")
        copyfile(f"{pre_mrdna}.exb", f"./{prefix}.exb")

        # NAME.pdb is start not finish. use output/NAME.coor instead
        coor = f"../{directory}/output/{prefix}-{export_idx}-1.coor"
        if not Path(coor).is_file():
            logger.error(
                "erngMD stopped with error. no coor file generated. Abort!")
            sys.exit(1)
        universe = mda.Universe(f"./{prefix}.psf", coor)

        mrc_shift = recenter_mrc(mrc_file, apply=False)
        translation = mrc_shift - universe.atoms.center_of_geometry()
        universe.atoms.translate(translation)
        universe.atoms.write(f"./{prefix}.pdb")

    except Exception as exc:
        logger.exception(
            "mrdna: failed to copy mrdna files to working directory %s", os.getcwd())
        raise Exception from exc
    finally:
        os.chdir(home_directory)
        logger.debug("changing directory to: %s", os.getcwd())


def recenter_conf(top: Path, conf: Path, to_position=np.array([0.0, 0.0, 0.0])) -> None:
    """ recenter a set of atoms on a specific position by its center of geomerty.
    """
    universe = mda.Universe(str(top), str(conf))
    translation = to_position - universe.atoms.center_of_geometry()
    universe.atoms.translate(translation)
    universe.atoms.write(str(conf))
