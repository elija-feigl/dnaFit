#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
""" custom mrDNA setup
"""
import logging
import os
import sys
from pathlib import Path
from shutil import copyfile

import MDAnalysis as mda
import numpy as np

from ..data.mrc import get_mrc_box
from ..data.mrc import recenter_mrc
from ..data.mrc import write_mrc_from_atoms
from .utils import _exec
from .utils import _get_executable

logger = logging.getLogger(__name__)


def run_mrdna(
    cad_file: Path,
    seq_file: Path,
    prefix: str,
    directory="mrdna",
    gpu: int = 0,
    multidomain: bool = False,
    coarse_steps: float = 5e7,
    bond_cutoff: int = 300,
):
    """running a mrdna simulation by executing mrdna externally
    all files are automatically written into a folder "mrdna"
    checks of completion of mrdna run
    """
    try:
        mrdna = _get_executable("mrdna")
    except OSError as exc:
        logger.exception("Abort. %s", exc)
        sys.exit(1)

    cmd = [
        str(mrdna),
        "-o",
        str(prefix),
        "-d",
        str(directory),
        "-g",
        str(gpu),
        "--run-enrg-md",
        "--enrg-md-steps",
        "1e5",
    ]
    if multidomain:
        cmd += [
            "--coarse-steps",
            str(coarse_steps),
            "--crossover-to-intrahelical-cutoff",
            "25",
            "--coarse-bond-cutoff",
            str(bond_cutoff),
        ]
    export_idx = 4 if multidomain else 3
    cmd += ["--sequence-file", str(seq_file), str(cad_file)]
    logger.info("starting mrdna: creates folder ./%s.", directory)
    logfile = Path("mrdna-interal.log")
    _exec(cmd, logfile)

    logger.info("mrdna: finished. Checking for final files.")
    is_ok = (
        Path(f"./{directory}/{prefix}-{export_idx}.psf").is_file()
        and Path(f"./{directory}/{prefix}-{export_idx}.pdb").is_file()
        and Path(f"./{directory}/{prefix}-{export_idx}.exb").is_file()
    )
    if not is_ok:
        logger.error("At least one of mrdna's *-%s. files was not created", export_idx)
        sys.exit(1)


def prep_cascaded_fitting(
    prefix: str,
    cad_file: Path,
    seq_file: Path,
    mrc_file: Path,
    directory="mrdna",
    multidomain=False,
):
    """prep cascaded fitting:
    prepares a new folder "prep" with files from a finished mrdna run
    in subfolder "mrdna" and copies necessary files.
    """
    home_directory = Path.cwd()
    try:
        Path("prep").mkdir(parents=True, exist_ok=True)
        os.chdir("prep")
        logger.debug("changing directory to: %s", Path.cwd())

        copyfile(cad_file, f"./{prefix}.json")
        copyfile(seq_file, f"./{prefix}.seq")
        export_idx = 4 if multidomain else 3
        pre_mrdna = f"../{directory}/{prefix}-{export_idx}"
        copyfile(f"{pre_mrdna}.psf", f"./{prefix}.psf")
        copyfile(f"{pre_mrdna}.exb", f"./{prefix}.exb")

        # NAME.pdb is start not finish. use output/NAME.coor instead
        coor = f"../{directory}/output/{prefix}-{export_idx}-1.coor"
        if not Path(coor).is_file():
            logger.error("erngMD stopped with error. no coor file generated. Abort!")
            sys.exit(1)
        universe = mda.Universe(f"./{prefix}.psf", coor)

        logger.info("cropping mrc to minimal box to improve simulation efficiency")
        mrc_boxed = mrc_file.with_name(f"{mrc_file.stem}-boxed.mrc")
        # create and write boxed mrc
        write_mrc_from_atoms(
            path=mrc_file,
            atoms=universe.atoms,
            path_out=mrc_boxed,
            context=50.0,
            cut_box=True,
            keep_data=True,
        )
        copyfile(mrc_boxed, f"./{prefix}.mrc")

        # determine shift
        mrc_shift = recenter_mrc(mrc_boxed, apply=False)
        # shift atoms and write pdb
        translation = mrc_shift - universe.atoms.center_of_geometry()
        universe.atoms.translate(translation)
        # determine pdb box size from mrc
        universe.dimensions = get_mrc_box(mrc_boxed) + [90.0, 90.0, 90.0]
        universe.atoms.write(f"./{prefix}.pdb")

    except FileNotFoundError as exc:
        logger.error("mrdna: failed to copy mrdna files to working directory %s", Path.cwd())
        logger.exception(exc)
    finally:
        os.chdir(home_directory)
        logger.debug("changing directory to: %s", Path.cwd())


def recenter_conf(top: Path, conf: Path, to_position=np.array([0.0, 0.0, 0.0])) -> None:
    """recenter a set of atoms on a specific position by its center of geomerty."""
    universe = mda.Universe(str(top), str(conf))
    translation = to_position - universe.atoms.center_of_geometry()
    universe.atoms.translate(translation)
    universe.atoms.write(str(conf))
