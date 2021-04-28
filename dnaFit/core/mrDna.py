import logging
import os
import sys
from pathlib import Path
from shutil import copyfile, copytree

import MDAnalysis as mda
import numpy as np

from ..data.mrc import recenter_mrc
from .utils import _exec, _get_executable

logger = logging.getLogger(__name__)


def run_mrDNA(cad_file: Path, seq_file: Path, prefix: str, directory: str = "mrDNA", gpu: int = 0, multidomain=False):
    """ running a mrDNA simulation by executing mrDNA externally
            all files are automatically written into a folder "mrDNA"
            checks of completion of mrDNA run
    """
    home_directory = os.getcwd()
    try:
        Path(directory).mkdir(parents=True, exist_ok=True)
        os.chdir(directory)
        logger.debug(f"changing directory to: {os.getcwd()}")

        mrDNA = _get_executable("mrdna")

        input_mrDNA = f"-o {prefix} -d {directory} -g {gpu} --run-enrg-md"
        if multidomain:
            input_mrDNA += " --coarse-steps 3e7 --crossover-to-intrahelical-cutoff 25 --coarse-bond-cutoff 300 "
        input_mrDNA += f" --sequence-file {seq_file} {cad_file}"

        cmd = (str(mrDNA), input_mrDNA)
        _exec(cmd)

        logger.info("mrDNA: finished. Checking for final files.")
        is_ok = (os.path.isfile(f"./{prefix}-3.psf")
                 and os.path.isfile(f"./{prefix}-3.pdb")
                 and os.path.isfile(f"./{prefix}-3.exb"))
        if not is_ok:
            logger.error("At least one of mrdna's *-3. files was not created")
            sys.exit(1)

    finally:
        os.chdir(home_directory)
        logger.debug(f"changing directory to: {os.getcwd()}")


def prep_cascaded_fitting(prefix: str, cad_file: Path, seq_file: Path, mrc_file: Path):
    """ prep cascaded fitting:
            prepares a new folder "dnaFit" with files from a finished mrDNA run
            in subfolder "mrDNA" and copies necessary files.
    """
    home_directory = os.getcwd()
    try:
        Path("dnaFit").mkdir(parents=True, exist_ok=True)
        os.chdir("dnaFit")
        logger.debug(f"changing directory to: {os.getcwd()}")

        copyfile(cad_file, f"./{prefix}.json")
        copyfile(seq_file, f"./{prefix}.seq")

        copyfile(f"../mrDNA/{prefix}-3.psf", f"./{prefix}.psf")
        copyfile(f"../mrDNA/{prefix}-3.exb", f"./{prefix}.exb")
        copyfile(f"../mrDNA/{prefix}-3.pdb", f"./{prefix}.pdb")

        mrc_shift = recenter_mrc(mrc_file, apply=False)
        recenter_conf(top=Path(f"./{prefix}.psf"),
                      conf=Path(f"./{prefix}.pdb"), to_position=mrc_shift)

        copytree("../mrDNA/charmm36.nbfix", "./charmm36.nbfix")
    except Exception:
        logger.exception(
            f"mrDNA: failed to copy mrDNA files to working directory {os.getcwd()}")
        raise Exception
    finally:
        os.chdir(home_directory)
        logger.debug(f"changing directory to: {os.getcwd()}")


def recenter_conf(top: Path, conf: Path, to_position=np.array([0.0, 0.0, 0.0])) -> None:
    u = mda.Universe(str(top), str(conf))
    translation = to_position - u.atoms.center_of_geometry()
    u.atoms.translate(translation)
    u.atoms.write(str(conf))
