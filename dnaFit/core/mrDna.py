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


def run_mrDNA(cad_file: Path, seq_file: Path, prefix: str, directory="mrdna", gpu: int = 0, multidomain=False):
    """ running a mrDNA simulation by executing mrDNA externally
            all files are automatically written into a folder "mrDNA"
            checks of completion of mrDNA run
    """
    mrDNA = _get_executable("mrdna")

    cmd = [str(mrDNA), "-o", str(prefix), "-d",
           str(directory), "-g", str(gpu), "--run-enrg-md"]
    if multidomain:
        cmd += ["--coarse-steps", "3e7", "--crossover-to-intrahelical-cutoff",
                "25", "--coarse-bond-cutoff", "300"]
    export_idx = 4 if multidomain else 3
    cmd += ["--sequence-file", str(seq_file), str(cad_file)]
    logger.info(f"starting mrDNA: creates folder ./{directory}.")
    _exec(cmd)

    logger.info("mrDNA: finished. Checking for final files.")
    is_ok = (os.path.isfile(f"./{directory}/{prefix}-{export_idx}.psf")
             and os.path.isfile(f"./{directory}/{prefix}-{export_idx}.pdb")
             and os.path.isfile(f"./{directory}/{prefix}-{export_idx}.exb"))
    if not is_ok:
        logger.error(
            f"At least one of mrdna's *-{export_idx}. files was not created")
        sys.exit(1)


def prep_cascaded_fitting(prefix: str, cad_file: Path, seq_file: Path,
                          mrc_file: Path, directory="mrdna", multidomain=False):
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
        copyfile(mrc_file, f"./{prefix}.mrc")
        export_idx = 4 if multidomain else 3
        pre_mrdna = f"../{directory}/{prefix}-{export_idx}"
        copyfile(f"{pre_mrdna}.psf", f"./{prefix}.psf")
        copyfile(f"{pre_mrdna}.exb", f"./{prefix}.exb")

        # NOTE: pdb is start not finish. use output/NAME.coor instead
        coor = f"../{directory}/output/{prefix}-{export_idx}-1.coor"
        u = mda.Universe(f"./{prefix}.psf", coor)
        u.atoms.write(f"./{prefix}.pdb")

        mrc_shift = recenter_mrc(mrc_file, apply=False)
        recenter_conf(top=Path(f"./{prefix}.psf"),
                      conf=Path(f"./{prefix}.pdb"), to_position=mrc_shift)

        copytree(f"../{directory}/charmm36.nbfix", "./charmm36.nbfix")
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
