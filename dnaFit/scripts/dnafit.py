#! /usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os
import subprocess
import sys
from pathlib import Path
from shutil import copyfile, copytree

import click
from dnaFit import get_resource
from dnaFit.core.utils import _get_executable
from dnaFit.data.mrc import write_mrc_from_atoms
from dnaFit.fit.atomic_model_fit import AtomicModelFit
from dnaFit.fit.cascade import Cascade
from dnaFit.pdb.structure import Structure
from dnaFit.version import get_version

""" cascaded mrDNA-driven MD flexible fitting:
"""
logger = logging.getLogger(__name__)


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(get_version())
    ctx.exit()


@click.group()
@click.option('--version', is_flag=True, callback=print_version,
              expose_value=False, is_eager=True)
def cli():
    pass


def run_mrDNA(cad_file: Path, seq_file: Path, prefix: str, directory: str = "mrDNA", gpu: int = 0, multidomain=False):
    """ running a mrDNA simulation by executing mrDNA externally
            all files are automatically written into a foler "mrDNA"
            checks of completion of mrDNA run
    """
    home_directory = os.getcwd()
    try:
        Path(directory).mkdir(parents=True, exist_ok=True)
        os.chdir(directory)
        logger.debug(f"changing directory to: {os.getcwd()}")

        mrDNA = _get_executable("mrdna")

        input_mrDNA = f"-o {prefix} -d {directory} -g {gpu} --run-enrg-md --sequence-file "
        if multidomain:
            input_mrDNA += "--coarse-steps 3e7 --crossover-to-intrahelical-cutoff 25 --coarse-bond-cutoff 300 "
        input_mrDNA += f"{seq_file} {cad_file}"

        input_mrDNA
        cmd = (str(mrDNA), input_mrDNA)

        logger.info(f"mrDNA: with {cmd}")
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, universal_newlines=True)
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()

        logger.info("mrDNA: finished. Checking for final files.")
        is_ok = (os.path.isfile(f"./{prefix}-3.psf")
                 and os.path.isfile(f"./{prefix}-3.pdb")
                 and os.path.isfile(f"./{prefix}-3.exb"))
        if not is_ok:
            logger.error("At least one of mrdna's *-3. files was not created")
            raise Exception("mrDNA incomplete")
    finally:
        os.chdir(home_directory)
        logger.debug(f"changing directory to: {os.getcwd()}")


def prep_cascaded_fitting(prefix: str, cad_file: Path, seq_file: Path):
    """ prep cascaded fitting:
            prepares a new folder "dnaFit" with files from a finished mrDNA run
            in fubfolder "mrDNA" and copyies necessary files.
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
        copyfile(f"../mrDNA/{prefix}-3.pdb", f"./{prefix}-undocked.pdb")
        copytree("../mrDNA/charmm36.nbfix", "./charmm36.nbfix")
    except:
        logger.exception(
            f"mrDNA: failed to copy mrDNA files to working directory {os.getcwd()}")
        raise Exception("Failed to copy mrDNA files")
    finally:
        os.chdir(home_directory)
        logger.debug(f"changing directory to: {os.getcwd()}")


@cli.command()
@click.argument('cadnano', type=click.Path(exists=True))
@click.argument('mrc', type=click.Path(exists=True))
@click.argument('sequence', type=click.Path(exists=True))
@click.option('-g', '--gpu', type=int, default=0, help='GPU used for simulation', show_default=True)
@click.option('-o', '--output-prefix', 'prefix', type=str, default=None,
              help="short design name, default to json name")
@click.option('--timesteps', type=int, default=12000,
              help='timesteps per cascade (multiple of 12)')
@click.option('--resolution', type=float, default=10.0,
              help='mrc map resolution in Angstrom')
@click.option('--multidomain', is_flag=True,
              help='multidomain structures require different settings for equilibration')
def mrDnaFit(cadnano, mrc, sequence, gpu, prefix, timesteps, resolution, multidomain):
    """mrDNA simulation of CADNANO design file followed by cascaded
        mrDNA-driven MD flexible fitting to MRC cryo data

        CADNANO is the name of the design file [.json]\n
        SEQUENCE is the scaffold strand sequence file [.txt, .seq]\n
        MRC is the name of the cryo EM volumetric data file [.mrc]\n
    """
    # NOTE: click will drop python2 support soon and actually return a Path
    cad_file = Path(cadnano).resolve()
    mrc_file = Path(mrc).resolve()
    seq_file = Path(sequence).resolve()

    prefix = cad_file.stem if prefix is None else prefix

    run_mrDNA(cad_file, seq_file, prefix, directory="mrDNA",
              gpu=gpu, multidomain=multidomain)
    prep_cascaded_fitting(prefix, cad_file, seq_file)

    # TODO: -low- parse & set additional fitting parameters
    copyfile(mrc_file, f"./dnaFit/{prefix}.mrc")

    # NOTE changing design and sequence file to copy in folder dnaFit to ensure consistency
    mrc_file = Path(f"./dnaFit/{prefix}.mrc").resolve()
    cad_file = Path(f"./dnaFit/{prefix}.json").resolve()
    seq_file = Path(f"./dnaFit/{prefix}.seq").resolve()
    top = Path(f"./dnaFit/{prefix}.psf").resolve()
    conf = Path(f"./dnaFit/{prefix}-undocked.pdb").resolve()
    exb = Path(f"./dnaFit/{prefix}.exb").resolve()

    home_directory = os.getcwd()
    try:
        os.chdir("dnaFit")
        logger.debug(f"changing directory to: {os.getcwd()}")

        # NOTE: creating Cascade object triggers external docking prompt
        cascade = Cascade(conf=conf, top=top, mrc=mrc_file,
                          exb=exb, recenter=True, is_docked=False)
        # NOTE: fit is moved back to original mrc position, recentering invisible to user
        dnaFit = cascade.run_cascaded_fitting(
            base_time_steps=timesteps, resolution=resolution)
        dnaFit.write_linkage(cad_file, seq_file)
        dnaFit.write_output(dest=home_directory,
                            write_mmCif=True, crop_mrc=True)

    finally:
        os.chdir(home_directory)
        logger.debug(f"changing directory to: {os.getcwd()}")


@cli.command()
@click.argument('cadnano', type=click.Path(exists=True))
@click.argument('sequence', type=click.Path(exists=True))
@click.argument('mrc', type=click.Path(exists=True))
@click.argument('top', type=click.Path(exists=True))
@click.argument('conf', type=click.Path(exists=True))
@click.option('-g', '--gpu', type=int, default=0, help='GPU used for simulation', show_default=True)
@click.option('-o', '--output-prefix', 'prefix', type=str, default=None,
              help="short design name, default to json name")
@click.option('--timesteps', type=int, default=12000,
              help='timesteps per cascade (multiple of 12)')
@click.option('--resolution', type=float, default=10.0,
              help='mrc map resolution in Angstrom')
@click.option('--pdb-docked', is_flag=True,
              help='add if pdb has already been docked to the mrc data')
def fit(cadnano, sequence, mrc, top, conf, gpu, prefix, timesteps, resolution, pdb_docked):
    """Cascaded mrDNA-driven MD flexible fitting to MRC cryo data, creates dnaFit folder

        CADNANO is the name of the design file [.json]\n
        SEQUENCE is the scaffold strand sequence file [.txt, .seq]\n
        MRC is the name of the cryo EM volumetric data file [.mrc]\n
        TOP is the name of the namd topology file [.top]\n
        CONF is the name of the namd configuration file [.pdb, .coor]\n
    """
    # NOTE: click will drop python2 support soon and actually return a Path
    cad_file = Path(cadnano).resolve()
    seq_file = Path(sequence).resolve()
    mrc_file = Path(mrc).resolve()
    top = Path(top).resolve()
    conf = Path(conf).resolve()
    prefix = cad_file.stem if prefix is None else prefix

    # create duplicates of input files in dnaFit folder
    Path("dnaFit").mkdir(parents=True, exist_ok=True)
    copyfile(cadnano, f"./dnaFit/{prefix}.json")
    copyfile(sequence, f"./dnaFit/{prefix}.seq")
    copyfile(mrc, f"./dnaFit/{prefix}.mrc")
    pdb_tag = "" if pdb_docked else "-undocked"
    copyfile(conf, f"./dnaFit/{prefix}{pdb_tag}.pdb")
    copyfile(top, f"./dnaFit/{prefix}.psf")

    home_directory = os.getcwd()
    try:
        os.chdir("dnaFit")
        logger.debug(f"changing directory to: {os.getcwd()}")
        mrc_file = Path(f"./{prefix}.mrc").resolve()
        cad_file = Path(f"./{prefix}.json").resolve()
        seq_file = Path(f"./{prefix}.seq").resolve()
        top_file = Path(f"./{prefix}.psf").resolve()
        con_file = Path(f"./{prefix}{pdb_tag}.pdb").resolve()
        exb_file = Path(f"./{prefix}.exb").resolve()
        if not Path("./charmm36.nbfix").exists():
            copytree(get_resource("charmm36.nbfix"), "./charmm36.nbfix")

        # NOTE: creating Cascade object triggers external docking prompt
        cascade = Cascade(conf=con_file, top=top_file, mrc=mrc_file,
                          exb=exb_file, recenter=True, is_docked=False)
        dnaFit = cascade.run_cascaded_fitting(
            base_time_steps=timesteps, resolution=resolution)
        dnaFit.write_linkage(cad_file, seq_file)
        dnaFit.write_output(dest=home_directory,
                            write_mmCif=True, crop_mrc=True)

    finally:
        os.chdir(home_directory)
        logger.debug(f"changing directory to: {os.getcwd()}")


@cli.command()
@click.argument('cadnano', type=click.Path(exists=True))
@click.argument('sequence', type=click.Path(exists=True))
@click.argument('mrc', type=click.Path(exists=True))
@click.argument('top', type=click.Path(exists=True))
@click.argument('conf', type=click.Path(exists=True))
@click.option('--enrgMD_server', is_flag=True,
              help='add if pdb has already been docked to the mrc data')
def link(cadnano, sequence, mrc, top, conf, enrgMD_server):
    """ links structural information of the cadnano designfile[design.json] to
         fitted atomic model[design.psf, design.dcd].
        * linkage information ist stored in human readable csv format

        CADNANO is the name of the design file [.json]\n
        SEQUENCE is the scaffold strand sequence file [.txt, .seq]\n
        MRC is the name of the cryo EM volumetric data file [.mrc]\n
        TOP is the name of the namd topology file [.top]\n
        CONF is the name of the namd configuration file [.pdb, .coor]\n
    """
    # NOTE: click will drop python2 support soon and actually return a Path
    mrc = Path(mrc).resolve()
    top = Path(top).resolve()
    conf = Path(conf).resolve()
    cadnano = Path(cadnano).resolve()
    sequence = Path(sequence).resolve()

    dnaFit = AtomicModelFit(conf=conf, top=top, mrc=mrc,
                            generated_with_mrdna=enrgMD_server)
    dnaFit.write_linkage(json=cadnano, seq=sequence)


@cli.command()
@click.argument('mrc', type=click.Path(exists=True))
@click.argument('top', type=click.Path(exists=True))
@click.argument('conf', type=click.Path(exists=True))
@click.option('--enrgMD_server', is_flag=True,
              help='add if pdb has already been docked to the mrc data')
def mask(mrc, top, conf, enrgMD_server):
    """ links structural information of the cadnano designfile[design.json] to
         fitted atomic model[design.psf, design.dcd].
         used to mask mrc map to fit

        MRC is the name of the cryo EM volumetric data file [.mrc]\n
        TOP is the name of the namd topology file [.top]\n
        CONF is the name of the namd configuration file [.pdb, .coor]\n
    """
    # NOTE: click will drop python2 support soon and actually return a Path
    mrc = Path(mrc).resolve()
    top = Path(top).resolve()
    conf = Path(conf).resolve()

    dnaFit = AtomicModelFit(conf=conf, top=top, mrc=mrc,
                            generated_with_mrdna=enrgMD_server)
    u = dnaFit.get_universe()
    mrc_masked = mrc.with_name(f"{mrc.stem}-masked.mrc")
    write_mrc_from_atoms(path=mrc, atoms=u.atoms,
                         path_out=mrc_masked, context=10., cut_box=True)


@cli.command()
@click.argument('pdb', type=click.Path(exists=True))
@click.option('--remove-H', is_flag=True,
              help='remove hydrogen atoms')
def pdb2CIF(pdb, remove_H):
    """ generate mmCIF from namd pdb

        PDB is the name of the namd configuration file [.pdb]\n
    """
    # NOTE: click will drop python2 support soon and actually return a Path
    pdb = Path(pdb).resolve()

    structure = Structure(path=pdb, remove_H=remove_H)
    structure.parse_pdb()
    # TODO: -low- ask for additional info (name, author, etc)
    output_name = pdb.with_suffix(".cif")
    structure.write_cif(output_name)
