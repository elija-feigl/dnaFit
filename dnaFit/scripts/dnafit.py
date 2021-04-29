#! /usr/bin/env python
# -*- coding: utf-8 -*-

import logging
import os

from pathlib import Path
from shutil import copyfile, copytree

import click
from dnaFit import get_resource
from dnaFit.core.mrDna import run_mrDNA, prep_cascaded_fitting, recenter_conf
from dnaFit.data.mrc import write_mrc_from_atoms, recenter_mrc
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
def mrDna(cadnano, mrc, sequence, gpu, prefix, timesteps, resolution, multidomain):
    """ mrDNA simulation of CADNANO design file followed by prep of cascaded
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

    prep_cascaded_fitting(prefix, cad_file, seq_file, mrc_file)

    # TODO: -low- parse & set additional fitting parameters
    # TODO: move to separate object
    copyfile(mrc_file, f"./dnaFit/{prefix}.mrc")

    logger.info("Config file is moved to center of mass with mrc map but still \
        has to be rotated before fitting. execute vmd_info for additional info")

# TODO: enrgMD only setup


@cli.command()
@click.argument('mrc', type=click.Path(exists=True))
@click.argument('top', type=click.Path(exists=True))
@click.argument('conf', type=click.Path(exists=True))
def center_on_map(mrc, top, conf):
    """ recenter atomic model on mrc cryo map

        MRC is the name of the cryo EM volumetric data file [.mrc]\n
        TOP is the name of the namd topology file [.top]\n
        CONF is the name of the namd configuration file, docked to map (VMD)  [.pdb, .coor]\n
    """
    # NOTE: click will drop python2 support soon and actually return a Path
    mrc = Path(mrc).resolve()
    top = Path(top).resolve()
    conf = Path(conf).resolve()
    conf_docked = conf.with_name(f"{conf.stem}-docked.pdb")
    copyfile(conf, conf_docked)

    mrc_shift = recenter_mrc(mrc, apply=False)
    recenter_conf(top=top, conf=conf_docked, to_position=mrc_shift)


@ cli.command()
def vmd_info():
    """ Print VMD command for rotation of pdb around center of mass
    """
    logger.info(get_resource("vmd_rot.txt").read_text())


@ cli.command()
@ click.argument('cadnano', type=click.Path(exists=True))
@ click.argument('sequence', type=click.Path(exists=True))
@ click.argument('mrc', type=click.Path(exists=True))
@ click.argument('top', type=click.Path(exists=True))
@ click.argument('conf', type=click.Path(exists=True))
@ click.argument('exb', type=click.Path(exists=True))
@ click.option('-g', '--gpu', type=int, default=0, help='GPU used for simulation', show_default=True)
@ click.option('-o', '--output-prefix', 'prefix', type=str, default=None,
               help="short design name, default to json name")
@ click.option('--timesteps', type=int, default=12000,
               help='timesteps per cascade (multiple of 12)')
@ click.option('--resolution', type=float, default=10.0,
               help='mrc map resolution in Angstrom')
@ click.option('--SR-fitting', is_flag=True,
               help='retain Sr bonds troughout fitting cascade.')
def fit(cadnano, sequence, mrc, top, conf, exb, gpu, prefix, timesteps, resolution, sr_fitting):
    """Cascaded mrDNA-driven MD flexible fitting to MRC cryo data, creates dnaFit folder

        CADNANO is the name of the design file [.json]\n
        SEQUENCE is the scaffold strand sequence file [.txt, .seq]\n
        MRC is the name of the cryo EM volumetric data file [.mrc]\n
        TOP is the name of the namd topology file [.top]\n
        CONF is the name of the namd configuration file, docked to map (VMD)  [.pdb, .coor]\n
        EXB is the name of the enrgMD extrabond file (expect mrDNA > march 2021) [.exb]
    """
    # NOTE: click will drop python2 support soon and actually return a Path
    cad_file = Path(cadnano).resolve()
    seq_file = Path(sequence).resolve()
    mrc_file = Path(mrc).resolve()
    top = Path(top).resolve()
    conf = Path(conf).resolve()
    exb = Path(exb).resolve()
    prefix = cad_file.stem if prefix is None else prefix

    # create duplicates of input files in dnaFit folder
    Path("dnaFit").mkdir(parents=True, exist_ok=True)
    copyfile(cadnano, f"./dnaFit/{prefix}.json")
    copyfile(sequence, f"./dnaFit/{prefix}.seq")
    copyfile(mrc, f"./dnaFit/{prefix}.mrc")
    copyfile(conf, f"./dnaFit/{prefix}.pdb")
    copyfile(top, f"./dnaFit/{prefix}.psf")
    copyfile(exb, f"./dnaFit/{prefix}.exb")

    home_directory = os.getcwd()
    try:
        os.chdir("dnaFit")
        logger.debug(f"changing directory to: {os.getcwd()}")
        mrc_file = Path(f"./{prefix}.mrc").resolve()
        cad_file = Path(f"./{prefix}.json").resolve()
        seq_file = Path(f"./{prefix}.seq").resolve()
        top_file = Path(f"./{prefix}.psf").resolve()
        con_file = Path(f"./{prefix}.pdb").resolve()
        exb_file = Path(f"./{prefix}.exb").resolve()
        if not Path("./charmm36.nbfix").exists():
            copytree(get_resource("charmm36.nbfix"), "./charmm36.nbfix")

        cascade = Cascade(conf=con_file, top=top_file,
                          mrc=mrc_file, exb=exb_file)
        dnaFit = cascade.run_cascaded_fitting(
            base_time_steps=timesteps, resolution=resolution, is_SR=sr_fitting)
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
def pdb2CIF(pdb, remove_h):
    """ generate mmCIF from namd pdb

        PDB is the name of the namd configuration file [.pdb]\n
    """
    # NOTE: click will drop python2 support soon and actually return a Path
    pdb = Path(pdb).resolve()

    structure = Structure(path=pdb, remove_H=remove_h)
    structure.parse_pdb()
    # TODO: -low- ask for additional info (name, author, etc)
    output_name = pdb.with_suffix(".cif")
    structure.write_cif(output_name)
