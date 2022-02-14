#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
""" cascaded mrDNA-driven MD flexible fitting script module:
    dnaFit commands:
        main pipeline:
        1. mrdna:   mrDNA structure prediction with custom settings for cadnano file.
                    (includes creation of prep folder for fitting)
        2. vmd_info:    print info for rigid body docking with VMD
        3. fit:     shrink wrap fitting if of rigid body docked mrDNA prediction
                    includes mask, pdb2cif, and link

        additional exposed commands:
        * center_on_map:    write new PDB with COM of the configuration at .mrc center
        * link:     create a linkage file for cadnano-design and atomic model
                    (required by FitViewer)
        * mask:     mask a .mrc map using an atomic model
        * pdb2cif:  create .cif from .pdb file
 """
import logging
import os
from pathlib import Path
from shutil import copyfile
from shutil import copytree

import click
from dnaFit import get_resource
from dnaFit import get_version
from dnaFit.core.predict import prep_cascaded_fitting
from dnaFit.core.predict import recenter_conf
from dnaFit.core.predict import run_mrdna
from dnaFit.core.utils import _check_path
from dnaFit.data.mrc import get_mrc_center
from dnaFit.data.mrc import write_mrc_from_atoms
from dnaFit.fit.atomic_model_fit import AtomicModelFit
from dnaFit.fit.cascade import Cascade
from pdb2cif.pdb.structure import Structure

logger = logging.getLogger(__name__)


def print_version(ctx, _, value):
    """click print version."""
    if not value or ctx.resilient_parsing:
        return
    click.echo(get_version())
    ctx.exit()


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "-v",
    "--version",
    is_flag=True,
    help="Show __version__ and exit.",
    callback=print_version,
    expose_value=False,
    is_eager=True,
)
def cli():
    """Cascaded mrDNA-driven MD flexible fitting script module:

    \b
    dnaFit commands:
        main pipeline:
        1. mrdna:   mrDNA structure prediction with custom settings for cadnano file.
                    (includes creation of prep folder for fitting)
        2. vmd_info:    print info for rigid body docking with VMD
        3. fit:     shrink wrap fitting if of rigid body docked mrDNA prediction
                    includes mask, pdb2cif, and link
    """


@cli.command()
@click.argument("cadnano", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("sequence", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("mrc", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.option(
    "-g",
    "--gpu",
    type=int,
    default=0,
    show_default=True,
    help="NOT SUPPORTED. GPU used for simulation.",
)
@click.option(
    "-o",
    "--output-prefix",
    "prefix",
    type=str,
    default=None,
    help="short design name, default to json name",
)
@click.option(
    "--coarse-steps",
    "coarse_steps",
    type=float,
    default=5e7,
    show_default=True,
    help='multidomain only. see "mrdna --help" --coarse-steps',
)
@click.option(
    "--bond-cutoff",
    "bond_cutoff",
    type=int,
    default=300,
    show_default=True,
    help='multidomain only. see "mrdna --help" --coarse-bond-cutoff',
)
@click.option("--multidomain", is_flag=True, help="use multidomain structures settings.")
@click.option("--no-prep", is_flag=True, help="do not perform fitting prep.")
def mrdna(cadnano, mrc, sequence, gpu, prefix, multidomain, coarse_steps, bond_cutoff, no_prep):
    """\b
    MrDNA simulation of CADNANO design file with custom settings.
        followed by preperation of files for "dnaFit fit"
    \b
    Note1:  includes centering of model and masking of map
    Note2:  map and model are not rigid body docked "dnaFit vmd_info"
    \b
    CADNANO is the name of the design file [.json]
    SEQUENCE is the scaffold strand sequence file [.txt, .seq]
    MRC is the name of the cryo EM volumetric data file [.mrc]
    """
    file_types = {
        cadnano: [".json"],
        mrc: [".mrc"],
        sequence: [".txt", ".seq"],
    }
    for file, types in file_types.items():
        _check_path(file, types)
    prefix = cadnano.stem if prefix is None else prefix

    run_mrdna(
        cad_file=cadnano,
        seq_file=sequence,
        prefix=prefix,
        gpu=gpu,
        multidomain=multidomain,
        coarse_steps=coarse_steps,
        bond_cutoff=bond_cutoff,
    )
    if not no_prep:
        prep_cascaded_fitting(prefix, cadnano, sequence, mrc, multidomain=multidomain)

    logger.info(
        "Config file is moved to center of mass with mrc map but still \
        has to be rotated before fitting. execute vmd_info for additional info"
    )


@cli.command()
@click.argument("cadnano", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("sequence", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("mrc", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.option(
    "-o",
    "--output-prefix",
    "prefix",
    type=str,
    default=None,
    help="short design name, default to json name",
)
@click.option(
    "--multidomain", is_flag=True, help="set true multidomain structures settings used for mrDNA."
)
def prep(cadnano, mrc, sequence, prefix, multidomain):
    """\b
    Prepare mrDNA results for fitting "dnaFit fit".
    \b
    Note1:  includes centering of model and masking of map
    Note2:  map and model are not rigid body docked "dnaFit vmd_info"
    \b
    CADNANO is the name of the design file [.json]
    SEQUENCE is the scaffold strand sequence file [.txt, .seq]
    MRC is the name of the cryo EM volumetric data file [.mrc]
    """
    file_types = {
        cadnano: [".json"],
        mrc: [".mrc"],
        sequence: [".txt", ".seq"],
    }
    for file, types in file_types.items():
        _check_path(file, types)
    prefix = cadnano.stem if prefix is None else prefix
    prep_cascaded_fitting(prefix, cadnano, sequence, mrc, multidomain=multidomain)


@cli.command()
@click.argument("mrc", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("top", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("conf", type=click.Path(exists=True, resolve_path=True, path_type=Path))
def center_on_map(mrc, top, conf):
    """\b
    Recenter atomic model center-of-mass on mrc cryo map center and write new .pdb.
    \b
    MRC is the name of the cryo EM volumetric data file [.mrc]
    TOP is the name of the namd topology file [.psf]
    CONF is the name of the namd configuration, docked to the map [.pdb, .coor]
    """
    file_types = {
        mrc: [".mrc"],
        top: [".psf"],
        conf: [".pdb", ".coor"],
    }
    for file, types in file_types.items():
        _check_path(file, types)
    conf_docked = conf.with_name(f"{conf.stem}-docked.pdb")
    copyfile(conf, conf_docked)

    mrc_shift = get_mrc_center(mrc)
    recenter_conf(top=top, conf=conf_docked, to_position=mrc_shift)


@cli.command()
def vmd_info():
    """\b
    Print VMD command for rotation of pdb around center of mass.
    \b
    Note:   The current pipeline does not include an automated rigid-body-docking step.
            Docking of the model has to be performed manually. Using VMD is recommended.
    """
    logger.info(get_resource("vmd_rot.txt").read_text())


@cli.command()
@click.argument("cadnano", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("sequence", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("mrc", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("top", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("conf", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("exb", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.option(
    "-o",
    "--output-prefix",
    "prefix",
    type=str,
    default=None,
    help="short design name, default to json name",
)
@click.option(
    "--timesteps",
    type=int,
    default=12000,
    show_default=True,
    help="timesteps per cascade (multiple of 12)",
)
@click.option(
    "--resolution",
    type=float,
    default=10.0,
    show_default=True,
    help="mrc map resolution in Angstrom",
)
@click.option(
    "--SR-fitting",
    "sr_fitting",
    is_flag=True,
    help="retain SR-elastic-network troughout fitting cascade.",
)
@click.option(
    "--exclude-ss",
    "exclude_ss",
    is_flag=False,
    help="dont include single stranded DNA from flexible fitting.",
)
@click.option(
    "--grid-pdb",
    "grid_pdb",
    type=click.Path(exists=True),
    default=None,
    help="use custom grid.pdb for segment exclusion.",
)
def fit(
    cadnano,
    sequence,
    mrc,
    top,
    conf,
    exb,
    prefix,
    timesteps,
    resolution,
    sr_fitting,
    exclude_ss,
    grid_pdb,
):
    """\b
    Cascaded mrDNA-driven MD flexible fitting (shrink wrap fitting) to MRC cryo data.
        creates dnaFit folder. This folder will contain final results prefix-last.cif
    \b
    Note1:  run after mrDNA ("dnaFit mrdna") and manual rigid-body-docking ("dnaFit vmd_info")
    Note2:  SR-fitting (Short Range). Recommended for maps with resolution >> 10 Angstrom
    \b
    CADNANO is the name of the design file [.json]
    SEQUENCE is the scaffold strand sequence file [.txt, .seq]
    MRC is the name of the cryo EM volumetric data file [.mrc]
    TOP is the name of the namd topology file [.psf]
    CONF is the name of the namd configuration file, docked to map (VMD)  [.pdb, .coor]
    EXB is the name of the enrgMD extrabond file (expect mrDNA > march 2021) [.exb]
    """
    file_types = {
        cadnano: [".json"],
        sequence: [".seq", ".txt"],
        mrc: [".mrc"],
        top: [".psf"],
        conf: [".pdb", ".coor"],
        exb: [".exb"],
    }
    for file, types in file_types.items():
        _check_path(file, types)
    prefix = cadnano.stem if prefix is None else prefix

    # TODO-low: gpu support
    logger.info("GPU currently not supported for fitting. using CPU only.")

    # create duplicates of input files in dnaFit folder
    Path("dnaFit").mkdir(parents=True, exist_ok=True)
    for file, types in file_types.items():
        copyfile(file, f"./dnaFit/{prefix}{types[0]}")

    home_directory = Path.cwd()
    try:
        os.chdir("dnaFit")
        logger.debug("changing directory to: %s", Path.cwd())
        # use duplicates instead of input
        mrc = Path(f"./{prefix}.mrc").resolve()
        cadnano = Path(f"./{prefix}.json").resolve()
        sequence = Path(f"./{prefix}.seq").resolve()
        top = Path(f"./{prefix}.psf").resolve()
        conf = Path(f"./{prefix}.pdb").resolve()
        exb = Path(f"./{prefix}.exb").resolve()
        if not Path("./charmm36.nbfix").exists():
            copytree(get_resource("charmm36.nbfix"), "./charmm36.nbfix")

        cascade = Cascade(
            conf=conf,
            top=top,
            mrc=mrc,
            exb=exb,
            json=cadnano,
            seq=sequence,
            custom_grid_pdb=grid_pdb,
        )
        model = cascade.run_cascaded_fitting(
            time_steps=timesteps,
            resolution=resolution,
            is_sr=sr_fitting,
            exclude_ss=exclude_ss,
        )
        model.write_linkage(cadnano, sequence)
        model.write_output(dest=Path(home_directory), write_mmcif=True, mask_mrc=True)
    finally:
        os.chdir(home_directory)
        logger.debug("changing directory to: %s", Path.cwd())


@cli.command()
@click.argument("cadnano", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("sequence", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("top", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("conf", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.option(
    "--enrgmd-server",
    "enrgmd_server",
    is_flag=True,
    help="add if initial .pdb was generated with enrgMD web-server.",
)
def link(cadnano, sequence, top, conf, enrgmd_server):
    """\b
    Links structural information of the CADNANO designfile to
        fitted atomic model files TOP, CONF.
    \b
    Note: linkage information ist stored in human readable csv format
    \b
    CADNANO is the name of the design file [.json]
    SEQUENCE is the scaffold strand sequence file [.txt, .seq]
    TOP is the name of the namd topology file [.psf]
    CONF is the name of the namd configuration file [.pdb, .coor]
    """
    file_types = {
        cadnano: [".json"],
        sequence: [".txt", ".seq"],
        top: [".psf"],
        conf: [".pdb", ".coor"],
    }
    for file, types in file_types.items():
        _check_path(file, types)

    model = AtomicModelFit(conf=conf, top=top, mrc=None, generated_with_mrdna=(not enrgmd_server))
    model.write_linkage(json=cadnano, seq=sequence)


@cli.command()
@click.argument("mrc", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("top", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.argument("conf", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.option(
    "-e",
    "--enrgmd-server",
    "enrgmd_server",
    is_flag=True,
    help="add if pdb has been generated with enrgMD server",
)
@click.option(
    "--dont-cut-box", "no_cut_box", is_flag=True, help="do not reduce grid to minimum box."
)
@click.option("--keep-full", "keep_full", is_flag=True, help="retain all data within minimum box.")
def mask(mrc, top, conf, enrgmd_server, no_cut_box, keep_full):
    """\b
    Mask mrc map to fitted atomic model.
        Used to:
            * make minimum box by keeping full data (for manual fit prep)
            * create mask for atomic model
    \b
    MRC is the name of the cryo EM volumetric data file [.mrc]
    TOP is the name of the namd topology file [.top]
    CONF is the name of the namd configuration file [.pdb, .coor]
    """
    file_types = {
        mrc: [".mrc"],
        top: [".psf"],
        conf: [".pdb", ".coor"],
    }
    for file, types in file_types.items():
        _check_path(file, types)

    model = AtomicModelFit(conf=conf, top=top, mrc=mrc, generated_with_mrdna=enrgmd_server)
    universe = model.get_universe()
    mrc_masked = mrc.with_name(f"{mrc.stem}-masked.mrc")
    write_mrc_from_atoms(
        path=mrc,
        atoms=universe.atoms,
        path_out=mrc_masked,
        context=40.0,
        cut_box=no_cut_box,
        keep_data=keep_full,
    )


@cli.command()
@click.argument("pdb", type=click.Path(exists=True, resolve_path=True, path_type=Path))
@click.option("--remove-H", "remove_h", is_flag=True, help="remove hydrogen atoms")
@click.option("--snupi", "is_snupi", is_flag=True, help="convert from SNUPI pdb")
@click.option(
    "--flip-fields",
    "flip_fields",
    is_flag=True,
    help="flip the values of occupancy and temperature",
)
def pdb2cif(pdb, remove_h, is_snupi, flip_fields):
    """\b
    Generate atomic model in mmCIF format from namd PDB
    \b
    PDB is the name of the namd configuration file [.pdb]
    """
    _check_path(pdb, [".pdb"])
    structure = Structure(path=pdb, remove_H=remove_h, is_snupi=is_snupi, flip_fields=flip_fields)
    structure.parse_pdb()
    output_name = pdb.with_suffix(".cif")
    structure.write_cif(output_name)
