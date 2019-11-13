#!/usr/bin/env python
import MDAnalysis as mda
import argparse
import os
import sys


from pathlib import Path

from project import Project
from utils import ignored
from linker import get_linkage
from segmentation import categorise, mrc_segment


def mask_minimal_box(u, project):
    path_in = project.input / "{}.mrc".format(project.name)
    path_out = project.output / "{}-masked.mrc".format(project.name)
    mrc_segment(
        atoms=u.atoms,
        path_in=path_in,
        path_out=path_out,
        context=project.context,
    )


def check_abort() -> None:
    yes = {"yes", "y", "ye"}
    no = {"no", "n"}
    choice = ""
    print("abort without segmentation?")
    while choice not in yes.union(no):
        choice = input().lower()
        if choice in yes:
            sys.exit(0)
        elif choice in no:
            return
        else:
            print("Please respond with 'yes' or 'no'")


def get_description() -> str:
    return """cut subset from map according to atoms belongign to sepcific
              motif. Also produces minimal box map. can also segment halfmaps
              and evaluate local-resolution per residue -> dict and pdb
              """


def proc_input() -> Project:
    parser = argparse.ArgumentParser(
        description=get_description(),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--folder",
                        help="input folder",
                        type=str,
                        default="./",
                        )
    parser.add_argument("--name",
                        help="name of design and files",
                        type=str,
                        required=True,
                        default=argparse.SUPPRESS,
                        )
    parser.add_argument("--context",
                        help="context aroud atom in Angstrom",
                        type=int,
                        default=5,
                        )
    parser.add_argument("--range",
                        help="number of addit. basepairs in helix of motif",
                        type=int,
                        default=10,
                        )
    parser.add_argument("--halfmap",
                        help="also segment halfmaps",
                        action="store_true"
                        )
    parser.add_argument("--star",
                        help="create starfile",
                        action="store_true"
                        )
    parser.add_argument("--relink",
                        help="force relink fit",
                        action="store_true"
                        )
    args = parser.parse_args()
    project = Project(input=Path(args.folder),
                      output=Path(args.folder) / "analysis",
                      name=args.name,
                      context=args.context,
                      range=args.range,
                      halfmap=args.halfmap,
                      star=args.star,
                      relink=args.relink,
                      )
    return project


def main():
    H1 = "_unfil_half1"
    H2 = "_unfil_half2"

    project = proc_input()
    link = get_linkage(project)
    link.u.trajectory[-1]

    print("mask minimal box")
    mask_minimal_box(link.u, project)

    if project.halfmap:
        print("segmenting halfmaps")

    check_abort()
    motifs = categorise(link=link, project=project)
    for motif_name, motif in motifs.items():
        path_motif = project.output / motif_name
        print("output to ", path_motif)
        with ignored(FileExistsError):
            os.mkdir(path_motif)

        for subset in motif:
            base_selection, key, typ = subset
            atoms_select = mda.AtomGroup([], link.u)
            for resindex in base_selection:
                atoms_select += link.u.residues[resindex].atoms

            if project.halfmap:
                specs = {"-segment": "", H1: "h1-", H2: "h2-"}
            else:
                specs = {"-segment": ""}
            for halfmap_inp, halfmap_out in specs.items():
                in_suffix = "{}{}.mrc".format(project.name,
                                              halfmap_inp,
                                              )
                out_suffix = "{}__{}{}_{}.mrc".format(project.name,
                                                      halfmap_out,
                                                      typ,
                                                      key,
                                                      )
                path_in = project.input / in_suffix
                path_out = path_motif / out_suffix
                mrc_segment(
                    atoms=atoms_select,
                    path_in=path_in,
                    path_out=path_out,
                    context=project.context,
                    star=project.star,
                )


if __name__ == "__main__":
    main()
