#!/usr/bin/env python3#
import os
import contextlib
import argparse

from pathlib import Path

from project import Project
from linker import Linker
from elastic_network import ElaticNetwortModifier


@contextlib.contextmanager
def ignored(*exceptions):
    try:
        yield
    except exceptions:
        pass


def get_description() -> str:
    return """links structural information of the cadnano designfile
              [design.json] to fitted atomic model [design.psf, design.dcd].
              stores dictionaries as pickles containg mapping for motifs,
              residue-id, lattice position and base-pairing."""


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
    parser.add_argument("--ENmodify",
                        help="reset nomenclature enrgMD",
                        action="store_true"
                        )
    parser.add_argument("--EN",
                        help="""bool-string:
                                long,strand,Hbond,xstack,nick,co,ss,dhdrl""",
                        type=str,
                        default="11111110",
                        )
    args = parser.parse_args()
    project = Project(input=Path(args.folder),
                      output=Path(args.folder) / "analysis",
                      name=args.name,
                      ENmodify=args.ENmodify,
                      EN=args.EN,
                      )
    with ignored(FileExistsError):
        os.mkdir(project.output)
    return project


def main():
    project = proc_input()

    linker = Linker(project)
    linkage = linker.create_linkage()

    if not project.ENmodify:
        print("linkage output to {}".format(project.output))
        linkage.dump_linkage(project=project)

    else:
        print("modifying extrabonds")
        en = ElaticNetwortModifier(linker)
        en.write_en()


if __name__ == "__main__":
    main()
