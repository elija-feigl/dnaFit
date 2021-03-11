import os
import contextlib
import argparse
from pathlib import Path

from ..core.project import Project
from ..link.linker_project import Linker
from ..link.elastic_network import ElaticNetwortModifier
from ..version import get_version, get_authors


__descr__ = """
    links structural information of the cadnano designfile
    [design.json] to fitted atomic model [design.psf, design.dcd].
    stores dictionaries as pickles containg mapping for motifs,
    residue-id, lattice position and base-pairing.
    allows modification of Elastic Network of enrgMD simulation
"""


@contextlib.contextmanager
def ignored(*exceptions):
    try:
        yield
    except exceptions:
        pass


def proc_input() -> Project:
    def get_description() -> str:
        return "{}\n {}\n {}".format(__descr__, get_version, get_authors)
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
    print("linkage output to {}".format(project.output))
    linkage.dump_linkage(project=project)

    if project.ENmodify:
        print("modifying extrabonds")
        en = ElaticNetwortModifier(linker)
        en.write_en()


if __name__ == "__main__":
    main()
