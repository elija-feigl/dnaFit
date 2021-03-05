import attr
import os
from pathlib import Path
#from typing import Dict, Tuple, Optional, List, Set

from ..core.utils import _get_executable

""" DESCR:
    mrDNA driven cascade fitting simulation class.
"""


@attr.s
class mrCascade(object):
    """ mrCascade simulation class
    """
    conf: Path = attr.ib()
    top: Path = attr.ib()
    json: Path = attr.ib()
    mrc: Path = attr.ib()
    seq: Path = attr.ib()

    prefix: str = attr.ib()

    def __attrs_post_init__(self) -> None:
        # TODO: recenter pdb and map
        #       NOTE: centered at origin (0,0,0 simplifies vmd rotation)
        #       NOTE: save translation:
        ...

        self.charmrun = _get_executable("charmrun")
        self.namd2 = _get_executable("namd2")
        self.vmd = _get_executable("vmd")

    def run_cascaded_fitting(self, time_step: int, resolution: float):
        # TODO: allow additional parameter changes
        home_directory = os.getcwd()
        try:
            os.chdir("dnaFit")

            # TODO: either prep & call .sh script or call namd directly
            ...
            # TODO: setup fitting script and call it
            #       NOTE: pure enrgMD not required if mrDNA default
            # TODO: reset pdb and map loaction to orginal mrc using transl
            ...

        finally:
            os.chdir(home_directory)

    def create_linkage(self):
        # TODO: create persistent linkage, masked mrc
        # NOTE: convert to cif better in dnaFit (package independence)
        raise NotImplementedError
