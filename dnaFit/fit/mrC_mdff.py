import attr
from pathlib import Path
from typing import Any

from ..core.utils import _get_executable

""" DESCR:
    mrDNA driven cascade fitting simulation class.
"""


def external_docking_loop(prefix):
    # TODO: implement automatic alignmant or python UI based alignment
    # BREAK: reorient helix and fit still external
    conf = Path("./{}.pdb".format(prefix))
    while True:
        input(
            "You need to manualy align mrc and pdb using VMD. save as {} and press ENTER".format(conf))
        if conf.is_file():
            return conf


@attr.s
class Cascade(object):
    """ Cascade simulation class
    """
    conf: Path = attr.ib()
    top: Path = attr.ib()
    exb: Path = attr.ib()
    mrc: Path = attr.ib()
    recenter: bool = attr.ib()
    is_docked: bool = attr.ib()

    def __attrs_post_init__(self) -> None:
        self.prefix: str = self.top.stem

        self.mrc_shift = self._recenter_mrc()
        self._recenter_conf()

        if not self.is_docked:
            self.conf = external_docking_loop(self.prefix)

        self._split_exb_file()

        self.charmrun = _get_executable("charmrun")
        self.namd2 = _get_executable("namd2")
        self.vmd = _get_executable("vmd")

    def _recenter_mrc(self) -> Any:
        # TODO: recenter map
        #       NOTE: save translation:
        raise NotImplementedError

    def _recenter_conf(self) -> None:
        # TODO: centered conf at origin (0,0,0 simplifies vmd rotation)
        #       NOTE: edit original file and discard translation
        raise NotImplementedError

    def _split_exb_file(self) -> None:
        # TODO: split .exb file
        # NOTE: assuming modifed mrDNA > march 2021, (annotated, sorted .exb files)
        #   NOTE: splitting here allows starting from an old mrDNA or energMD setup. as long as .exb is orderer (extra script)
        ...
        raise NotImplementedError

    def run_cascaded_fitting(self, time_step: int, resolution: float):
        # TODO: allow additional parameter changes
        # NOTE: assuming we can put all files in the current folder!

        # TODO: either prep & call .sh script or call namd directly
        ...
        # TODO: setup fitting script and call it
        #       NOTE: pure enrgMD not required if mrDNA default
        # TODO: reset pdb and map loaction to orginal mrc using transl
        ...
        # TODO: returns dnaFit object
        return AtomicModelFit()


@attr.s
class AtomicModelFit(object):

    def write_linkage(self, cad_file: Path, seq_file: Path):
        # TODO: write persistent and human readable linkage (Fbp, FidDid)
        raise NotImplementedError

    def write_output(self, write_mmCif=True, crop_mrc=True):
        # TODO: save final pdb/mmcif together with mrc and masked mrc
        raise NotImplementedError
