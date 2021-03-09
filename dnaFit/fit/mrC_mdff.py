import attr
from pathlib import Path
from typing import Any

import numpy as np
import mrcfile
import MDAnalysis as mda

from ..core.utils import _get_executable

""" DESCR:
    mrDNA driven cascade fitting simulation class.
"""


def external_docking_loop(prefix):
    # TODO -low-: implement automatic alignmant or python UI based alignment
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

        _ = self._recenter_mrc()
        self._recenter_conf()

        if not self.is_docked:
            self.conf = external_docking_loop(self.prefix)

        self._split_exb_file()

        self.charmrun = _get_executable("charmrun")
        self.namd2 = _get_executable("namd2")
        self.vmd = _get_executable("vmd")

    def _recenter_mrc(self, reverse=False) -> Any:
        with mrcfile.open(self.mrc, mode='r+') as mrc:
            c = np.array(mrc.header["cella"])
            cell = np.array([c["x"], c["y"], c["z"]])
        if reverse:
            shift = cell * 0.5 if reverse else cell * - 0.5
        mrc.header["origin"] = tuple(shift)
        return shift

    def _recenter_conf(self, to_position=np.array([0.0, 0.0, 0.0])) -> None:
        u = mda.Universe(str(self.top), str(self.conf))
        translation = to_position - u.atoms.center_of_geometry
        u.atoms.translate(translation)
        u.atoms.write(str(self.conf))

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

        mrc_shift = self._recenter_mrc(reverse=True)
        self._recenter_conf(to_position=mrc_shift)

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
