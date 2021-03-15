import datetime

from pathlib import Path
from typing import Any, List, Optional
from shutil import copyfile

import attr
import mrcfile
import numpy as np
import MDAnalysis as mda

from ..link.linker import Linker
from ..link.linkage import Linkage
from ..data.mrc import write_mrc_from_atoms

""" class containing a minished atomic model for a specific cryo EM file.
        allow linking and corpping of the mrc map to a cadnano design
"""


@ attr.s
class AtomicModelFit(object):
    conf: Path = attr.ib()
    top: Path = attr.ib()
    mrc: Path = attr.ib()
    linkage: Optional[Linkage] = attr.ib(default=None)
    generated_with_mrdna: bool = attr.ib(default=True)

    def __attrs_post_init__(self) -> None:
        self.logger = logging.getLogger(__name__)

    def _write_logfile(self, prefix: str, dest: Path, **kwargs):
        log_file = dest / f"{prefix}_log.txt"
        with log_file.open(mode='w') as f:
            f.write(f"{datetime.datetime.now()}")
            f.write(f"topology: {self.top}")
            f.write(f"coordinate file: {self.conf}")
            f.write(f"fitted to cryo-EM data: {self.conf}")
            for key, value in kwargs.items():
                f.write(f"{key}: {value}")

    def _get_linkage(self, json: Path, seq: Path) -> Linkage:
        """ call Linker calss to create linkage of design and atomic model"""
        linker = Linker(conf=self.conf, top=self.top, json=json, seq=seq,
                        generated_with_mrdna=self.generated_with_mrdna)
        return linker.create_linkage()

    def write_linkage(self, json: Path, seq: Path):
        """ write persistent and human readable linkage (Fbp, FidDhps)"""
        out_link = json.parent / "dnaLink"

        Path(out_link).mkdir(parents=True, exist_ok=True)
        self.linkage = self._get_linkage(json=json, seq=seq)
        self.linkage.write_linkage(prefix=json.stem, dest=out_link)
        self.logger.debug(f"writing linkage log file to {out_link}")
        self._write_logfile(
            prefix=json.stem, dest=out_link, json=json, seq=seq)

    def write_output(self, dest: Path, write_mmCif=True, crop_mrc=True):
        def _copyfile(f, dest):
            copyfile(f, dest / f"{f.stem}-AtomicModelFit{f.suffix}")

        if crop_mrc:
            mrc_masked = self.mrc.with_name(f"{self.mrc.stem}-masked.mrc")
            if self.linkage is None:
                u = mda.Universe(str(self.top), str(self.conf))
            else:
                u = self.linkage.u
            write_mrc_from_atoms(path=self.mrc, atoms=u.atoms,
                                 path_out=mrc_masked, context=10., cut_box=True)

            copyfile(mrc_masked, dest / mrc_masked.name)

        if write_mmCif:
            # TODO: create mmcif for fit
            mmcif = self.conf
            raise NotImplementedError
            _copyfile(mmcif, dest)
        else:
            _copyfile(self.top, dest)
            _copyfile(self.conf, dest)
