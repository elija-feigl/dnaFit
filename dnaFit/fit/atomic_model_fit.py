#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2021  Elija Feigl
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see https://www.gnu.org/licenses/gpl-3.0.html.

import datetime
import logging
from dataclasses import dataclass
from pathlib import Path
from shutil import copyfile
from typing import Optional

import MDAnalysis as mda

from ..data.mrc import write_mrc_from_atoms
from ..link.linkage import Linkage
from ..link.linker import Linker
from ..pdb.structure import Structure

""" class containing a minished atomic model for a specific cryo EM file.
        allow linking and corpping of the mrc map to a cadnano design
"""


@dataclass
class AtomicModelFit(object):
    conf: Path
    top: Path
    mrc: Optional[Path] = None
    json: Optional[Path] = None
    seq: Optional[Path] = None
    linkage: Optional[Linkage] = None
    generated_with_mrdna: bool = True

    def __post_init__(self) -> None:
        self.logger = logging.getLogger(__name__)
        if self.linkage is not None:
            self.logger.debug("Using existing Linkage.")
        elif self.json is not None and self.seq is not None:
            self.logger.debug(
                "Auto-generate linkage using seq and json attributes of AtomicModel.")
            self.linkage = self._get_linkage(self.json, self.seq)

    def _write_logfile(self, prefix: str, dest: Path, **kwargs):
        log_file = dest / f"{prefix}_log.txt"
        with log_file.open(mode='w') as f:
            f.write(f"{datetime.datetime.now()}\n")
            f.write(f"topology: {self.top}\n")
            f.write(f"coordinate file: {self.conf}\n")
            f.write(f"fitted to cryo-EM data: {self.mrc}\n")
            for key, value in kwargs.items():
                f.write(f"{key}: {value}\n")

    def _get_linkage(self, json: Path, seq: Path) -> Linkage:
        """ call Linker class to create linkage of design and atomic model"""
        linker = Linker(conf=self.conf, top=self.top, json=json, seq=seq,
                        generated_with_mrdna=self.generated_with_mrdna)
        return linker.create_linkage()

    def write_linkage(self, json: Path, seq: Path):
        """ write persistent and human readable linkage (Fbp, FidDhps)"""
        out_link = json.parent / "dnaLink"

        Path(out_link).mkdir(parents=True, exist_ok=True)
        if self.linkage is None:
            self.linkage = self._get_linkage(json=json, seq=seq)

        self.linkage.write_linkage(prefix=json.stem, dest=out_link)
        self.logger.debug(f"writing linkage log file to {out_link}")
        self._write_logfile(
            prefix=json.stem, dest=out_link, json=json, seq=seq)

    def write_output(self, dest: Path, mrc=None, write_mmCif=True, crop_mrc=True):
        def _copyfile(f, dest):
            copyfile(f, dest / f"{f.stem}-AtomicModelFit{f.suffix}")

        if crop_mrc:
            if mrc is None:
                if self.mrc is None:
                    self.logger.warn("no mrc file specified.")
                    return
                else:
                    self.logger.debug("Using mrc attribute of AtomicModel.")
                    mrc = self.mrc

            mrc_masked = mrc.with_name(f"{mrc.stem}-masked.mrc")
            if self.linkage is None:
                u = self.get_universe()
            else:
                u = self.linkage.u
            write_mrc_from_atoms(path=mrc, atoms=u.atoms,
                                 path_out=mrc_masked, context=10., cut_box=True)

            copyfile(mrc_masked, dest / mrc_masked.name)

        if write_mmCif:
            # NOTE: Hydrogen removed is standart for RCSB upload
            structure = Structure(path=self.conf, remove_H=True)
            structure.parse_pdb()
            mmcif = self.conf.with_suffix(".cif")
            structure.write_cif(mmcif)
            _copyfile(mmcif, dest)
        else:
            _copyfile(self.top, dest)
            _copyfile(self.conf, dest)

    def get_universe(self):
        return mda.Universe(str(self.top), str(self.conf))
