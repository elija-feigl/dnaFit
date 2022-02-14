#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
""" class containing a minished atomic model for a specific cryo EM file.
        allow linking and corpping of the mrc map to a cadnano design
"""
import datetime
import logging
from dataclasses import dataclass
from pathlib import Path
from shutil import copyfile
from typing import Optional

import MDAnalysis as mda
from pdb2cif.pdb.structure import Structure

from ..data.mrc import write_mrc_from_atoms
from ..link.linkage import Linkage
from ..link.linker import Linker


@dataclass
class AtomicModelFit:
    """Atomic Model class generated form a specific cadnano design and sequence file,
    built to fit a specific cryp em mrc map.
    design positions and model indices can be linked
    """

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
            self.logger.debug("Auto-generate linkage using seq and json attributes of AtomicModel.")
            self.linkage = self._get_linkage(self.json, self.seq)

    def _write_logfile(self, prefix: str, dest: Path, **kwargs):
        log_path = dest / f"{prefix}_log.txt"
        with log_path.open(mode="w") as log_file:
            log_file.write(f"{datetime.datetime.now()}\n")
            log_file.write(f"topology: {self.top}\n")
            log_file.write(f"coordinate file: {self.conf}\n")
            log_file.write(f"fitted to cryo-EM data: {self.mrc}\n")
            for key, value in kwargs.items():
                log_file.write(f"{key}: {value}\n")

    def _get_linkage(self, json: Path, seq: Path) -> Linkage:
        """call Linker class to create linkage of design and atomic model"""
        linker = Linker(
            conf=self.conf,
            top=self.top,
            json=json,
            seq=seq,
            generated_with_mrdna=self.generated_with_mrdna,
        )
        return linker.create_linkage()

    def get_universe(self):
        """return MDAnalysis universe"""
        return mda.Universe(str(self.top), str(self.conf))

    def write_linkage(self, json: Path, seq: Path):
        """write persistent and human readable linkage (Fbp, FidDhps)"""
        out_link = json.parent / "dnaLink"

        Path(out_link).mkdir(parents=True, exist_ok=True)
        if self.linkage is None:
            self.linkage = self._get_linkage(json=json, seq=seq)

        self.linkage.write_linkage(prefix=json.stem, dest=out_link)
        self.logger.debug("writing linkage log file to %s", out_link)
        self._write_logfile(prefix=json.stem, dest=out_link, json=json, seq=seq)

    def write_output(self, dest: Path, write_mmcif=True, mask_mrc=True):
        """write atomic model and masked mrc file from."""

        def _copyfile(filepath: Path, dest):
            copyfile(filepath, dest / f"{filepath.stem}-AtomicModelFit{filepath.suffix}")

        if write_mmcif:
            # NOTE: Hydrogen removed is standart for RCSB upload
            structure = Structure(path=self.conf, remove_H=True)
            structure.parse_pdb()
            mmcif = self.conf.with_suffix(".cif")
            structure.write_cif(mmcif)
            _copyfile(mmcif, dest)
        else:
            _copyfile(self.top, dest)
            _copyfile(self.conf, dest)

        if mask_mrc:
            if self.mrc is None:
                self.logger.warning("no mrc file specified. Omitting mrc output")
            else:
                mrc_masked = self.mrc.with_name(f"{self.mrc.stem}-masked.mrc")
                if self.linkage is None:
                    universe = self.get_universe()
                else:
                    universe = self.linkage.u
                    write_mrc_from_atoms(
                        path=self.mrc,
                        atoms=universe.atoms,
                        path_out=mrc_masked,
                        context=10.0,
                        cut_box=True,
                    )
            copyfile(mrc_masked, dest / mrc_masked.name)
