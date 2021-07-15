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

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Dict, List, TextIO, Tuple

if TYPE_CHECKING:
    from .structure import Structure

from .. import get_resource
from .types import ChainID, ResName


@dataclass
class PDB(object):
    struct: Structure

    def __post_init__(self):
        # TODO: get box from struct
        self.box = "1000.000 1000.000 1000.000"

    def _header(self) -> List[str]:
        return ["HEADER   "]

    def _authors(self) -> List[str]:
        return ["AUTHOR   "]

    def _remarks(self) -> List[str]:
        return ["REMARK   "]

    def _box(self) -> List[str]:
        angles = "90.00  90.00  90.00"
        space_group = "P 1"
        z_value = "1"
        return [f"CRYST1 {self.box}  {angles} {space_group}           {z_value}"]

    def _atoms(self) -> List[str]:
        return [atom.asPdb() for atom in self.struct.atoms]

    def write(self, outfile: Path) -> None:
        with open(outfile, mode="w+") as of:
            for part in [
                    self._header(),
                    self._authors(),
                    self._remarks(),
                    self._box(),
                    self._atoms(),
            ]:
                of.writelines(part)


@dataclass
class CIF(object):
    struct: Structure

    def __post_init__(self):
        self.atoms: List[str] = self._set_atoms()
        self.chains, self.seqs = self._get_chains_seqs()

    def _set_atoms(self) -> List[str]:
        return [atom.asCif() for atom in self.struct.atoms]

    def _get_chains_seqs(self) -> Tuple[
            Dict[int, int], Dict[int, List[ResName]]]:
        chains: Dict[int, int] = dict()
        seqs: Dict[int, List[ResName]] = dict()
        C3ps = [a for a in self.struct.atoms if a.atom_name.as_str() == "C3\'"]
        for a in C3ps:
            chain_id = a.chain_id.n
            res_number = a.res_number.n

            if chain_id not in chains:
                chains[chain_id] = res_number
            elif res_number > chains[chain_id]:
                chains[chain_id] = res_number

            if chain_id not in seqs:  # only res not atom
                seqs[chain_id] = [a.res_name]
            else:
                seqs[chain_id].append(a.res_name)

        return chains, seqs

    def _write_atoms(self, fo: TextIO) -> None:
        atom_header = get_resource("cif_templates/atom_header.txt").read_text()
        fo.write(atom_header)
        fo.writelines(self.atoms)

    def _write_header(self, fo: TextIO) -> None:
        header = get_resource("cif_templates/header.txt").read_text()
        fo.write(f"data_{self.struct.name}\n")
        fo.write(header)

    def _write_pdbx_struct(self, fo: TextIO) -> None:
        pdbx_struct = get_resource("cif_templates/pdbx_struct.txt").read_text()
        pdbx_struct = pdbx_struct.replace("NCHAINS", str(len(self.chains)))
        fo.write(pdbx_struct)

    def _write_entity(self, fo: TextIO) -> None:
        entity = get_resource("cif_templates/entity.txt").read_text()
        fo.write(entity)
        for chain_id, length in self.chains.items():
            is_staple = (length < 500)
            n = str(chain_id).ljust(4)
            src = "syn" if is_staple else "?  "
            typ = "\'STAPLE STRAND\'  " if is_staple else "\'SCAFFOLD STRAND\'"
            fo.write(f"{n} polymer {src} {typ} ?   1 ? ? ? ?\n")

    def _write_entity_src(self, fo: TextIO) -> None:
        entity_src = get_resource("cif_templates/entity_src.txt").read_text()
        fo.write(entity_src)
        for chain_id, length in self.chains.items():
            is_staple = (length < 500)
            n = str(chain_id).ljust(4)
            typ = "\'synthetic construct\'" if is_staple else "?".ljust(21)
            tax = "32630" if is_staple else "?".ljust(5)
            fo.write(
                f"{n}   1 sample 1 {str(length).ljust(5)} {typ} ? {tax} ?\n")

    def _write_entity_poly(self, fo: TextIO) -> None:
        entity_poly = get_resource("cif_templates/entity_poly.txt").read_text()
        fo.write(entity_poly)
        for chain_id, seq in self.seqs.items():
            n = str(chain_id).ljust(4)
            cid = ChainID(chain_id).as_chimera()
            seq1 = "".join([f"({s.as_str()})" for s in seq])
            seq2 = "".join([s.as_X() for s in seq])
            fo.write(
                f"{n} polydeoxyribonucleotide no no\n;{seq1}\n;\n{seq2} {cid} ?\n")

    def write(self, outfile: Path) -> None:
        with open(outfile, mode="w+") as fo:
            self._write_header(fo)
            self._write_pdbx_struct(fo)
            self._write_atoms(fo)
            self._write_entity(fo)
            self._write_entity_src(fo)
            self._write_entity_poly(fo)
