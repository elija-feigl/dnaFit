#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
from typing import List, Tuple

import attr

from .atom import Atom
from .files import CIF, PDB


@attr.s
class Structure(object):
    path: Path = attr.ib()
    remove_H: bool = attr.ib(default=True)
    atoms: List[Atom] = list()
    # TODO: default and other attributes: sequences etc

    def __attrs_post_init__(self):
        self.name = self.filename.stem  # TODO pass name option
        self.keep_resID: bool = True
        self.previous_atm: Tuple[str, int] = ("", 0)
        self.previous_res: Tuple[str, int] = ("", 0)
        self.previous_chain: Tuple[str, int, bool] = (
            "", 0, True)  # str, int, new

    def add_atom(self, atom: Atom) -> None:
        self.atoms.append(atom)

    def _parse_cif_info(self):
        raise NotImplementedError

    def _parse_cif_atom(self):
        raise NotImplementedError

    def parse_cif(self) -> None:
        self._parse_cif_info()
        self._parse_cif_atom()
        raise NotImplementedError

    def _parse_pdb_info(self, line: str):
        # TODO: collect extra data
        return

    def _eval_atm_number(self, string: str, ) -> int:
        atm_number = self.previous_atm[1] + 1
        self.previous_atm = (string, atm_number)
        return atm_number

    def _eval_res_number(self, string: str) -> int:
        s = string.strip()
        if s.isdigit() and self.keep_resID:
            res_number = int(s)
        else:
            self.keep_resID = False
            if s == self.previous_res[0]:
                res_number = self.previous_res[1]
            else:
                res_number = self.previous_res[1] + 1

        if res_number == 1 and self.previous_res[1] != 1:
            self.previous_chain = (
                self.previous_chain[0], self.previous_chain[1], True)
        self.previous_res = (s, res_number)
        return res_number

    def _eval_chain_id(self, string: str) -> int:
        if string != self.previous_chain[0] or self.previous_chain[2]:
            chain_id = self.previous_chain[1] + 1
            self.previous_chain = (string, chain_id, False)
        else:
            chain_id = self.previous_chain[1]
        return chain_id

    def _parse_pdb_atom(self, line: str) -> None:
        atom_name = line[12:16]
        if "H" in atom_name and self.remove_H:
            return
        atm_number = self._eval_atm_number(line[6:11])
        res_number = self._eval_res_number(line[22:26])
        chain_id = self._eval_chain_id(line[20:22])

        self.add_atom(
            Atom(i_atom_coor=[line[30:38], line[38:46], line[46:54]],
                 i_atom_number=atm_number,
                 i_atom_name=atom_name,
                 i_res_name=line[17:20],
                 i_chain_id=chain_id,
                 i_res_number=res_number,
                 i_opacity=line[54:60].strip(),
                 i_temperature=line[60:66].strip(),
                 )
        )

    def _parse_pdb_generate_info(self):
        # generate sequences etc
        raise NotImplementedError

    def parse_pdb(self) -> None:
        with self.path.open() as fi:
            for line in fi.readlines():
                lineType = line[0:6].strip()
                if lineType == "ATOM":
                    self._parse_pdb_atom(line)
                else:
                    self._parse_pdb_info(line)

        # TODO: generate implicit data: sequences
        # self._parse_pdb_generate_info()

    def write_pdb(self, outfile: Path) -> None:
        pdb = PDB(struct=self)
        pdb.write(outfile=outfile)

    def write_cif(self, outfile: Path) -> None:
        # generate sequences here (from atom)
        cif = CIF(struct=self)
        cif.write(outfile=outfile)
