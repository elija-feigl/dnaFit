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

from dataclasses import dataclass
from typing import Union

from .utils import h36_2_int, int_2_chimeraSegID, int_2_cifSegID, int_2_h36


@dataclass
class Number:
    """ read from hybrid36 string string or int """
    _input: Union[int, str]

    def __post_init__(self):
        self.n: int = self._convert_input(self._input)

    def _convert_input(self, inpt: Union[int, str]):
        if isinstance(inpt, int):
            return inpt
        elif isinstance(inpt, str):
            inpt.strip()
            if inpt.isdigit():
                return int(inpt)
            else:
                return h36_2_int(inpt)
        else:
            raise NotImplementedError

    def as_pdb4namd(self, width: int) -> str:
        return self.as_h36(width=width)

    def as_h36(self, width: int) -> str:
        return int_2_h36(number=self.n, width=width)

    def as_str(self) -> str:
        return str(self.n)


@dataclass
class AtomName:
    """ use cif standart """
    _input: str

    def __post_init__(self):
        self.NAMD = {"O1P": "OP1", "O2P": "OP2", "C5M": "C7", }
        self.name: str = self._convert_input(ipt=self._input)

    def _convert_input(self, ipt: str) -> str:
        ipt = ipt.strip()
        return self.NAMD[ipt] if ipt in self.NAMD.keys() else ipt

    def element_name(self) -> str:
        if len(self.name) == 1:
            return self.name
        else:
            # TODO: review double-letter elements
            return ''.join(filter(str.isalpha, self.name[:2]))[0]

    def as_pdb4namd(self, width: int) -> str:
        opt = self.name
        DMNA = {v: k for k, v in iter(self.NAMD.items())}
        namePDB = DMNA[opt] if opt in DMNA.keys() else opt
        if len(namePDB) == 1:
            return namePDB.ljust(2, " ").rjust(width, " ")
        elif len(namePDB) == 2:
            return namePDB.ljust(1, " ").rjust(width, " ")
        else:
            return namePDB.rjust(width, " ")

    def as_cif(self) -> str:
        if "'" in self.name:
            return f"\"{self.name}\""
        return self.name

    def as_str(self) -> str:
        return self.name


@dataclass
class ResName:
    """resname [DA, DC, DG, DT]"""
    _input: str

    def __post_init__(self):
        self.NAMD = {"CYT": "DC", "GUA": "DG", "THY": "DT", "ADE": "DA", }
        self.name: str = self._convert_input(inpt=self._input)

    def _convert_input(self, inpt: str) -> str:
        inpt = inpt.strip()
        if "5" in inpt:  # TODO: dirty fix for tacoxDNA
            inpt = inpt.replace("5", "")
        elif "3" in inpt:
            inpt = inpt.replace("3", "")
        return self.NAMD[inpt] if inpt in self.NAMD.keys() else inpt

    def as_pdb4namd(self, width: int) -> str:
        """returns resname [ADE, CYT, GUA, THY], len=3 """
        DMNA = {v: k for k, v in iter(self.NAMD.items())}
        return DMNA[self.name].rjust(width, " ")

    def as_str(self) -> str:
        return self.name

    def as_X(self) -> str:
        return self.name[-1]


@dataclass
class ChainID:
    _input: Union[int, str]

    def __post_init__(self):
        self.n: int = self._convert_input(self._input)

    def _convert_input(self, inpt: Union[int, str]):
        if isinstance(inpt, int):
            return inpt
        elif isinstance(inpt, str):
            inpt.strip()
            if inpt.isdigit():
                return int(inpt)
            else:
                return h36_2_int(inpt)
        else:
            raise NotImplementedError

    def as_pdb4namd(self, width: int) -> str:
        nlast = int(str(self.n)[0])
        n2char = "ABCDEFGHIJ"[nlast]
        return n2char.rjust(width, " ")

    def as_cif(self):
        """1-2 letter all caps"""
        return int_2_cifSegID(self.n - 1).rjust(2, " ")

    def as_chimera(self):
        """2 letter h36?"""
        return int_2_chimeraSegID(self.n - 1)

    def as_segName4namd(self, width: int) -> str:
        return "{:0{width}d}".format(self.n, width=width)

    def as_str(self) -> str:
        return str(self.n)
