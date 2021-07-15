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
from typing import List, Tuple, Union

import numpy as np

from .types import AtomName, ChainID, Number, ResName


@dataclass
class Atom(object):
    i_atom_coor: Union[np.ndarray, Tuple[str, str, str],
                       Tuple[float, float, float]]
    i_atom_number: Union[int, str]
    i_atom_name: str
    i_res_name: str
    i_chain_id: Union[int, str]
    i_res_number: Union[int, str]
    i_opacity: Union[float, str] = 0.0
    i_temperature: Union[float, str] = 0.0

    def __post_init__(self):
        self.atom_coor: np.ndarray = self._convert_coor_input(self.i_atom_coor)
        self.atom_number: Number = Number(self.i_atom_number)
        self.atom_name: AtomName = AtomName(self.i_atom_name)
        self.res_name: ResName = ResName(self.i_res_name)
        self.chain_id: ChainID = ChainID(self.i_chain_id)
        self.res_number: Number = Number(self.i_res_number)
        self.opacity: float = (
            float(self.i_opacity) if float(self.i_opacity) != 0.0 else 1.0)
        self.temperature: float = float(self.i_temperature)
        self.element: str = self.atom_name.element_name()

    def _convert_coor_input(
            self, inpt: Union[np.ndarray, List[Union[str, float]]]
    ) -> np.ndarray:
        if isinstance(inpt, np.ndarray):
            return inpt
        elif isinstance(inpt, list):
            coor = [(float(c) if isinstance(c, str) else c) for c in inpt]
            return np.array(coor)
        else:
            raise NotImplementedError

    def asCif(self) -> str:
        return "".join([
            "ATOM ",
            self.atom_number.as_str().ljust(12, " "),
            self.element.ljust(2, " "),
            self.atom_name.as_cif().ljust(7, " "),
            ". ",
            self.res_name.as_str().ljust(3, " "),
            self.chain_id.as_cif().ljust(3, " "),
            self.chain_id.as_str().ljust(4, " "),
            self.res_number.as_str().ljust(6, " "),
            "? ",
            "{: .3f}".format(self.atom_coor[0]).ljust(9, " "),
            "{: .3f}".format(self.atom_coor[1]).ljust(9, " "),
            "{: .3f}".format(self.atom_coor[2]).ljust(9, " "),
            "{: .2f}".format(self.opacity).ljust(6, " "),
            "{: .2f}".format(self.temperature).ljust(8, " "),
            "? ",
            self.res_number.as_str().ljust(6, " "),
            self.res_name.as_str().ljust(3, " "),
            self.chain_id.as_chimera().ljust(3, " "),
            self.atom_name.as_cif().ljust(7, " "),
            "1",
            "\n",
        ])

    def asPdb(self) -> str:
        return "".join([
            "ATOM  ",
            self.atom_number.as_pdb4namd(width=5),
            self.atom_name.as_pdb4namd(width=5),
            self.res_name.as_pdb4namd(width=4),
            self.chain_id.as_pdb4namd(width=2),
            self.res_number.as_pdb4namd(width=4),
            (4 * " "),
            "{: .3f}".format(self.atom_coor[0]).rjust(8, " "),
            "{: .3f}".format(self.atom_coor[1]).rjust(8, " "),
            "{: .3f}".format(self.atom_coor[2]).rjust(8, " "),
            "{: .2f}".format(self.opacity).rjust(6, " "),
            "{: .2f}".format(self.temperature).rjust(6, " "),
            (6 * " "),
            self.chain_id.as_segName4namd(width=4),
            self.element,
            (2 * " "),  # charge
            "\n",
        ])
