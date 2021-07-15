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
from typing import Dict, Tuple

import numpy as np
from MDAnalysis.core.groups import Residue

from ..core.utils import _norm

""" BasePair Class represents a watson-crick baspair of two nanodesign base
    object. Important Attributes are their position in the design-file and
    their spatial orientation in real space (BasePlane and BasePairPlane class)

    COMMENTS:
"""


@dataclass(frozen=True)
class BasePairPlane(object):
    """n0: plane-normal vector. always pointing in scaffold 5'->3' direction
    """
    __slots__ = ["P", "a", "n0"]
    P: Dict[str, np.ndarray]
    a: Dict[str, np.ndarray]
    n0: np.ndarray


@dataclass(frozen=True)
class BasePlane(object):
    """n0: plane-normal vector. always pointing in scaffold 5'->3' direction
    """
    __slots__ = ["P", "n0"]
    P: Dict[str, np.ndarray]
    n0: np.ndarray


@dataclass
class BasePair(object):
    """ every square of the JSON can be represented as BP
    """
    sc: Residue
    st: Residue
    hp: Tuple[int, int]

    def __post_init__(self):
        if self.sc is None or self.st is None:
            self.is_ds = False
        else:
            self.is_ds = True

    def calculate_baseplanes(self):
        self.sc_plane = (self._get_base_plane(res=self.sc, is_scaf=True)
                         if self.sc is not None else None)
        self.st_plane = (self._get_base_plane(res=self.st, is_scaf=False)
                         if self.st is not None else None)
        self.plane = (self._get_bp_plane(sc=self.sc_plane, st=self.st_plane)
                      if self.is_ds else None)

    def _get_base_plane(self, res: Residue, is_scaf: bool) -> BasePlane:
        P = dict()
        atom = []
        for atom_name in ["C2", "C4", "C6"]:
            A, = res.atoms.select_atoms("name " + atom_name)
            atom.append(A.position)

        n0 = _norm(np.cross((atom[1] - atom[0]), (atom[2] - atom[0])))
        if res.resname in ["ADE", "GUA"] and is_scaf:
            n0 = -n0
        elif res.resname in ["THY", "CYT"] and not is_scaf:
            n0 = -n0

        P["diazine"] = sum(atom) / 3.

        C6C8 = "C8" if res.resname in ["ADE", "GUA"] else "C6"
        P["C6C8"] = res.atoms.select_atoms("name " + C6C8)[0].position

        P["C1'"] = res.atoms.select_atoms("name C1'")[0].position

        return BasePlane(n0=n0, P=P)

    def _get_bp_plane(self, sc, st) -> BasePairPlane:
        a, P = dict(), dict()
        n0 = (sc.n0 + st.n0) * 0.5
        for n in sc.P:
            P[n] = (sc.P[n] + st.P[n]) * 0.5
            a[n] = st.P[n] - sc.P[n]

        return BasePairPlane(n0=n0, a=a, P=P)
