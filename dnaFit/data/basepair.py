#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
""" BasePair Class represents a watson-crick baspair of two nanodesign base
    object. Important Attributes are their position in the design-file and
    their spatial orientation in real space (BasePlane and BasePairPlane class)
"""
from dataclasses import dataclass
from typing import Any
from typing import Dict
from typing import Optional
from typing import Tuple

import numpy as np
import numpy.typing as npt
from MDAnalysis.core.groups import Residue

from ..core.utils import _norm


@dataclass(frozen=True)
class BasePairPlane:
    """plane_versor: plane-normal vector. always pointing in scaffold 5'->3' direction
    wc_vector: vector pointing from scaffold to staple
    """

    __slots__ = ["positions", "wc_vectors", "plane_versor"]
    positions: Dict[str, Any]
    wc_vectors: Dict[str, Any]
    plane_versor: npt.NDArray[np.float64]


@dataclass(frozen=True)
class BasePlane:
    """plane_versor: plane-normal vector. always pointing in scaffold 5'->3' direction"""

    __slots__ = ["positions", "plane_versor"]
    positions: Dict[str, Any]
    plane_versor: npt.NDArray[np.float64]


@dataclass
class BasePair:
    """every square of the JSON can be represented as BP"""

    scaffold: Residue
    staple: Residue
    hp: Tuple[int, int]

    sc_plane: Optional[BasePlane] = None
    st_plane: Optional[BasePlane] = None
    plane: Optional[BasePairPlane] = None

    def __post_init__(self):
        if self.scaffold is None or self.staple is None:
            self.is_ds = False
        else:
            self.is_ds = True

    def calculate_baseplanes(self):
        """calculate base planes and basepair planes"""
        self.sc_plane = (
            self._get_base_plane(res=self.scaffold, is_scaf=True)
            if self.scaffold is not None
            else None
        )
        self.st_plane = (
            self._get_base_plane(res=self.staple, is_scaf=False)
            if self.staple is not None
            else None
        )
        self.plane = (
            self._get_bp_plane(scaffold=self.sc_plane, staple=self.st_plane) if self.is_ds else None
        )

    def _get_base_plane(self, res: Residue, is_scaf: bool) -> BasePlane:
        positions = dict()
        atom = []
        for atom_name in ["C2", "C4", "C6"]:
            (atom_select,) = res.atoms.select_atoms("name " + atom_name)
            atom.append(atom_select.position)

        plane_versor = _norm(np.cross((atom[1] - atom[0]), (atom[2] - atom[0])))
        if res.resname in ["ADE", "GUA"] and is_scaf:
            plane_versor = -plane_versor
        elif res.resname in ["THY", "CYT"] and not is_scaf:
            plane_versor = -plane_versor

        positions["diazine"] = sum(atom) / 3.0

        c6c8 = "C8" if res.resname in ["ADE", "GUA"] else "C6"
        positions["C6C8"] = res.atoms.select_atoms("name " + c6c8)[0].position

        positions["C1'"] = res.atoms.select_atoms("name C1'")[0].position

        return BasePlane(plane_versor=plane_versor, positions=positions)

    def _get_bp_plane(self, scaffold, staple) -> BasePairPlane:
        wc_vectors, positions = dict(), dict()
        plane_versor = (scaffold.plane_versor + staple.plane_versor) * 0.5
        for pos in scaffold.positions:
            positions[pos] = (scaffold.positions[pos] + staple.positions[pos]) * 0.5
            wc_vectors[pos] = staple.positions[pos] - scaffold.positions[pos]

        return BasePairPlane(plane_versor=plane_versor, wc_vectors=wc_vectors, positions=positions)
