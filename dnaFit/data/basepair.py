#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
""" BasePair Class represents a watson-crick baspair of two nanodesign base
    object. Important Attributes are their position in the design-file and
    their spatial orientation in real space.
    Olson et al. (2001). A standard reference frame for the description of nucleic acid base-pair geometry.
    """
import logging
from dataclasses import dataclass
from dataclasses import field
from typing import Tuple

import numpy as np
import numpy.typing as npt
from MDAnalysis.core.groups import Residue

from ..core.utils import _norm
from ..core.utils import _project_v2plane


@dataclass(frozen=True)
class Plane:
    """position: base or basepair center point
    normal: plane-normal vector. always pointing in scaffold 5'->3' direction
    direction: vector pointing away from C1' or from scaffold to staple
    """

    __slots__ = ["origin", "direction", "normal"]
    origin: npt.NDArray[np.float64]
    direction: npt.NDArray[np.float64]  # y-vector in  Olson et al. (2001), norm = C1'-C1'
    normal: npt.NDArray[np.float64]  # z-versor in  Olson et al. (2001)

    @property
    def vector3(self) -> npt.NDArray[np.float64]:
        """x-versor in Olson et al. (2001)"""
        return _norm(np.cross(self.normal, self.direction))


@dataclass
class BasePair:
    """every square of the JSON can be represented as BP"""

    scaffold: Residue
    staple: Residue
    hp: Tuple[int, int]

    sc_plane: Plane = field(init=False)
    st_plane: Plane = field(init=False)
    plane: Plane = field(init=False)

    def __post_init__(self):
        self.logger = logging.getLogger(__name__)

    def calculate_baseplanes(self):
        """calculate base planes and basepair planes"""
        self.plane = self._get_bp_plane()

        self.sc_plane = self._get_base_plane(res=self.scaffold, is_scaf=True)

        self.st_plane = self._get_base_plane(res=self.staple, is_scaf=False)

    @staticmethod
    def _atom_position(res: Residue, name: str) -> npt.NDArray[np.float64]:
        return res.atoms.select_atoms(f"name {name}")[0].position

    def _get_bp_plane(self) -> Plane:
        def c6c8_position(res: Residue) -> npt.NDArray[np.float64]:
            atom_name = "C8" if res.resname in ["ADE", "GUA"] else "C6"
            return self._atom_position(res=res, name=atom_name)

        st_c1p = self._atom_position(self.staple, "C1'")
        sc_c1p = self._atom_position(self.scaffold, "C1'")
        direction = sc_c1p - st_c1p
        dyad_point = (sc_c1p + st_c1p) * 0.5

        sc_c6c8 = c6c8_position(self.scaffold)
        origin_direction = sc_c6c8 - c6c8_position(self.staple)

        # intersect c6-c8 line with pseudo-dyad plane of C1'-C1'
        projection_on_dyad_plane = np.inner(direction, origin_direction)
        w = dyad_point - sc_c6c8
        si = np.inner(direction, w) / projection_on_dyad_plane
        origin = dyad_point - w + si * direction

        normal = _norm(np.cross((origin - dyad_point), direction))
        return Plane(origin=origin, direction=direction, normal=normal)

    def _get_base_plane(self, res: Residue, is_scaf: bool) -> Plane:
        """projection on pyrimidine plane.
        center in pyrimidine"""
        pyr_names = ["N3", "N1", "C5"] if res.resname in ["ADE", "GUA"] else ["N1", "N3", "C5"]

        atom = [self._atom_position(res, name) for name in pyr_names]
        normal = _norm(np.cross((atom[1] - atom[0]), (atom[2] - atom[0])))
        normal = normal if is_scaf else -normal
        origin = sum(atom) / 3.0

        direction = _project_v2plane(self.plane.direction, normal)
        return Plane(origin=origin, direction=direction, normal=normal)
