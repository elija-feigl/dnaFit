#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
"""
    BDna Class handles evaluation of B-DNA properties and local-resolution of
    cryo em reconstruction of specific bases
    Olson et al. (2001). A standard reference frame for the description of nucleic acid base-pair geometry.
"""
from dataclasses import dataclass
from dataclasses import field
from typing import Any
from typing import Dict
from typing import List
from typing import Optional
from typing import Tuple

import numpy as np
import numpy.typing as npt
from MDAnalysis.core.groups import AtomGroup

from ..core.utils import _norm_proj
from ..core.utils import _proj
from ..core.utils import _save_arccos_deg
from ..link.linkage import Linkage
from .basepair import BasePair

# from ..core.utils import _norm
# from ..core.utils import _v_proj

# from crossover import Crossover
# import MDAnalysis as mda
# from ..core.utils import (
#     C1P_BASEDIST, WC_HBONDS, WC_HBONDS_DIST, BB_ATOMS,
#     PUR_ATOMS, PYR_ATOMS, DH_ATOMS,
#     _dh_angle, _proj2plane
# )


@dataclass
class BDna:
    link: Linkage

    bps: Dict[Tuple[int, int], BasePair] = field(default_factory=dict)

    bp_trace: Dict[int, npt.NDArray[np.float64]] = field(default_factory=dict)

    bp_twist: Dict[int, float] = field(default_factory=dict)
    bp_rise: Dict[int, float] = field(default_factory=dict)
    bp_tilt: Dict[int, float] = field(default_factory=dict)
    bp_roll: Dict[int, float] = field(default_factory=dict)
    bp_shift: Dict[int, float] = field(default_factory=dict)
    bp_slide: Dict[int, float] = field(default_factory=dict)

    # b_stretch: Dict[int, float] = field(default_factory=dict) TODO-low: implement base plane defintion first
    # b_stagger: Dict[int, float] = field(default_factory=dict)
    # b_shear: Dict[int, float] = field(default_factory=dict)
    # b_open: Dict[int, float] = field(default_factory=dict)
    # b_buckle: Dict[int, float] = field(default_factory=dict)
    # b_propeller: Dict[int, float] = field(default_factory=dict)

    dh_backbone: Dict[int, Any] = field(default_factory=dict)

    # distances: Dict[int, Any] = field(default_factory=dict) TODO-low: any value/interest?
    co_angles: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.bps = self._create_bps()
        self.eval_basepair()
        # self.eval_base()

        # self.eval_distances()
        # self.eval_dh()
        # self.eval_co_angles() TODO-low: wait for relevance

    def _create_bps(self) -> Dict[Tuple[int, int], BasePair]:
        bps = dict()
        for scaffold_resid, staple_resid in self.link.Fbp.items():
            scaffold = self.link.u.residues[scaffold_resid]
            staple = self.link.u.residues[staple_resid]
            h, p, _ = self.link.DidDhps[self.link.FidDid[scaffold_resid]]

            bp = BasePair(scaffold=scaffold, staple=staple, hp=(h, p))
            bp.calculate_baseplanes()
            bps[(h, p)] = bp
        return bps

    def _get_n_bp(self, bp: BasePair, steps: int = 1) -> Optional[BasePair]:
        if steps == 0:
            return bp
        helix, position = bp.hp
        if (helix % 2) == 1:  # and local:
            steps = -steps
        direct = np.sign(steps)
        n_position = position + steps

        # first check the number of skips passed
        n_skips = 0
        for n in range(direct, direct * (steps + 1), direct):
            n_position = position + n
            if (helix, n_position) not in self.bps.keys():
                n_skips += 1
        # move one position further if on skip
        n_position = position + direct * (steps + n_skips)
        if (helix, n_position) not in self.bps.keys():
            n_position += direct

        return self.bps.get((helix, n_position), None)

    def eval_basepair(self) -> None:
        """Affects
        -------
            self.bp_...
        """
        for bp in self.bps.values():

            resindex = bp.scaffold.resindex
            self.bp_trace[resindex] = bp.plane.origin

            # scaffold 5'->3'
            n_bp = self._get_n_bp(bp=bp)  # , local=True)
            if n_bp is not None:
                self.bp_rise[resindex] = self.get_bp_rise(bp, n_bp)
                self.bp_shift[resindex] = self.get_bp_shift(bp, n_bp)
                self.bp_slide[resindex] = self.get_bp_slide(bp, n_bp)

                self.bp_twist[resindex] = self.get_bp_twist(bp, n_bp)
                self.bp_roll[resindex] = self.get_bp_roll(bp, n_bp)
                self.bp_tilt[resindex] = self.get_bp_tilt(bp, n_bp)

    def eval_base(self) -> None:
        raise NotImplementedError

    # def _get_bp_quality(self, bp: BasePair) -> Dict[str, Any]:
    #     quality = dict()
    #     # TODO: cleanup (loop?)
    #     atom_name = "C1'"
    #     bnd = atom_name * 2
    #     C1p = []
    #     atoms = []
    #     for b in [bp.sc, bp.st]:
    #         C1p.append(b.atoms.select_atoms("name {}".format(atom_name))[0])
    #         atoms.append(
    #             b.atoms.select_atoms(
    #                 "name " + ' or name '.join(map(str, WC_HBONDS[b.resname]))
    #             )
    #         )

    #     c1p_dist = mdamath.norm(C1p[0].position - C1p[1].position)
    #     c1p_dev = abs(c1p_dist - C1P_BASEDIST) / C1P_BASEDIST
    #     quality[bnd] = c1p_dev

    #     # TODO: review hbond-dev, move to function
    #     for idx, a1 in enumerate(atoms[0]):
    #         hbond_dist = mdamath.norm(atoms[1][idx].position - a1.position)
    #         hbond_should = WC_HBONDS_DIST[bp.sc.resname][idx]
    #         hbond_dev = abs(hbond_dist - hbond_should) / hbond_should
    #         bnd = WC_HBONDS[bp.sc.resname][idx] + WC_HBONDS[bp.st.resname][idx]
    #         quality[bnd] = (hbond_dev if hbond_dist > hbond_should
    #                         else -hbond_dev)
    #     return quality

    @staticmethod
    def get_bp_twist(bp: BasePair, n_bp: BasePair) -> float:
        cos_phi = _norm_proj(bp.plane.direction, n_bp.plane.direction)
        return _save_arccos_deg(cos_phi)

    @staticmethod
    def get_bp_tilt(bp: BasePair, n_bp: BasePair) -> float:
        cos_phi = _norm_proj(bp.plane.normal, n_bp.plane.normal)
        return _save_arccos_deg(cos_phi)

    @staticmethod
    def get_bp_roll(bp: BasePair, n_bp: BasePair) -> float:
        cos_phi = _norm_proj(bp.plane.vector3, n_bp.plane.vector3)
        return _save_arccos_deg(cos_phi)

    @staticmethod
    def get_bp_rise(bp: BasePair, n_bp: BasePair) -> float:
        trace = n_bp.plane.origin - bp.plane.origin  # TODO: check abs
        return np.abs(_proj(trace, bp.plane.normal))  # type: ignore

    @staticmethod
    def get_bp_slide(bp: BasePair, n_bp: BasePair) -> float:
        trace = n_bp.plane.origin - bp.plane.origin
        return _proj(trace, bp.plane.direction)  # type: ignore

    @staticmethod
    def get_bp_shift(bp: BasePair, n_bp: BasePair) -> float:
        trace = n_bp.plane.origin - bp.plane.origin
        return _proj(trace, bp.plane.vector3)  # type: ignore

    # def eval_dh(self) -> None:
    #     """ Affects
    #         -------
    #             self.dh_quality
    #     """
    #     for res in self.link.u.residues:
    #         self.dh_quality[res.resindex] = self._get_dihedrals(res)

    # # TODO: improve
    # def _get_dihedrals(self, res: "mda.residue") -> Dict[str, float]:
    #     def _get_residue_BB(res: "mda.residue"
    #                         ) -> Tuple[Dict[str, "np.ndarray"],
    #                                    Tuple[bool, bool, bool]
    #                                    ]:
    #         iniSeg, terSeg, ter5 = False, False, False

    #         atoms = dict()
    #         try:
    #             P = res.atoms.select_atoms("name P")[0]
    #             atoms["P"] = P.position
    #         except (KeyError, IndexError):
    #             ter5 = True

    #         for x in BB_ATOMS[:-1]:
    #             atom_name = "name {}".format(x)
    #             atoms[x] = res.atoms.select_atoms(atom_name)[0].position
    #         if res.resname in ["ADE", "GUA"]:
    #             Y = PUR_ATOMS
    #         else:
    #             Y = PYR_ATOMS
    #         for y in Y:
    #             atom_name = "name {}".format(y)
    #             atoms[y] = res.atoms.select_atoms(atom_name)[0].position

    #         try:
    #             n_res = self.link.u.residues[res.resindex + 1]
    #             if res.segindex == n_res.segindex:
    #                 n_P = n_res.atoms.select_atoms("name P")[0]
    #                 n_O5p = n_res.atoms.select_atoms("name O5'")[0]
    #                 atoms["P +"] = n_P.position
    #                 atoms["O5' +"] = n_O5p.position
    #             else:
    #                 terSeg = True
    #         except (KeyError, IndexError):
    #             terSeg = True

    #         try:
    #             p_res = self.link.u.residues[res.resindex - 1]
    #             if res.segindex == p_res.segindex:
    #                 p_O3p = p_res.atoms.select_atoms("name O3'")[0]
    #                 atoms["O3' -"] = p_O3p.position
    #             else:
    #                 iniSeg = True
    #         except (KeyError, IndexError):
    #             iniSeg = True

    #         return atoms, (ter5, terSeg, iniSeg)

    #     def _get_valid_DH(ter5: bool, terSeg: bool, iniSeg: bool) -> List[str]:
    #         dh_valid = ["gamma", "delta", "xi"]
    #         if terSeg is False:
    #             dh_valid.extend(["epsilon", "zeta"])
    #         if iniSeg is False:
    #             dh_valid.append("alpha")
    #         if ter5 is False:
    #             dh_valid.append("beta")
    #         return dh_valid

    #     def _get_dh_for_res(atoms: Dict[str, "np.ndarray"],
    #                         pyr: bool,
    #                         dh_valid: List[str],
    #                         ) -> Dict[str, float]:
    #         dh = {}
    #         for dh_name in DH_ATOMS:
    #             if dh_name in dh_valid:
    #                 angle = _get_dhangle(atoms, pyr, dh_name)
    #             else:
    #                 angle = np.nan
    #             dh[dh_name] = angle
    #         return dh

    #     def _get_dhangle(atoms: Dict[str, "np.ndarray"],
    #                      pyr: bool,
    #                      dh_name: str,
    #                      ) -> float:
    #         p = []
    #         if dh_name == "xi":
    #             if pyr:
    #                 for i in DH_ATOMS[dh_name]["pyr"]:
    #                     p.append(atoms[i])
    #             else:
    #                 for i in DH_ATOMS[dh_name]["pur"]:
    #                     p.append(atoms[i])
    #         else:
    #             for i in DH_ATOMS[dh_name]:
    #                 p.append(atoms[i])
    #         return _dh_angle(p)

    #     atoms, logic = _get_residue_BB(res)
    #     dh_valid = _get_valid_DH(*logic)
    #     if res.resname in ["ADE", "GUA"]:
    #         pyr = False
    #     else:
    #         pyr = True

    #     dh = _get_dh_for_res(atoms, pyr, dh_valid)
    #     return dh

    # def eval_distances(self) -> None:
    #     """ Affects
    #         -------
    #             self.distances
    #     """
    #     for bp in self.bps.values():
    #         if bp.sc is not None:
    #             self.distances[bp.sc.resindex] = dict()
    #         if bp.st is not None:
    #             self.distances[bp.st.resindex] = dict()
    #         for atom in ["C1'", "P"]:
    #             sc_dist, st_dist = self._get_distance(bp=bp, name=atom)
    #             if bp.sc is not None:
    #                 self.distances[bp.sc.resindex][atom] = sc_dist
    #             if bp.st is not None:
    #                 self.distances[bp.st.resindex][atom] = st_dist

    # def _get_distance(self, bp: BasePair, name: str
    #                   ) -> Tuple[Dict[str, float], Dict[str, float]]:
    #     def try_pos(res: "mda.residue",
    #                 name: str,
    #                 ) -> Optional["np.ndarray"]:
    #         try:
    #             return res.atoms.select_atoms(name)[0].position
    #         except (AttributeError, IndexError):
    #             return None

    #     sc_dist: Dict[str, float] = dict()
    #     st_dist: Dict[str, float] = dict()
    #     atom_name = "name {}".format(name)

    #     n_bp = self._get_n_bp(bp=bp, steps=1)
    #     p_bp = self._get_n_bp(bp=bp, steps=-1)
    #     for s, x_bp, dist in [("sc", n_bp, sc_dist), ("st", p_bp, st_dist)]:
    #         if x_bp is None:
    #             X = try_pos(res=bp.sc, name=atom_name)
    #             Y = try_pos(res=bp.st, name=atom_name)
    #             if X is not None and Y is not None:
    #                 dist["pair"] = np.linalg.norm(Y - X)
    #                 dist["stack"] = np.nan
    #                 dist["crossstack"] = np.nan
    #             else:
    #                 for typ in ["pair", "stack", "crossstack"]:
    #                     dist[typ] = np.nan
    #         else:
    #             SC = [try_pos(res=x.sc, name=atom_name) for x in [bp, x_bp]]
    #             ST = [try_pos(res=x.st, name=atom_name) for x in [bp, x_bp]]
    #             if s == "st":
    #                 SC, ST = ST, SC
    #             for X, Y, typ in [(SC[0], ST[0], "pair"),
    #                               (SC[0], SC[1], "stack"),
    #                               (SC[0], ST[1], "crossstack")
    #                               ]:
    #                 if X is not None and Y is not None:
    #                     dist[typ] = np.linalg.norm(Y - X)
    #                 else:
    #                     dist[typ] = np.nan
    #     return sc_dist, st_dist

    # def eval_co_angles(self) -> None:
    #     """ Definition: Bai, X. (2012).  doi: 10.1073/pnas.1215713109
    #         Affects
    #         -------
    #             self.co_angles
    #     """
    #     def get_co_angles(co: Crossover, typ="C6C8"):
    #         def P_from_p(p: Optional[BasePair], typ: str) -> "np.ndarray":
    #             if p is None:
    #                 return np.zeros(3)
    #             if p.plane is not None:
    #                 return p.plane.P[typ]
    #             elif p.sc is not None:
    #                 return p.sc_plane.P[typ]
    #             else:
    #                 return p.st_plane.P[typ]

    #         resindices = list()
    #         for P in co.Ps:
    #             if P is None:
    #                 continue
    #             elif P.sc is not None:
    #                 resindices.append(P.sc.resindex)
    #             elif P.st is not None:
    #                 resindices.append(P.st.resindex)

    #         X = [P_from_p(p, typ) for p in co.Ps]
    #         X_ = [P_from_p(l, typ) for l in co.Ls]

    #         center = sum(X) / len(X)
    #         x = [(p_ - p) for p, p_ in zip(X, X_)]
    #         n0 = _norm(X[0] + X[1] - X[2] - X[3])
    #         proj = _proj2plane(x, n0)

    #         d1 = _proj(proj[0], proj[2])
    #         a_n0 = _save_arccos_deg(_proj(x[0], n0))

    #         if co.typ != "end":
    #             b_n0 = _save_arccos_deg(_proj(x[1], n0))
    #             d2 = _proj(proj[3], proj[1])
    #             a1 = _proj(proj[0], [-x for x in proj[1]])
    #             a2 = _proj(proj[2], [-x for x in proj[3]])

    #             gamma1 = _save_arccos_deg(d1)
    #             gamma2 = _save_arccos_deg(d2)
    #             alpha1 = _save_arccos_deg(a1)
    #             alpha2 = _save_arccos_deg(a2)
    #         else:
    #             b_n0 = 90.
    #             gamma2 = np.nan
    #             alpha1 = np.nan
    #             alpha2 = np.nan

    #         gamma1 = np.rad2deg(np.arccos(d1))
    #         beta = 180. - abs(a_n0) - abs(b_n0)

    #         return {"angles": {"co_beta": beta,
    #                            "co_gamma1": gamma1,
    #                            "co_gamma2": gamma2,
    #                            "co_alpha1": alpha1,
    #                            "co_alpha2": alpha2,
    #                            },
    #                 "center-co": center,
    #                 "plane": n0,
    #                 "resindices": resindices,
    #                 }

    #     for key, co in self.link.Fco.items():
    #         co_data = get_co_angles(co=co)

    #         self.co_angles[key] = {
    #             "co": key,
    #             "type": co.typ,
    #             "is_scaffold": co.is_scaf,
    #             "angles": co_data["angles"],
    #             "center-co": co_data["center-co"],
    #             "resindices": co_data["resindices"],
    #         }

    def pair_resid_selection(self, atoms=None) -> List[int]:
        residues = self.link.u.residues if atoms is None else atoms.residues

        atoms_sc = AtomGroup([], self.link.u)
        for residue in residues:
            if residue.resindex in self.link.Fbp.keys():
                atoms_sc += residue.atoms
        return list(atoms_sc.resids)
