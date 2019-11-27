#!/usr/bin/env python
# -*- coding: utf-8 -*-3
# TODO: abkürzungsverzeichnis

import numpy as np
import MDAnalysis as mda
import mrcfile as mrc
from MDAnalysis.lib import mdamath

import attr
from typing import Dict, Tuple, Any, Optional, List, Set

from linker import Linkage
from basepair import BasePair
from crossover import Crossover
from utils import (C1P_BASEDIST, WC_HBONDS, WC_HBONDS_DIST, BB_ATOMS,
                   PUR_ATOMS, PYR_ATOMS, DH_ATOMS,
                   _proj, _norm, _v_proj, _save_arccos_deg,
                   _dh_angle, _proj2plane
                   )

_author__ = "Elija Feigl"
__copyright__ = "Copyright 2019, Dietzlab (TUM)"
__credits__ = ["Autodesk: Nanodesign", "MDAnalysis", "mrcfile"]
__license__ = "None"
__version__ = "0.4"
__maintainer__ = "Elija Feigl"
__email__ = "elija.feigl@tum.de"
__status__ = "Development"


@attr.s
class BDna(object):
    link: Linkage = attr.ib()

    def __attrs_post_init__(self) -> None:
        self.bps: Dict[Tuple[int, int], BasePair] = self._get_pot_bp()
        self.link.relink_crossover_basepairs(self.bps)
        self.bp_quality: Dict[int, Any] = {}
        self.bp_geometry_local: Dict[int, Any] = {}
        self.bp_geometry_global: Dict[int, Any] = {}
        self.dh_quality: Dict[int, Any] = {}
        self.distances: Dict[int, Any] = {}
        self.co_angles: Dict[str, Any] = {}

    def sample(self) -> None:
        for bp in self.bps.values():
            bp.calculate_baseplanes()

        self.eval_bp()
        self.eval_distances()
        self.eval_dh()
        self.eval_co_angles()

    def _get_bp(self, resindex: int) -> Tuple[Tuple[int, int], BasePair]:
        h, p, is_scaf = self.link.DidDhps[self.link.FidDid[resindex]]
        res = self.link.u.residues[resindex]
        wcindex = self.link.Fbp_full.get(resindex, None)
        wc = None if wcindex is None else self.link.u.residues[wcindex]

        if is_scaf:
            sc, st = res, wc
        else:
            sc, st = wc, res

        bp = BasePair(sc=sc, st=st, hp=(h, p))
        return (h, p), bp

    def _get_pot_bp(self) -> Dict[Tuple[int, int], BasePair]:
        bps = dict()
        done: Set[int] = set()
        for resindex in self.link.u.residues.resindices:
            if resindex in done:
                continue
            pos, bp = self._get_bp(resindex)
            bps[pos] = bp

            done.add(resindex)
            if resindex in self.link.Fbp_full:
                done.add(self.link.Fbp_full[resindex])
        return bps

    def _get_n_bp(self, bp: BasePair, steps: int = 1, local=True,
                  ) -> Optional[BasePair]:
        if steps == 0:
            return bp
        helix, position = bp.hp
        if (helix % 2) == 1 and local:
            steps = -steps
        direct = np.sign(steps)
        n_position = position + steps

        # first check the number of skips passed
        n_skips = 0
        for n in range(direct, direct * (steps + 1), direct):
            n_position = position + n
            if (helix, n_position) in self.link.Dhp_skips:
                n_skips += 1
        # move one position further if on skip
        n_position = position + direct * (steps + n_skips)
        if (helix, n_position) in self.link.Dhp_skips:
            n_position += direct

        return self.bps.get((helix, n_position), None)

    def eval_bp(self) -> None:
        """ Affects
            -------
                self.bp_quality
                self.bp_geometry_local
                self.bp_geometry_global
        """
        for bp in self.bps.values():
            if bp.sc is None or bp.st is None:
                continue
            bp_qual = self._get_bp_quality(bp=bp)
            for resindex in [bp.sc.resindex, bp.st.resindex]:
                self.bp_quality[resindex] = bp_qual

            # local: scaffold 5'->3'
            n_bp = self._get_n_bp(bp=bp, local=True)
            if n_bp is not None:
                if n_bp.is_ds:
                    bp_geom = self._get_bp_geometry(bp=bp, n_bp=n_bp)
                    for resindex in [bp.sc.resindex, bp.st.resindex]:
                        self.bp_geometry_local[resindex] = bp_geom

            # global: even helix scaffold 5'->3', odd helix scaffold 3'->5'
            n_bp = self._get_n_bp(bp=bp, local=False)
            if n_bp is not None:
                if n_bp.is_ds:
                    bp_geom = self._get_bp_geometry(bp=bp, n_bp=n_bp)
                    for resindex in [bp.sc.resindex, bp.st.resindex]:
                        self.bp_geometry_global[resindex] = bp_geom

    def _get_bp_quality(self, bp: BasePair) -> Dict[str, Any]:
        quality = dict()
        # TODO: cleanup (loop?)
        atom_name = "C1'"
        bnd = atom_name * 2
        C1p = []
        atoms = []
        for b in [bp.sc, bp.st]:
            C1p.append(b.atoms.select_atoms("name {}".format(atom_name))[0])
            atoms.append(
                b.atoms.select_atoms(
                    "name " + ' or name '.join(map(str, WC_HBONDS[b.resname]))
                )
            )

        c1p_dist = mdamath.norm(C1p[0].position - C1p[1].position)
        c1p_dev = abs(c1p_dist - C1P_BASEDIST) / C1P_BASEDIST
        quality[bnd] = c1p_dev

        # TODO: review hbond-dev, move to function
        for idx, a1 in enumerate(atoms[0]):
            hbond_dist = mdamath.norm(atoms[1][idx].position - a1.position)
            hbond_should = WC_HBONDS_DIST[bp.sc.resname][idx]
            hbond_dev = abs(hbond_dist - hbond_should) / hbond_should
            bnd = WC_HBONDS[bp.sc.resname][idx] + WC_HBONDS[bp.st.resname][idx]
            quality[bnd] = (hbond_dev if hbond_dist > hbond_should
                            else -hbond_dev)

        return quality

    def _get_bp_geometry(self, bp: BasePair, n_bp: BasePair,
                         ) -> Dict[str, Any]:
        """ Returns
            -------
                bp_quality
        """
        def _get_twist(bp: BasePair, n_bp: BasePair) -> Dict[str, float]:
            twist = dict()
            for key, a in bp.plane.a.items():
                n_a = n_bp.plane.a[key]
                twist[key] = _save_arccos_deg(_proj(a, n_a))
            return twist

        def _get_rise(bp: BasePair, n_bp: BasePair) -> Dict[str, float]:
            rise = dict()
            n0 = bp.plane.n0
            for key, P in bp.plane.P.items():
                n_P = n_bp.plane.P[key]
                rise[key] = np.abs(np.inner((n_P - P), n0))
            return rise

        def _get_shift(bp: BasePair, n_bp: BasePair) -> Dict[str, float]:
            shift = dict()
            n0 = bp.plane.n0
            for key, P, in bp.plane.P.items():
                a = bp.plane.a[key]
                n_P = n_bp.plane.P[key]
                m0 = _norm(np.cross(n0, a))
                shift[key] = np.inner((n_P - P), m0)
            return shift

        def _get_slide(bp: BasePair, n_bp: BasePair) -> Dict[str, float]:
            slide = dict()
            for key, P in bp.plane.P.items():
                a = bp.plane.a[key]
                n_P = n_bp.plane.P[key]
                slide[key] = np.inner((n_P - P), _norm(a))
            return slide

        def _get_tilt(bp: BasePair, n_bp: BasePair) -> Dict[str, float]:
            tilt = dict()
            n0 = bp.plane.n0
            n_n0 = n_bp.plane.n0
            for key, a in bp.plane.a.items():
                rot_axis = np.cross(a, n0)
                projn0 = n0 - _v_proj(n0, rot_axis)
                projn_n0 = n_n0 - _v_proj(n_n0, rot_axis)
                proj = _proj(projn0, projn_n0)
                tilt[key] = _save_arccos_deg(proj)
            return tilt

        def _get_roll(bp: BasePair, n_bp: BasePair) -> Dict[str, float]:
            roll = dict()
            n0 = bp.plane.n0
            n_n0 = n_bp.plane.n0
            for key, rot_axis in bp.plane.a.items():
                projn0 = n0 - _v_proj(n0, rot_axis)
                projn_n0 = n_n0 - _v_proj(n_n0, rot_axis)
                proj = _proj(projn0, projn_n0)
                roll[key] = _save_arccos_deg(proj)
            return roll

        geometry = {
            "rise": _get_rise(bp=bp, n_bp=n_bp),
            "slide": _get_slide(bp=bp, n_bp=n_bp),
            "shift": _get_shift(bp=bp, n_bp=n_bp),
            "twist": _get_twist(bp=bp, n_bp=n_bp),
            "tilt": _get_tilt(bp=bp, n_bp=n_bp),
            "roll": _get_roll(bp=bp, n_bp=n_bp),
        }
        return geometry

    def eval_dh(self) -> None:
        """ Affects
            -------
                self.dh_quality
        """
        for res in self.link.u.residues:
            self.dh_quality[res.resindex] = self._get_dihedrals(res)

    # TODO: improve
    def _get_dihedrals(self, res: "mda.residue") -> Dict[str, float]:
        def _get_residue_BB(res: "mda.residue"
                            ) -> Tuple[Dict[str, "np.ndarray"],
                                       Tuple[bool, bool, bool]
                                       ]:
            iniSeg, terSeg, ter5 = False, False, False

            atoms = dict()
            try:
                P = res.atoms.select_atoms("name P")[0]
                atoms["P"] = P.position
            except (KeyError, IndexError):
                ter5 = True

            for x in BB_ATOMS[:-1]:
                atom_name = "name {}".format(x)
                atoms[x] = res.atoms.select_atoms(atom_name)[0].position
            if res.resname in ["ADE", "GUA"]:
                Y = PUR_ATOMS
            else:
                Y = PYR_ATOMS
            for y in Y:
                atom_name = "name {}".format(y)
                atoms[y] = res.atoms.select_atoms(atom_name)[0].position

            try:
                n_res = self.link.u.residues[res.resindex + 1]
                if res.segindex == n_res.segindex:
                    n_P = n_res.atoms.select_atoms("name P")[0]
                    n_O5p = n_res.atoms.select_atoms("name O5'")[0]
                    atoms["P +"] = n_P.position
                    atoms["O5' +"] = n_O5p.position
                else:
                    terSeg = True
            except (KeyError, IndexError):
                terSeg = True

            try:
                p_res = self.link.u.residues[res.resindex - 1]
                if res.segindex == p_res.segindex:
                    p_O3p = p_res.atoms.select_atoms("name O3'")[0]
                    atoms["O3' -"] = p_O3p.position
                else:
                    iniSeg = True
            except (KeyError, IndexError):
                iniSeg = True

            return atoms, (ter5, terSeg, iniSeg)

        def _get_valid_DH(ter5: bool, terSeg: bool, iniSeg: bool) -> List[str]:
            dh_valid = ["gamma", "delta", "xi"]
            if terSeg is False:
                dh_valid.extend(["epsilon", "zeta"])
            if iniSeg is False:
                dh_valid.append("alpha")
            if ter5 is False:
                dh_valid.append("beta")
            return dh_valid

        def _get_dh_for_res(atoms: Dict[str, "np.ndarray"],
                            pyr: bool,
                            dh_valid: List[str],
                            ) -> Dict[str, float]:
            dh = {}
            for dh_name in DH_ATOMS:
                if dh_name in dh_valid:
                    angle = _get_dhangle(atoms, pyr, dh_name)
                else:
                    angle = np.nan
                dh[dh_name] = angle
            return dh

        def _get_dhangle(atoms: Dict[str, "np.ndarray"],
                         pyr: bool,
                         dh_name: str,
                         ) -> float:
            p = []
            if dh_name == "xi":
                if pyr:
                    for i in DH_ATOMS[dh_name]["pyr"]:
                        p.append(atoms[i])
                else:
                    for i in DH_ATOMS[dh_name]["pur"]:
                        p.append(atoms[i])
            else:
                for i in DH_ATOMS[dh_name]:
                    p.append(atoms[i])
            return _dh_angle(p)

        atoms, logic = _get_residue_BB(res)
        dh_valid = _get_valid_DH(*logic)
        if res.resname in ["ADE", "GUA"]:
            pyr = False
        else:
            pyr = True

        dh = _get_dh_for_res(atoms, pyr, dh_valid)
        return dh

    def eval_distances(self) -> None:
        """ Affects
            -------
                self.distances
        """
        for bp in self.bps.values():
            if bp.sc is not None:
                self.distances[bp.sc.resindex] = dict()
            if bp.st is not None:
                self.distances[bp.st.resindex] = dict()
            for atom in ["C1'", "P"]:
                sc_dist, st_dist = self._get_distance(bp=bp, name=atom)
                if bp.sc is not None:
                    self.distances[bp.sc.resindex][atom] = sc_dist
                if bp.st is not None:
                    self.distances[bp.st.resindex][atom] = st_dist

    def _get_distance(self, bp: BasePair, name: str
                      ) -> Tuple[Dict[str, float], Dict[str, float]]:
        def try_pos(res: "mda.residue",
                    name: str,
                    ) -> Optional["np.ndarray"]:
            try:
                return res.atoms.select_atoms(name)[0].position
            except (AttributeError, IndexError):
                return None

        sc_dist: Dict[str, float] = dict()
        st_dist: Dict[str, float] = dict()
        atom_name = "name {}".format(name)

        n_bp = self._get_n_bp(bp=bp, steps=1)
        p_bp = self._get_n_bp(bp=bp, steps=-1)
        for s, x_bp, dist in [("sc", n_bp, sc_dist), ("st", p_bp, st_dist)]:
            if x_bp is None:
                X = try_pos(res=bp.sc, name=atom_name)
                Y = try_pos(res=bp.st, name=atom_name)
                if X is not None and Y is not None:
                    dist["pair"] = np.linalg.norm(Y - X)
                    dist["stack"] = np.nan
                    dist["crossstack"] = np.nan
                else:
                    for typ in ["pair", "stack", "crossstack"]:
                        dist[typ] = np.nan
            else:
                SC = [try_pos(res=x.sc, name=atom_name) for x in [bp, x_bp]]
                ST = [try_pos(res=x.st, name=atom_name) for x in [bp, x_bp]]
                if s == "st":
                    SC, ST = ST, SC
                for X, Y, typ in [(SC[0], ST[0], "pair"),
                                  (SC[0], SC[1], "stack"),
                                  (SC[0], ST[1], "crossstack")
                                  ]:
                    if X is not None and Y is not None:
                        dist[typ] = np.linalg.norm(Y - X)
                    else:
                        dist[typ] = np.nan
        return sc_dist, st_dist

    def eval_co_angles(self) -> None:
        """ Definition: Bai, X. (2012).  doi: 10.1073/pnas.1215713109
            Affects
            -------
                self.co_angles
        """
        def get_co_angles(co: Crossover, typ="C6C8"):
            def P_from_p(p: Optional[BasePair], typ: str) -> "np.ndarray":
                if p is None:
                    return np.zeros(3)
                if p.plane is not None:
                    return p.plane.P[typ]
                elif p.sc is not None:
                    return p.sc_plane.P[typ]
                else:
                    return p.st_plane.P[typ]

            resindices = list()
            for P in co.Ps:
                if P is None:
                    continue
                elif P.sc is not None:
                    resindices.append(P.sc.resindex)
                elif P.st is not None:
                    resindices.append(P.st.resindex)

            X = [P_from_p(p, typ) for p in co.Ps]
            X_ = [P_from_p(l, typ) for l in co.Ls]

            center = sum(X) / len(X)
            x = [(p_ - p) for p, p_ in zip(X, X_)]
            n0 = _norm(X[0] + X[1] - X[2] - X[3])
            proj = _proj2plane(x, n0)

            d1 = _proj(proj[0], proj[2])
            a_n0 = _save_arccos_deg(_proj(x[0], n0))

            if co.typ != "end":
                b_n0 = _save_arccos_deg(_proj(x[1], n0))
                d2 = _proj(proj[3], proj[1])
                a1 = _proj(proj[0], [-x for x in proj[1]])
                a2 = _proj(proj[2], [-x for x in proj[3]])

                gamma1 = _save_arccos_deg(d1)
                gamma2 = _save_arccos_deg(d2)
                alpha1 = _save_arccos_deg(a1)
                alpha2 = _save_arccos_deg(a2)
            else:
                b_n0 = 90.
                gamma2 = np.nan
                alpha1 = np.nan
                alpha2 = np.nan

            gamma1 = np.rad2deg(np.arccos(d1))
            beta = 180. - abs(a_n0) - abs(b_n0)

            return {"angles": {"co_beta": beta,
                               "co_gamma1": gamma1,
                               "co_gamma2": gamma2,
                               "co_alpha1": alpha1,
                               "co_alpha2": alpha2,
                               },
                    "center-co": center,
                    "plane": n0,
                    "resindices": resindices,
                    }

        for key, co in self.link.Fco.items():
            co_data = get_co_angles(co=co)

            self.co_angles[key] = {
                "co": key,
                "type": co.typ,
                "is_scaffold": co.is_scaf,
                "angles": co_data["angles"],
                "center-co": co_data["center-co"],
                "resindices": co_data["resindices"],
            }

    def _mrc_localres(self, path_in: str) -> Dict[int, float]:
        def get_localres(atoms: "mda.atomgroup",
                         data: np.ndarray,
                         origin: np.ndarray,
                         voxel_size: np.ndarray,
                         ) -> float:
            locres = 0.
            for atom in atoms:
                atom_voxel = np.rint((atom.position - origin) / voxel_size)
                locres += data[tuple(atom_voxel.astype(int))]
            locres /= len(atoms)
            return locres

        self.link.u.trajectory[-1]
        with mrc.open(path_in, mode='r') as mrc_in:
            o = np.array(mrc_in.header["origin"])
            origin = np.array([o["x"], o["y"], o["z"]])
            c = np.array(mrc_in.header["cella"])
            cellA = np.array([c["x"], c["y"], c["z"]])
            shape = np.array([mrc_in.header["nx"],
                              mrc_in.header["ny"],
                              mrc_in.header["nz"],
                              ])
            voxel_size = cellA / shape
            data = np.swapaxes(mrc_in.data, 0, 2)

        dict_localres = {}
        for res in self.link.u.atoms.residues:
            localres = get_localres(
                atoms=res.atoms,
                data=data,
                origin=origin,
                voxel_size=voxel_size,
            )
            dict_localres[res.resindex] = localres
        return dict_localres
