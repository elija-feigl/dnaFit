#!/usr/bin/env python3
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib import mdamath

import attr
from typing import Dict, Tuple, Any, Optional

from linker import Linkage
from basepair import BasePair
from utils import (C1P_BASEDIST, WC_HBONDS, WC_HBONDS_DIST, BB_ATOMS,
                   PUR_ATOMS, PYR_ATOMS, DH_ATOMS,
                   _proj, _norm, _v_proj, _save_arccos_deg,
                   _dh_angle, _proj2plane
                   )


@attr.s
class BDna(object):
    u: "mda.universe" = attr.ib()
    link: Linkage = attr.ib()

    def __attrs_post_init__(self):
        self.link.reverse()
        self.bps: Dict[Tuple[int, int], BasePair] = self._get_pot_bp()
        self.bp_quality: Dict[int, Any] = {}
        self.bp_geometry: Dict[int, Any] = {}
        self.dh_quality: Dict[int, Any] = {}
        self.distances: Dict[int, Any] = {}
        self.co_angles: Dict[int, Any] = {}

    def sample(self):
        for bp in self.bps.values():
            bp.calculate_baseplanes()

        self.eval_bp()
        self.eval_distances()
        self.eval_dh()
        self.eval_co_angles()
        return

    def _get_bp(self, resindex: int) -> Tuple[Tuple[int, int], BasePair]:
        h, p, is_scaf = self.link.DidDhps[self.link.FidDid[resindex]]
        res = self.u.residues[resindex]
        if resindex in self.link.Fbp_full:
            wcindex = self.link.Fbp_full[resindex]
            wc = self.u.residues[wcindex]
        else:
            wcindex = None
            wc = None

        if is_scaf:
            sc, st = res, wc
        else:
            sc, st = wc, res

        bp = BasePair(sc=sc, st=st, hp=(h, p))
        return (h, p), bp

    def _get_pot_bp(self) -> Dict[Tuple[int, int], BasePair]:
        bps = dict()
        done = set()
        # TODO: -low- loop over residues not resindex?
        for resindex in self.u.residues.resindices:
            if resindex in done:
                continue
            pos, bp = self._get_bp(resindex)
            bps[pos] = bp

            done.add(resindex)
            if resindex in self.link.Fbp_full.keys():
                done.add(self.link.Fbp_full[resindex])
        return bps

    def _get_n_bp(self, bp: BasePair, steps: int = 1
                  ) -> Optional[BasePair]:
        """ get next residue and its complemt wc pair whithin a helix
            check if next exists.
            check has bp (ony scaffold residues apply here)
            else returns None
        """
        if steps == 0:
            return bp
        helix, position = bp.hp
        if (helix % 2) == 1:  # scaffold 5->3: even l->r odd r->l
            steps = -steps

        stp = np.sign(steps)
        n_skips = 0

        # check if we passed skips
        for n in range(stp, steps + stp, stp):
            n_helix, n_position = (helix, position + n)
            if (n_helix, n_position, True) in self.link.Dskips:
                n_skips += 1

        # check if land on skip
        n_position = position + steps + n_skips
        while (n_helix, n_position, True) in self.link.Dskips:
            n_position += 1

        if (helix, n_position) in self.bps.keys():
            return self.bps[(helix, n_position)]
        else:
            return None

    def eval_bp(self):
        """ Affects
            -------
                self.bp_quality
                self.bp_geometry
        """
        for bp in self.bps.values():
            if bp.sc is None or bp.st is None:
                continue
            bp_qual = self._get_bp_quality(bp=bp)
            for resindex in [bp.sc.resindex, bp.st.resindex]:
                self.bp_quality[resindex] = bp_qual

            n_bp = self._get_n_bp(bp=bp)
            if n_bp is not None:
                if n_bp.is_ds:
                    bp_geom = self._get_bp_geometry(bp=bp, n_bp=n_bp)
                    for resindex in [bp.sc.resindex, bp.st.resindex]:
                        self.bp_geometry[resindex] = bp_geom

        return

    def _get_bp_quality(self, bp: BasePair) -> Tuple[dict, dict]:
        quality = {}

        # C1p distance
        # TODO: cleanup (loop?)
        atom_name = "C1'"
        bnd = atom_name * 2
        C1p = []
        atoms = []
        for b in [bp.sc, bp.st]:
            C1p.append(b.atoms.select_atoms("name {}".format(atom_name))[0])
            atoms.append(b.atoms.select_atoms(
                "name " + ' or name '.join(map(str, WC_HBONDS[b.resname]))))

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

    def _get_bp_geometry(self, bp: BasePair, n_bp: BasePair
                         ) -> Tuple[dict, dict]:
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

        geometry = {"rise": _get_rise(bp=bp, n_bp=n_bp),
                    "slide": _get_slide(bp=bp, n_bp=n_bp),
                    "shift": _get_shift(bp=bp, n_bp=n_bp),
                    "twist": _get_twist(bp=bp, n_bp=n_bp),
                    "tilt": _get_tilt(bp=bp, n_bp=n_bp),
                    "roll": _get_roll(bp=bp, n_bp=n_bp),
                    "next_seq": n_bp.sc.resname[0] + n_bp.st.resname[0]
                    }
        return geometry

    def eval_dh(self):
        """ Affects
            -------
                dh_quality
        """
        self.dh_quality = {}
        for res in self.u.residues:
            self.dh_quality[res.resindex] = self._get_dihedrals(res)
        return

    def _get_dihedrals(self, res):
        def _get_residue_BB(res):  # slow
            iniSeg, terSeg, ter5 = False, False, False

            atoms = {}
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
                n_res = self.u.residues[res.resindex + 1]
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
                p_res = self.u.residues[res.resindex - 1]
                if res.segindex == p_res.segindex:
                    p_O3p = p_res.atoms.select_atoms("name O3'")[0]
                    atoms["O3' -"] = p_O3p.position
                else:
                    iniSeg = True
            except (KeyError, IndexError):
                iniSeg = True

            return atoms, (ter5, terSeg, iniSeg)

        def _get_valid_DH(ter5, terSeg, iniSeg):
            dh_valid = ["gamma", "delta", "xi"]
            if terSeg is False:
                dh_valid.extend(["epsilon", "zeta"])
            if iniSeg is False:
                dh_valid.append("alpha")
            if ter5 is False:
                dh_valid.append("beta")
            return dh_valid

        def _get_dh_for_res(atoms, pyr, dh_valid):
            dh = {}
            for dh_name in DH_ATOMS.keys():
                if dh_name in dh_valid:
                    angle = _get_dhangle(atoms, pyr, dh_name)
                else:
                    angle = None
                dh[dh_name] = angle
            return dh

        def _get_dhangle(atoms, pyr, dh_name):  # slow
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

    def eval_distances(self):
        """ Affects
            -------
                distances
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
        return

    def _get_distance(self, bp: BasePair, name: str
                      ) -> Tuple[Dict[str, float], Dict[str, float]]:
        def try_position(res: "mda.residue", atom_name: str
                         ) -> Optional["np.ndarray"]:
            try:
                return res.atoms.select_atoms(atom_name)[0].position
            except (AttributeError, IndexError):
                return None

        sc_dist, st_dist = dict(), dict()
        atom_name = "name {}".format(name)

        n_bp = self._get_n_bp(bp)
        SC, ST = list(), list()
        if n_bp is None:
            X = try_position(res=bp.sc, atom_name=atom_name)
            Y = try_position(res=bp.st, atom_name=atom_name)
            if X is not None and Y is not None:
                for dist in [sc_dist, st_dist]:
                    dist["pair"] = np.linalg.norm(Y - X)
                    dist["stack"], dist["crossstack"] = None, None
            else:
                for dist in [sc_dist, st_dist]:
                    for typ in ["pair", "stack", "crossstack"]:
                        dist[typ] = None
        else:
            for x in [bp, n_bp]:
                SC.append(try_position(res=x.sc, atom_name=atom_name))
                ST.append(try_position(res=x.st, atom_name=atom_name))

            for A, B, dist in [(SC, ST, sc_dist), (ST, SC, st_dist)]:
                for X, Y, typ in [(A[0], B[0], "pair"),
                                  (A[0], A[1], "stack"),
                                  (A[0], B[1], "crossstack")
                                  ]:
                    if X is not None and Y is not None:
                        dist[typ] = np.linalg.norm(Y - X)
                    else:
                        dist[typ] = np.nan

        return sc_dist, st_dist

    def eval_co_angles(self):
        """ Definition: Bai, X. (2012).  doi: 10.1073/pnas.1215713109
            Affects
            -------
                co_angles
        """

        def get_co_angles_end(a, a_leg, c, c_leg, typ="C6C8"):
            
            if a.plane is not None:
                A = a.plane.P[typ]
            elif a.sc is not None:
                A = a.sc_plane.P[typ]
            else:
                A = a.st_plane.P[typ]

            if a_leg.plane is not None:
                A_ = a_leg.plane.P[typ]
            elif a.sc is not None:
                A_ = a_leg.sc_plane.P[typ]
            else:
                A_ = a_leg.st_plane.P[typ]

            if c.plane is not None:
                C = c.plane.P[typ]
            elif a.sc is not None:
                C = c.sc_plane.P[typ]
            else:
                C = c.st_plane.P[typ]

            if c_leg.plane is not None:
                C_ = c_leg.plane.P[typ]
            elif a.sc is not None:
                C_ = c_leg.sc_plane.P[typ]
            else:
                C_ = c_leg.st_plane.P[typ]

            center = (A + C) * 0.5
            a = A_ - A
            c = C_ - C
            n0 = _norm((A - C))

            proj_ac = _proj2plane([a, c], n0)
            dist = _proj(proj_ac[0], proj_ac[1])
            gamma = np.rad2deg(np.arccos(dist))
            beta = 90. - _save_arccos_deg(_proj(a, n0))

            return {"angles": {"co_beta": beta,
                               "co_gamma": gamma
                               },
                    "center-co": center,
                    "plane": n0,
                    }

        # TODO: -low- cleanup
        def get_co_angles_full(bpplanes, double_bpplanes, typ="C6C8"):
            points = []

            for x in [bpplanes[0:2], double_bpplanes[0:2], bpplanes[2:],
                      double_bpplanes[2:]]:  # abcd
                ins = x[0].plane.P[typ]
                out = x[1].plane.P[typ]
                points.append((ins, out))

            # abcd
            center = sum(p[0] for p in points) * 0.25
            n0 = _norm(points[0][0] + points[1][0] - points[2][0] - points[3][0])

            abcd = []
            for p_in, p_out in points:
                abcd.append(p_out - p_in)

            proj_abcd = _proj2plane(abcd, n0)

            d1 = _proj(proj_abcd[0], proj_abcd[2])
            gamma1 = _save_arccos_deg(d1)
            d2 = _proj(proj_abcd[3], proj_abcd[1])
            gamma2 = _save_arccos_deg(d2)
            a1 = _proj(proj_abcd[0], [-x for x in proj_abcd[1]])
            alpha1 = _save_arccos_deg(a1)
            a2 = _proj(proj_abcd[2], [-x for x in proj_abcd[3]])
            alpha2 = _save_arccos_deg(a2)

            a_n0 = _save_arccos_deg(_proj(abcd[0], n0))
            b_n0 = _save_arccos_deg(_proj(abcd[1], n0))
            beta = 180. - abs(a_n0) - abs(b_n0)

            return {"angles": {"co_beta": beta,
                               "co_gamma1": gamma1,
                               "co_gamma2": gamma2,
                               "co_alpha1": alpha1,
                               "co_alpha2": alpha2},
                    "center-co": center,
                    "plane": n0,
                    }

        co_done = set()
        for a_index, co in self.link.Fco.items():
            if a_index not in co_done:
                a_leg_index = co["leg"]
                c_index = co["co"]
                c_leg_index = self.link.Fco[c_index]["leg"]

                co_done.update([a_index, c_index])

                a = self.bps[self.link.DidDhps[self.link.FidDid[a_index]][:2]]
                a_leg = self.bps[self.link.DidDhps[self.link.FidDid[a_leg_index]][:2]]
                c = self.bps[self.link.DidDhps[self.link.FidDid[c_index]][:2]]
                c_leg = self.bps[self.link.DidDhps[self.link.FidDid[c_leg_index]][:2]]

                co_type = co["type"][0]
                co_type2 = self.link.Fco[a_index]["type"][0]
                is_full = (co_type == "double")
                is_half = (co_type == "single" and co_type2 != "end")

                if is_full:
                    b_index = co["type"][1]
                    b_leg_index = self.link.Fco[b_index]["leg"]
                    d_index = self.link.Fco[b_index]["co"]
                    d_leg_index = self.link.Fco[b_index]["leg"]
                    co_done.update([b_index, d_index])
                    b = self.bps[self.link.DidDhps[self.link.FidDid[b_index]][:2]]
                    b_leg = self.bps[self.link.DidDhps[self.link.FidDid[b_leg_index]][:2]]
                    d = self.bps[self.link.DidDhps[self.link.FidDid[d_index]][:2]]
                    d_leg = self.bps[self.link.DidDhps[self.link.FidDid[d_leg_index]][:2]]

                    co_data = get_co_angles_full([a, a_leg, c, c_leg],
                                                 [b, b_leg, d, d_leg])  # ac bd
                    crossover_ids = (a_index, b_index, c_index, d_index)
                elif is_half:
                    b_index = co["type"][1]
                    b_leg_index = co["type"][2]
                    d_index = self.link.Fco[c_index]["type"][1]
                    d_leg_index = self.link.Fco[c_index]["type"][2]
                    co_done.update([b_index, d_index])

                    b = self.bps[self.link.DidDhps[self.link.FidDid[b_index]][:2]]
                    b_leg = self.bps[self.link.DidDhps[self.link.FidDid[b_leg_index]][:2]]
                    d = self.bps[self.link.DidDhps[self.link.FidDid[d_index]][:2]]
                    d_leg = self.bps[self.link.DidDhps[self.link.FidDid[d_leg_index]][:2]]

                    co_data = get_co_angles_full([a, a_leg, c, c_leg],
                                                 [b, b_leg, d, d_leg])  # ac bd
                    crossover_ids = (a_index, b_index, c_index, d_index)
                else:
                    co_data = get_co_angles_end(a, a_leg, c, c_leg)
                    crossover_ids = (a_index, c_index)

                self.co_angles[co["co_index"]] = {
                    "ids(abcd)": crossover_ids, "type": co_type,
                    "is_scaffold": co["is_scaffold"],
                    "angles": co_data["angles"],
                    "center-co": co_data["center-co"]}
        return
