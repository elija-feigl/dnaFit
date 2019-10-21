#!/usr/bin/env python3
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib import mdamath

import attr
from typing import Dict, Tuple, Any, Optional

from linker import Linkage
from utils import (C1P_BASEDIST, WC_HBONDS, WC_HBONDS_DIST, TOL, BB_ATOMS,
                   PUR_ATOMS, PYR_ATOMS, DH_ATOMS,
                   _proj, _norm, _v_proj, _save_arccos_deg, _dh_angle
                   )


@attr.s(slots=True, frozen=True)
class BasePairPlane(object):
    P: Dict[str, "np.ndarray"] = attr.ib()
    a: Dict[str, "np.ndarray"] = attr.ib()
    n0: "np.ndarray" = attr.ib()


@attr.s(slots=True, frozen=True)
class BasePlane(object):
    P: Dict[str, "np.ndarray"] = attr.ib()
    n0: "np.ndarray" = attr.ib()


@attr.s
class BasePair(object):
    """ every square of the JSON can be represented as BP
    """
    sc: "mda.Residue" = attr.ib()
    st: "mda.Residue" = attr.ib()
    hp: Tuple[int, int] = attr.ib()

    def __attrs_post_init__(self):
        if self.sc is None or self.st is None:
            self.is_ds = False
        else:
            self.is_ds = True
        self.calculate_baseplanes()

    def calculate_baseplanes(self):
        self.sc_plane = (self._get_base_plane(self.sc)
                         if self.sc is not None else None)
        self.st_plane = (self._get_base_plane(self.st)
                         if self.st is not None else None)
        self.plane = (self._get_bp_plane(sc=self.sc_plane, st=self.st_plane)
                      if self.is_ds else None)
        return

    def _get_base_plane(self, res: "mda.Residue") -> BasePlane:
        P = dict()
        atom = []
        for atom_name in ["C2", "C4", "C6"]:
            A = res.atoms.select_atoms("name " + atom_name)[0]
            atom.append(A.position)

        n0 = _norm(np.cross((atom[1] - atom[0]), (atom[2] - atom[1])))
        P["diazine"] = sum(atom) / 3.

        C6C8 = "C8" if res.resname in ["ADE", "GUA"] else "C6"
        P["C6C8"] = res.atoms.select_atoms("name " + C6C8)[0].position

        P["C1'"] = res.atoms.select_atoms("name C1'")[0].position

        return BasePlane(n0=n0, P=P)

    def _get_bp_plane(self, sc, st) -> BasePairPlane:
        a, P = dict(), dict()
        n0 = (sc.n0 + st.n0) * 0.5
        for n in sc.P.keys():
            P[n] = (sc.P[n] + st.P[n]) * 0.5
            a[n] = sc.P[n] - st.P[n]

        return BasePairPlane(n0=n0, a=a, P=P)


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

    def _get_pot_bp(self) -> Dict[Tuple[int, int], BasePair]:
        bps = dict()
        done = set()
        # TODO: -low- loop over residues not resindex?
        for sc_resindex in self.u.residues.resindices:
            if sc_resindex in done:
                continue

            sc = self.u.residues[sc_resindex]
            h, p, is_scaf = self.link.DidDhps[self.link.FidDid[sc_resindex]]

            if sc_resindex in self.link.Fbp_full:
                st_resindex = self.link.Fbp_full[sc_resindex]
                st = self.u.residues[st_resindex]
            else:
                st_resindex = None
                st = None

            if not is_scaf:
                sc, st = st, sc
                sc_resindex, st_resindex = st_resindex, sc_resindex

            bps[(h, p)] = BasePair(sc=sc, st=st, hp=(h, p))
            done.update([sc_resindex, st_resindex])
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
        def _get_residue_BB(res):
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

        def _get_dhangle(atoms, pyr, dh_name):
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
        """ Definition of the angles enclosed by the four helical legs of a
            cross-over.
            Vectors are computed using the coordinates of base pair midpoints
            at the cross-over position and 2 bp away from  cross-over in each
            leg.
            The cross-over vector x is computed from the coordinates of the
            midpoints between the two base pairs in each of the two helices at
            the crossover position and is normal to what we call the cross-over
            plane. The subscript “jj” indicates vectorial projections into the
            cross-over plane. The angle β is also computed as indicated for
            vectors C and D.(Bai 2012)
        """

        def get_co_baseplanes(res_index, leg_index, co_index, coleg_index):
            bases = []
            for r in [res_index, leg_index, co_index, coleg_index]:
                try:
                    wc_r = self.link.Fbp_full[r]
                    bases.append(self._get_base_plane(self.u.residues[r]))
                    bases.append(self._get_base_plane(self.u.residues[wc_r]))
                except KeyError:
                    bases.append(self._get_base_plane(self.u.residues[r]))
                    bases.append(self._get_base_plane(self.u.residues[r]))
            bp_planes = []
            for i in [0, 2, 4, 6]:
                bp_planes.append(
                    {"center-C1p": ((bases[i]["C1p"] +
                                     bases[i + 1]["C1p"]) * 0.5),
                     "dir-C1p": (bases[i]["C1p"] - bases[i + 1]["C1p"]),
                     "center-diazine": ((bases[i]["diazine"] +
                                         bases[i + 1]["diazine"]) * 0.5),
                     "dir-diazine": (bases[i]["diazine"] -
                                     bases[i + 1]["diazine"]),
                     "center-C6C8": ((bases[i]["C6C8"] +
                                      bases[i + 1]["C6C8"]) * 0.5),
                     "dir-C6C8": (bases[i]["C6C8"] -
                                  bases[i + 1]["C6C8"]),
                     "n0": ((bases[i]["n0"] + bases[i + 1]["n0"]) * 0.5)})
            return bp_planes

        def get_co_angles_end(bpplanes):  # TODO: -low- cleanup
            """def project_to_plane(vect, n0):
                pro = []
                for x in vect:
                    d = np.inner(x, n0)
                    p = [d * n0[i] for i in range(len(n0))]
                    pro.append([x[i] - p[i] for i in range(len(x))])
                return pro"""

            a1 = bpplanes[0]["center-C6C8"]
            a2 = bpplanes[1]["center-C6C8"]
            c1 = bpplanes[2]["center-C6C8"]
            # c2 = bpplanes[3]["center-C6C8"]

            center = (a1 + c1) * 0.5
            a = a2 - a1
            # c = c2 - c1

            n0 = _norm((a1 - c1))
            # proj_ac = project_to_plane([a, c], n0)

            # dist = _proj(proj_ac[0], proj_ac[0])
            # if 1.0 < abs(dist) < 1.0 + TOL : dist = np.sign(dist)
            # gamma1 = np.rad2deg(np.arccos(dist))
            ang_temp = np.rad2deg(np.arccos(_proj(a, n0)))  # unprojected
            beta = 90. - ang_temp

            return {"angles": {"co_beta": beta}, "center-co": center,
                    "plane": n0}

        # TODO: -low- cleanup
        def get_co_angles_full(bpplanes, double_bpplanes):

            def project_to_plane(vect, n0):
                pro = []
                for x in vect:
                    d = np.inner(x, n0)
                    p = [d * n0[i] for i in range(len(n0))]
                    pro.append([x[i] - p[i] for i in range(len(x))])
                return pro

            points = []

            for x in [bpplanes[0:2], double_bpplanes[0:2], bpplanes[2:],
                      double_bpplanes[2:]]:  # abcd
                ins = x[0]["center-C6C8"]
                out = x[1]["center-C6C8"]
                points.append((ins, out))

            # abcd
            # ((a0 +b0)/2 +  (c0 +d0)/2)/2
            center = sum(p[0] for p in points) * 0.25
            # n0((a0 +b0)/2 -  (c0 +d0)/2)
            n0 = _norm(points[0][0] + points[1][0] -
                       points[2][0] - points[3][0])

            abcd = []
            for p_in, p_out in points:
                abcd.append(p_out - p_in)

            proj_abcd = project_to_plane(abcd, n0)

            # gamma1 = a|| c||
            d1 = _proj(proj_abcd[0], proj_abcd[2])
            if 1.0 < abs(d1) < 1.0 + TOL:
                d1 = np.sign(d1)
            gamma1 = np.rad2deg(np.arccos(d1))
            # gamma2 = d|| b||
            d2 = _proj(proj_abcd[3], proj_abcd[1])
            if 1.0 < abs(d2) < 1.0 + TOL:
                d2 = np.sign(d2)
            gamma2 = np.rad2deg(np.arccos(d2))

            # alpha1 = a|| -b||
            a1 = _proj(proj_abcd[0], [-x for x in proj_abcd[1]])
            if 1.0 < abs(a1) < 1.0 + TOL:
                a1 = np.sign(a1)
            alpha1 = np.rad2deg(np.arccos(a1))
            # alpha2 = c|| -d||
            a2 = _proj(proj_abcd[2], [-x for x in proj_abcd[3]])
            if 1.0 < abs(a2) < 1.0 + TOL:
                a2 = np.sign(a2)
            alpha2 = np.rad2deg(np.arccos(a2))

            # 180 - a n0 - b n0 , unprojected
            ang_temp1 = np.rad2deg(np.arccos(_proj(abcd[0], n0)))
            if ang_temp1 > 90.:
                ang_temp1 = 180. - ang_temp1
            ang_temp2 = np.rad2deg(np.arccos(_proj(abcd[1], n0)))
            if ang_temp2 > 90.:
                ang_temp2 = 180. - ang_temp2

            beta = 180. - ang_temp1 - ang_temp2

            return {"angles": {"co_beta": beta, "co_gamma1": gamma1,
                               "co_gamma2": gamma2, "co_alpha1": alpha1,
                               "co_alpha2": alpha2},
                    "center-co": center, "plane": n0}

        self.co_angles = {}

        co_done = set()
        for res_index, co in self.link.Fco.items():
            if res_index not in co_done:
                leg_index = co["leg"]
                co_index = co["co"]
                coleg_index = self.link.Fco[co_index]["leg"]

                co_done.update([res_index, co_index])
                # res -> a, co -> c
                bpplanes = get_co_baseplanes(
                    res_index, leg_index, co_index, coleg_index)

                co_type = co["type"][0]
                if co_type == "double":
                    double_res_index = co["type"][1]
                    double_leg_index = self.link.Fco[double_res_index]["leg"]
                    double_co_index = self.link.Fco[double_res_index]["co"]
                    double_coleg_index = self.link.Fco[double_co_index]["leg"]
                    co_done.update([double_res_index, double_co_index])
                    # double -> b, d_co -> d
                    double_bpplanes = get_co_baseplanes(
                        double_res_index, double_leg_index, double_co_index,
                        double_coleg_index)
                    co_data = get_co_angles_full(
                        bpplanes, double_bpplanes)  # ac bd
                    crossover_ids = (res_index, double_res_index,
                                     co_index, double_co_index)
                elif (co_type == "single" and
                        self.link.Fco[co_index]["type"][0] != "end"):
                    single_res_index = co["type"][1]
                    single_leg_index = co["type"][2]
                    single_co_index = self.link.Fco[co_index]["type"][1]
                    single_coleg_index = self.link.Fco[co_index]["type"][2]
                    co_done.update([single_res_index, single_co_index])
                    # single -> b, s_co -> d
                    single_bpplanes = get_co_baseplanes(
                        single_res_index, single_leg_index, single_co_index,
                        single_coleg_index)
                    co_data = get_co_angles_full(
                        bpplanes, single_bpplanes)  # ac bd
                    crossover_ids = (res_index, single_res_index,
                                     co_index, single_co_index)
                else:
                    co_data = get_co_angles_end(bpplanes)
                    crossover_ids = (res_index, co_index)

                self.co_angles[co["co_index"]] = {
                    "ids(abcd)": crossover_ids, "type": co_type,
                    "is_scaffold": co["is_scaffold"],
                    "angles": co_data["angles"],
                    "center-co": co_data["center-co"]}
        return

    # TODO delete
    def _get_base_plane(self, res):

        atom = []
        for atom_name in ["C2", "C4", "C6", "C1'"]:
            A = res.atoms.select_atoms("name {}".format(atom_name))[0]
            atom.append(A.position)

        n0 = _norm(np.cross((atom[1] - atom[0]), (atom[2] - atom[1])))
        diazine_center = sum(atom[:-1]) / 3.
        if res.resname in ["ADE", "GUA"]:
            ref = res.atoms.select_atoms("name C8")[0].position
        else:
            ref = res.atoms.select_atoms("name C6")[0].position

        plane = {"n0": n0, "C1p": atom[-1], "diazine": diazine_center,
                 "C6C8": ref}

        return plane