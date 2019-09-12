#!/usr/bin/env python3
import MDAnalysis as mda
import numpy as np
import sys
import os

from MDAnalysis.lib import mdamath

import pickle
import bpLinker
import contextlib
import argparse

from pathlib import Path
from itertools import chain
from operator import attrgetter
from typing import List, Set, Dict, Tuple, Optional
from collections import namedtuple

@contextlib.contextmanager
def ignored(*exceptions):
    try:
        yield
    except exceptions:
        pass

# TODO: -mid- DOC
# TODO: -mid- move CONSTANTS

C1P_BASEDIST = 10.7
TOL = 10e6
WC_DICT = {"DC": "DG", "DG": "DC", "DT": "DA", "DA": "DT",
           "C": "G", "G": "C", "T": "A", "A": "T",
           "CYT": "GUA", "GUA": "CYT", "THY": "ADE", "ADE": "THY"}
WC_HBONDS = {"DA": ("N6", "N1"), "A": ("N6", "N1"), "ADE": ("N6", "N1"),
             "DT": ("O4", "N3"), "T": ("O4", "N3"), "THY": ("O4", "N3"),
             "DG": ("O6", "N1", "N2"), "G": ("O6", "N1", "N2"),
             "GUA": ("O6", "N1", "N2"), "DC": ("N4", "N3", "O2"),
             "C": ("N4", "N3", "O2"), "CYT": ("N4", "N3", "O2")}
(N6O4, N1N3_AT) = (2.95, 2.88)
(O6N4, N1N3_GC, N2O2) = (2.85, 2.85, 2.85)
WC_HBONDS_DIST = {"DA": (N6O4, N1N3_AT), "A": (N6O4, N1N3_AT),
                  "ADE": (N6O4, N1N3_AT), "DT": (N6O4, N1N3_AT),
                  "T": (N6O4, N1N3_AT), "THY": (N6O4, N1N3_AT),
                  "DG": (O6N4, N1N3_GC, N2O2), "G": (O6N4, N1N3_GC, N2O2),
                  "GUA": (O6N4, N1N3_GC, N2O2), "DC": (O6N4, N1N3_GC, N2O2),
                  "C": (O6N4, N1N3_GC, N2O2), "CYT": (O6N4, N1N3_GC, N2O2)}
DH_ATOMS = {"alpha": ("O3' -", "P", "O5'", "C5'"),
            "beta": ("P", "O5'", "C5'", "C4'"),
            "gamma": ("O5'", "C5'", "C4'", "C3'"),
            "delta": ("C5'", "C4'", "C3'", "O3'"),
            "epsilon": ("C4'", "C3'", "O3'", "P +"),
            "zeta": ("C3'", "O3'", "P +", "O5' +"),
            "xi": {"pyr": ("C4'", "C1'", "N1", "C2"),
                   "pur": ("C4'", "C1'", "N9", "C4")}}
BB_ATOMS = ["C1'", "O3'", "C3'", "C4'", "O5'", "C5'", "P"]
PYR_ATOMS = ["N1", "C2"]
PUR_ATOMS = ["N9", "C4"]

WC_PROPERTIES = ["rise", "slide", "shift", "twist", "tilt", "roll"]


def _norm(vector):
    return vector / np.linalg.norm(vector)


def _proj(u, v):
    return np.inner(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))


def _v_proj(u, v):
    return np.inner(u, v) / (np.linalg.norm(v) * np.linalg.norm(v)) * v


class BDna(object):

    def __init__(self, universe, d_Fbp, d_DidFid, d_DhpsDid, d_Fco, l_Dskips):
        self.u = universe
        self.d_DidFid = d_DidFid
        self.d_FidDid = {v: k for k, v in d_DidFid.items()}
        self.d_Fbp = d_Fbp
        self.d_Fbp_full = {**d_Fbp, **{v: k for k, v in d_Fbp.items()}}
        self.d_DhpsDid = d_DhpsDid
        self.d_DidDhps = {v: k for k, v in d_DhpsDid.items()}
        self.d_Fco = d_Fco
        self.l_Dskips = l_Dskips

        self.wc_quality = None
        self.wc_geometry = None
        self.dh_quality = None
        self.distances = None
        self.co_angles = None

    def _get_next_wc(self, resindex, resindex_wc, steps=1):
        """ get next residue and its complemt wc pair whithin a helix
            check if next exists.
            check has bp (ony scaffold residues apply here)
            else returns None
        """
        if steps == 0:
            return resindex, resindex_wc

        n_resindex = resindex + steps
        h, p, is_scaffold = self.d_DidDhps[self.d_FidDid[resindex]]
        n_skips = 0
        # check if we passed skips
        for n in range(1, steps+np.sign(steps), np.sign(steps)):
            if (h, p+n, is_scaffold) in self.l_Dskips:
                n_skips += np.sign(steps)
        # check if land on skip
        while (h, p+steps+n_skips, is_scaffold) in self.l_Dskips:
            n_skips += np.sign(steps)

        try:
            n_resindex = self.d_DidFid[
                self.d_DhpsDid[(h, p+steps+n_skips, is_scaffold)]]
        except KeyError:
            n_resindex = None
        try:
            self.u.residues[n_resindex]
            try:
                n_resindex_wc = self.d_Fbp_full[n_resindex]
            except KeyError:
                n_resindex_wc = None
        except IndexError:
            n_resindex = None
            n_resindex_wc = None

        return n_resindex, n_resindex_wc

    def eval_wc(self):
        """  dict of wc-quality for each residue (key = resindex)
        """
        self.wc_quality = {}
        self.wc_geometry = {}
        for resindex, resindex_wc in self.d_Fbp.items():
            res = self.u.residues[resindex]
            res_wc = self.u.residues[resindex_wc]

            self.wc_quality[resindex], self.wc_quality[resindex_wc] = (
                self._get_wc_quality(res, res_wc))

            n_resindex, n_resindex_wc = self._get_next_wc(
                resindex, resindex_wc)

            if n_resindex is not None and n_resindex_wc is not None:
                n_res = self.u.residues[n_resindex]
                n_res_wc = self.u.residues[n_resindex_wc]
                self.wc_geometry[resindex], self.wc_geometry[resindex_wc] = (
                    self._get_wc_geometry(res, res_wc, n_res, n_res_wc))

        return

    def _get_wc_quality(self, res1, res2):
        quality = {}

        # C1p distance
        c1p1 = res1.atoms.select_atoms("name C1'")[0]
        c1p2 = res2.atoms.select_atoms("name C1'")[0]

        c1p_dist = mdamath.norm(c1p1.position - c1p2.position)
        bond = "C1'C1'"
        c1p_dev = abs(c1p_dist - C1P_BASEDIST)/C1P_BASEDIST
        quality[bond] = c1p_dev

        # H-bond distances
        atoms1 = res1.atoms.select_atoms(
            "name " + ' or name '.join(map(str, WC_HBONDS[res1.resname])))
        atoms2 = res2.atoms.select_atoms(
            "name " + ' or name '.join(map(str, WC_HBONDS[res2.resname])))

        for idx, a1 in enumerate(atoms1):
            hbond_dist = mdamath.norm(atoms2[idx].position - a1.position)
            hbond_should = WC_HBONDS_DIST[res1.resname][idx]
            hbond_dev = abs(hbond_dist - hbond_should)/hbond_should
            bond = WC_HBONDS[res1.resname][idx] + WC_HBONDS[res2.resname][idx]
            quality[bond] = (hbond_dev if hbond_dist > hbond_should
                             else -hbond_dev)

        return quality, quality

    def _get_wc_geometry(self, res, res_wc, n_res, n_res_wc):
        """  dict of geometry of basepairs for each wc-pair (key= resid-scaffo)
        """
        def _get_twist(basepair, n_basepair):
            twist = []
            for direct in ["dir-C1p", "dir-diazine", "dir-C6C8"]:
                v1 = basepair[direct]
                v2 = n_basepair[direct]
                twist.append(np.rad2deg(np.arccos(_proj(v1, v2))))

            return {"C1p": twist[0], "diazine": twist[1],  "C6C8": twist[2]}

        def _get_rise(basepair, n_basepair):
            rise = []
            n0 = basepair["n0"]
            for direct in ["center-C1p", "center-diazine", "center-C6C8"]:
                P1 = basepair[direct]
                P2 = n_basepair[direct]
                rise.append(np.abs(np.inner((P2 - P1), n0)))

            return {"C1p": rise[0], "diazine": rise[1], "C6C8": rise[2]}

        def _get_shift(basepair, n_basepair):
            shift = []
            for direct in [("center-C1p", "dir-C1p"),
                           ("center-diazine", "dir-diazine"),
                           ("center-C6C8", "dir-C6C8")]:
                n0 = _norm(np.cross(basepair["n0"], basepair[direct[1]]))
                P1 = basepair[direct[0]]
                P2 = n_basepair[direct[0]]
                shift.append(np.inner((P2 - P1), n0))

            return {"C1p": shift[0], "diazine": shift[1], "C6C8": shift[2]}

        def _get_slide(basepair, n_basepair):
            slide = []
            for direct in [("center-C1p", "dir-C1p"),
                           ("center-diazine", "dir-diazine"),
                           ("center-C6C8", "dir-C6C8")]:
                n0 = _norm(basepair[direct[1]])
                P1 = basepair[direct[0]]
                P2 = n_basepair[direct[0]]
                slide.append(np.inner((P2 - P1), n0))

            return {"C1p": slide[0], "diazine": slide[1], "C6C8": slide[2]}

        def _get_tilt(basepair, n_basepair):
            tilt = []
            n0 = basepair["n0"]
            n_n0 = n_basepair["n0"]

            for direct in ["dir-C1p", "dir-diazine", "dir-C6C8"]:
                rot_axis = np.cross(basepair[direct], n0)
                projn0 = n0 - _v_proj(n0, rot_axis)
                projn_n0 = n_n0 - _v_proj(n_n0, rot_axis)

                dist = _proj(projn0, projn_n0)
                if 1.0 < abs(dist) < 1.0 + TOL:
                    dist = np.sign(dist)
                if dist > 0:
                    tilt.append(np.rad2deg(np.arccos(dist)))
                else:
                    tilt.append(np.rad2deg(- np.arccos(abs(dist))))

            for value in tilt:
                if value > 90.:
                    value - 180.

            return {"C1p": tilt[0], "diazine": tilt[1], "C6C8": tilt[2]}

        def _get_roll(basepair, n_basepair):
            roll = []
            n0 = basepair["n0"]
            n_n0 = n_basepair["n0"]

            for direct in ["dir-C1p", "dir-diazine", "dir-C6C8"]:
                rot_axis = basepair[direct]
                projn0 = n0 - _v_proj(n0, rot_axis)
                projn_n0 = n_n0 - _v_proj(n_n0, rot_axis)

                dist = _proj(projn0, projn_n0)
                if 1.0 < abs(dist) < 1.0 + TOL:
                    dist = np.sign(dist)
                if dist > 0:
                    roll.append(np.rad2deg(np.arccos(dist)))
                else:
                    roll.append(np.rad2deg(- np.arccos(abs(dist))))

            return {"C1p": roll[0], "diazine": roll[1], "C6C8": roll[2]}

        bases = []
        for r in [res, res_wc, n_res, n_res_wc]:
            bases.append(self._get_base_plane(r))

        basepairs = []
        for i in [0, 2]:
            basepairs.append(
                {"center-C1p": ((bases[i]["C1p"] +
                                 bases[i+1]["C1p"]) * 0.5),
                 "dir-C1p": (bases[i]["C1p"] - bases[i+1]["C1p"]),
                 "center-diazine": ((bases[i]["diazine"] +
                                     bases[i+1]["diazine"]) * 0.5),
                 "dir-diazine": (bases[i]["diazine"] -
                                 bases[i+1]["diazine"]),
                 "center-C6C8": ((bases[i]["C6C8"] +
                                  bases[i+1]["C6C8"]) * 0.5),
                 "dir-C6C8": (bases[i]["C6C8"] -
                              bases[i+1]["C6C8"]),
                 "n0": ((bases[i]["n0"] + bases[i+1]["n0"]) * 0.5)})

        geometry = {"rise": _get_rise(*basepairs),
                    "slide": _get_slide(*basepairs),
                    "shift": _get_shift(*basepairs),
                    "twist": _get_twist(*basepairs),
                    "tilt": _get_tilt(*basepairs),
                    "roll": _get_roll(*basepairs)}
        return geometry, geometry

    def _get_base_plane(self, res):

        atom = []
        for atom_name in ["C2", "C4", "C6", "C1'"]:
            A = res.atoms.select_atoms("name " + atom_name)[0]
            atom.append(A.position)

        n0 = _norm(np.cross((atom[1]-atom[0]), (atom[2]-atom[1])))
        diazine_center = sum(atom[:-1]) / 3.
        if res.resname in ["ADE", "GUA"]:
            ref = res.atoms.select_atoms("name C8")[0].position
        else:
            ref = res.atoms.select_atoms("name C6")[0].position

        plane = {"n0": n0, "C1p": atom[-1], "diazine": diazine_center,
                 "C6C8": ref}

        return plane

    def eval_dh(self):
        """  dict of dihedrals for each residue (key= resid)
        """
        self.dh_quality = {}
        for res in self.u.residues:
            self.dh_quality[res.resindex] = self._get_dihedrals(res)
        return

    def _get_dihedrals(self, res):

        atoms, logic = self._get_residue_BB(res)
        dh_valid = self._get_valid_DH(*logic)
        if res.resname in ["ADE", "GUA"]:
            pyr = False
        else:
            pyr = True

        dh = self._get_dh_for_res(atoms, pyr, dh_valid)
        return dh

    def _get_residue_BB(self, res):
        iniSeg, terSeg, ter5 = False, False, False

        atoms = {}
        try:
            P = res.atoms.select_atoms("name P")[0]
            atoms["P"] = P.position

        except (KeyError, IndexError):
            ter5 = True

        for xx in BB_ATOMS[:-1]:
            atoms[xx] = res.atoms.select_atoms("name " + xx)[0].position
        if res.resname in ["ADE", "GUA"]:
            YY = PUR_ATOMS
        else:
            YY = PYR_ATOMS
        for yy in YY:
            atoms[yy] = res.atoms.select_atoms("name " + yy)[0].position

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

    def _get_valid_DH(self, ter5, terSeg, iniSeg):
        dh_valid = ["gamma", "delta", "xi"]
        if terSeg is False:
            dh_valid.extend(["epsilon", "zeta"])
        if iniSeg is False:
            dh_valid.append("alpha")
        if ter5 is False:
            dh_valid.append("beta")
        return dh_valid

    def _get_dh_for_res(self, atoms, pyr, dh_valid):
        dh = {}
        for dh_name in DH_ATOMS.keys():
            if dh_name in dh_valid:
                angle = self._get_angle(atoms, pyr, dh_name)
            else:
                angle = None
            dh[dh_name] = angle
        return dh

    def _get_angle(self, atoms, pyr, dh_name):

        def _dh_angle(p1, p2, p3, p4, as_rad=False):

            v1 = p2 - p1
            v2 = p3 - p2
            v3 = p4 - p3

            n1 = _norm(np.cross(v1, v2))
            n2 = _norm(np.cross(v2, v3))
            m1 = np.cross(n1, _norm(v2))

            x = np.dot(n1, n2)
            y = np.dot(m1, n2)

            angle = - np.arctan2(y, x)

            return angle if as_rad else np.rad2deg(angle)

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

        angle = _dh_angle(*p)

        return angle

    def eval_distances(self):
        """  dict of twist for each wc-pair (key= resid-scaffold)
        """
        self.distances = {}
        dist_C1p = self._get_A_distance("C1'")
        dist_P = self._get_A_distance("P")
        for resindex in range(self.u.residues.n_residues):
            self.distances[resindex] = {
                "C1'": dist_C1p[resindex], "P": dist_P[resindex]}

        return

    def _get_A_distance(self, atomname, n=2):
        d_c = {}
        for res in self.u.residues:
            resindex = res.resindex
            A = res.atoms.select_atoms("name " + atomname)
            dist_strand = []
            dist_compl = []
            try:
                resindex_wc = self.d_Fbp_full[resindex]
            except KeyError:
                resindex_wc = None

            for i in range(-n, n+1):
                resindex_x, resindex_x_wc = self._get_next_wc(resindex,
                                                              resindex_wc, i)

                if resindex_x is not None:
                    res_x = self.u.residues[resindex_x]
                    B = res_x.atoms.select_atoms("name " + atomname)
                    if len(A) + len(B) == 2:  # TODO use  try atomselect[0]
                        strand = np.linalg.norm(A.positions - B.positions)
                    else:
                        strand = None
                else:
                    strand = None
                dist_strand.append(strand)

                if resindex_x_wc is not None:
                    res_x_wc = self.u.residues[resindex_x_wc]
                    C = res_x_wc.atoms.select_atoms("name " + atomname)
                    if len(A) + len(C) == 2:
                        compl = np.linalg.norm(A.positions - C.positions)
                    else:
                        compl = None
                else:
                    compl = None
                dist_compl.append(compl)

            d_c[resindex] = {"strand": dist_strand, "compl": dist_compl}
        return d_c

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
                    wc_r = self.d_Fbp_full[r]
                    bases.append(self._get_base_plane(self.u.residues[r]))
                    bases.append(self._get_base_plane(self.u.residues[wc_r]))
                except KeyError:
                    bases.append(self._get_base_plane(self.u.residues[r]))
                    bases.append(self._get_base_plane(self.u.residues[r]))
            bp_planes = []
            for i in [0, 2, 4, 6]:
                bp_planes.append(
                    {"center-C1p": ((bases[i]["C1p"] +
                                     bases[i+1]["C1p"]) * 0.5),
                     "dir-C1p": (bases[i]["C1p"] - bases[i+1]["C1p"]),
                     "center-diazine": ((bases[i]["diazine"] +
                                         bases[i+1]["diazine"]) * 0.5),
                     "dir-diazine": (bases[i]["diazine"] -
                                     bases[i+1]["diazine"]),
                     "center-C6C8": ((bases[i]["C6C8"] +
                                      bases[i+1]["C6C8"]) * 0.5),
                     "dir-C6C8": (bases[i]["C6C8"] -
                                  bases[i+1]["C6C8"]),
                     "n0": ((bases[i]["n0"] + bases[i+1]["n0"]) * 0.5)})
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
        for res_index, co in self.d_Fco.items():
            if res_index not in co_done:
                leg_index = co["leg"]
                co_index = co["co"]
                coleg_index = self.d_Fco[co_index]["leg"]

                co_done.update([res_index, co_index])
                # res -> a, co -> c
                bpplanes = get_co_baseplanes(
                    res_index, leg_index, co_index, coleg_index)

                co_type = co["type"][0]
                if co_type == "double":
                    double_res_index = co["type"][1]
                    double_leg_index = self.d_Fco[double_res_index]["leg"]
                    double_co_index = self.d_Fco[double_res_index]["co"]
                    double_coleg_index = self.d_Fco[double_co_index]["leg"]
                    co_done.update([double_res_index, double_co_index])
                    # double -> b, d_co -> d
                    double_bpplanes = get_co_baseplanes(
                        double_res_index, double_leg_index, double_co_index,
                        double_coleg_index)
                    co_data = get_co_angles_full(
                        bpplanes, double_bpplanes)  # ac bd
                    crossover_ids = (res_index,  double_res_index,
                                     co_index, double_co_index)
                elif (co_type == "single" and
                        self.d_Fco[co_index]["type"][0] != "end"):
                    single_res_index = co["type"][1]
                    single_leg_index = co["type"][2]
                    single_co_index = self.d_Fco[co_index]["type"][1]
                    single_coleg_index = self.d_Fco[co_index]["type"][2]
                    co_done.update([single_res_index, single_co_index])
                    # single -> b, s_co -> d
                    single_bpplanes = get_co_baseplanes(
                        single_res_index, single_leg_index, single_co_index,
                        single_coleg_index)
                    co_data = get_co_angles_full(
                        bpplanes, single_bpplanes)  # ac bd
                    crossover_ids = (res_index,  single_res_index,
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


def get_description():
    return """initializes MDAnalysis universe and compoutes watson scric base pairs.
    they are returned as to dictionaries. this process is repeated for each
     Hbond-deviation criterion
    subsequently universe and dicts are stored into a pickle. each deviation
    criterion is stored in one pickle"""


def proc_input():
    parser = argparse.ArgumentParser(
        description=get_description(),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument("--folder",
                        help="input folder",
                        type=str,
                        default="./",
                        )
    parser.add_argument("--name",
                        help="name of design and files",
                        type=str,
                        required="True",
                        default=argparse.SUPPRESS,
                        )
    parser.add_argument("--frames",
                        help="number of frames samples",
                        type=int,
                        default=2,
                        )
    parser.add_argument("--dev",
                        help="deviation of H-bond",
                        type=float,
                        default=0.1,
                        )
    args = parser.parse_args()

    Project = namedtuple("Project", ["input",
                                     "output",
                                     "name",
                                     "frames",
                                     "dev"
                                     ]
                         )
    project = Project(input=Path(args.folder),
                      output=Path(args.folder) / "analysis",
                      name=args.name,
                      frames=args.frames,
                      dev=args.dev,
                      )
    with ignored(FileExistsError):
        os.mkdir(project.output)
    return project


def write_pdb(u, bDNA, PDBs):
    u.add_TopologyAttr(
        mda.core.topologyattrs.Tempfactors(np.zeros(len(u.atoms))))

    u.atoms.tempfactors = -1.
    for res in u.residues:
        try:
            res.atoms.tempfactors = bDNA.wc_quality[res.resindex]["C1'C1'"]
        except KeyError:
            pass
    PDBs["qual"].write(u.atoms)

    for cond in WC_PROPERTIES:
        u.atoms.tempfactors = -1.
        for res in u.residues:
            try:
                res.atoms.tempfactors = (
                    bDNA.wc_geometry[res.resindex][cond]["center-C6C8"])
            except KeyError:
                pass
        PDBs[cond].write(u.atoms)

    for dh in DH_ATOMS.keys():
        u.atoms.tempfactors = -1.
        for res in u.residues:
            try:
                res.atoms.tempfactors = bDNA.dh_quality[res.resindex][dh]
            except KeyError:
                pass
        PDBs[dh].write(u.atoms)

    u.atoms.tempfactors = -1.
    ing = 0.00
    for resindex, resindex_wc in bDNA.d_Fbp.items():
        u.residues[resindex].atoms.tempfactors = ing
        u.residues[resindex_wc].atoms.tempfactors = ing
        ing += 0.01
    PDBs["bp"].write(u.atoms)

    return


def main():
    project = proc_input()

    linker = bpLinker.Linker(project)
    linkage = linker.create_linkage()
    u = mda.Universe(linkage.universe)
    frames_step = int(len(u.trajectory) / project.frames)
    frames = list(range(len(u.trajectory)-1, 0, -frames_step))

    for name, link in linkage._asdict().items():
        pickle.dump(link, open(str(project.output) + "__" + name + ".p", "wb"))

    properties = []
    traj_out = project.output / "frames"
    try:
        os.mkdir(traj_out)
    except FileExistsError:
        pass

    # open PDB files
    PDBs = {}
    for pdb_name in [*WC_PROPERTIES, "bp", "qual"]:
        PDBs[pdb_name] = mda.Writer(
            str(project.output) + "__wc_" + pdb_name + ".pdb", multiframe=True)
    for pdb_name in DH_ATOMS.keys():
        PDBs[pdb_name] = mda.Writer(
            str(project.output) + "__dh_" + pdb_name + ".pdb", multiframe=True)

    # loop over selected frames
    for i, ts in enumerate([u.trajectory[i] for i in frames]):
        print(ts)
        bDNA = BDna(u,
                    linkage.Fbp,
                    linkage.DidFid,
                    linkage.DhpsDid,
                    linkage.Fco,
                    linkage.Dskips
                    )

        # perform analyis
        print("eval_wc", project.name)
        bDNA.eval_wc()
        print("eval_distances", project.name)
        bDNA.eval_distances()
        print("eval_dh", project.name)
        bDNA.eval_dh()
        print("eval_co_angles", project.name)
        bDNA.eval_co_angles()
        # ipdb.set_trace()

        print("write pickle", project.name)
        properties.append(bDNA)
        props_tuple = [
            (bDNA.wc_geometry, "wc_geometry"), (bDNA.wc_quality, "wc_quality"),
            (bDNA.dh_quality, "dh_quality"), (bDNA.distances, "distances"),
            (bDNA.co_angles, "co_angles")]
        for prop, prop_name in props_tuple:
            pickle.dump((ts, prop), open(
                traj_out + project.name + "__bDNA-" + prop_name + "-" + str(i) + ".p",
                "wb"))
        print("write pdbs", project.name)
        write_pdb(u, bDNA, PDBs)

    # close PDB files
    for _, PDB in PDBs.items():
        PDB.close()

    return


if __name__ == "__main__":
    main()
