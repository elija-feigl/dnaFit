#!/usr/bin/env python3

import numpy as np
from typing import List

C1P_BASEDIST: float = 10.7
TOL: float = 10e-6
WC_DICT: dict = {"DC": "DG", "DG": "DC", "DT": "DA", "DA": "DT",
                 "C": "G", "G": "C", "T": "A", "A": "T",
                 "CYT": "GUA", "GUA": "CYT", "THY": "ADE", "ADE": "THY"
                 }
WC_HBONDS: dict = {"DA": ("N6", "N1"), "A": ("N6", "N1"), "ADE": ("N6", "N1"),
                   "DT": ("O4", "N3"), "T": ("O4", "N3"), "THY": ("O4", "N3"),
                   "DG": ("O6", "N1", "N2"), "G": ("O6", "N1", "N2"),
                   "GUA": ("O6", "N1", "N2"), "DC": ("N4", "N3", "O2"),
                   "C": ("N4", "N3", "O2"), "CYT": ("N4", "N3", "O2")
                   }
(N6O4, N1N3_AT) = (2.95, 2.88)
(O6N4, N1N3_GC, N2O2) = (2.85, 2.85, 2.85)
WC_HBONDS_DIST: dict = {
    "DA": (N6O4, N1N3_AT), "A": (N6O4, N1N3_AT),
    "ADE": (N6O4, N1N3_AT), "DT": (N6O4, N1N3_AT),
    "T": (N6O4, N1N3_AT), "THY": (N6O4, N1N3_AT),
    "DG": (O6N4, N1N3_GC, N2O2), "G": (O6N4, N1N3_GC, N2O2),
    "GUA": (O6N4, N1N3_GC, N2O2), "DC": (O6N4, N1N3_GC, N2O2),
    "C": (O6N4, N1N3_GC, N2O2), "CYT": (O6N4, N1N3_GC, N2O2)
}
DH_ATOMS: dict = {
    "alpha": ("O3' -", "P", "O5'", "C5'"),
    "beta": ("P", "O5'", "C5'", "C4'"),
    "gamma": ("O5'", "C5'", "C4'", "C3'"),
    "delta": ("C5'", "C4'", "C3'", "O3'"),
    "epsilon": ("C4'", "C3'", "O3'", "P +"),
    "zeta": ("C3'", "O3'", "P +", "O5' +"),
    "xi": {"pyr": ("C4'", "C1'", "N1", "C2"),
           "pur": ("C4'", "C1'", "N9", "C4")
           }
}
BB_ATOMS: list = ["C1'", "O3'", "C3'", "C4'", "O5'", "C5'", "P"]
PYR_ATOMS: list = ["N1", "C2"]
PUR_ATOMS: list = ["N9", "C4"]

WC_PROPERTIES: list = ["rise", "slide", "shift", "twist", "tilt", "roll"]


class UnexpectedCaseError(Exception):
    """Raised when a case occurs that makes no sense in the programs context"""
    pass


def _norm(vector):
    return vector / np.linalg.norm(vector)


def _proj(u, v):
    return np.inner(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))


def _v_proj(u, v):
    return np.inner(u, v) / (np.linalg.norm(v) * np.linalg.norm(v)) * v


def _proj2plane(vect: List["np.ndarray"], n0: "np.ndarray"):
    pro = []
    for x in vect:
        d = np.inner(x, n0)
        p = [d * n0[i] for i in range(len(n0))]
        pro.append([x[i] - p[i] for i in range(len(x))])
    return pro


def _save_arccos_deg(dist):
    if 1.0 < abs(dist) < 1.0 + TOL:
        dist = np.sign(dist)
    if dist > 0:
        a = np.arccos(dist)
    else:
        a = - np.arccos(abs(dist))
    return np.rad2deg(a)


def _dh_angle(p: list, as_rad=False):  # slow

    v1 = p[1] - p[0]
    v2 = p[2] - p[1]
    v3 = p[3] - p[2]

    n1 = _norm(np.cross(v1, v2))
    n2 = _norm(np.cross(v2, v3))
    m1 = np.cross(n1, _norm(v2))

    x = np.dot(n1, n2)
    y = np.dot(m1, n2)

    angle = - np.arctan2(y, x)

    return angle if as_rad else np.rad2deg(angle)
