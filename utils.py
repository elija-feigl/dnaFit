#!/usr/bin/env python3

import numpy as np

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
