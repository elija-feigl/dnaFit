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

import os
import subprocess
from pathlib import Path
from typing import List

import numpy as np

C1P_BASEDIST: float = 10.7
TOL: float = 10e-6
WC_DICT: dict = {
    "DC": "DG", "DG": "DC", "DT": "DA", "DA": "DT",
    "C": "G", "G": "C", "T": "A", "A": "T",
    "CYT": "GUA", "GUA": "CYT", "THY": "ADE", "ADE": "THY"
}
WC_HBONDS: dict = {
    "DA": ("N6", "N1"), "A": ("N6", "N1"), "ADE": ("N6", "N1"),
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
    "xi": {
        "pyr": ("C4'", "C1'", "N1", "C2"),
        "pur": ("C4'", "C1'", "N9", "C4"),
    },
}
BB_ATOMS: list = ["C1'", "O3'", "C3'", "C4'", "O5'", "C5'", "P"]
PYR_ATOMS: list = ["N1", "C2"]
PUR_ATOMS: list = ["N9", "C4"]

WC_PROPERTIES: list = ["rise", "slide", "shift", "twist", "tilt", "roll"]


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


def _get_executable(name: str):
    # TODO: use pathlib only
    for path in os.environ["PATH"].split(os.pathsep):
        exe = os.path.join(path.strip('"'), name)
        if os.path.isfile(exe) and os.access(exe, os.X_OK):
            return exe
    raise Exception(f"{name} was not found")


def _exec(cmd, logfile: Path):
    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, universal_newlines=True)
    with logfile.open(mode='w') as log:
        for line in process.stdout:
            log.write(line)
            log.flush()


def _check_path(filepath: str, extensions: list):
    path = Path(filepath).resolve()
    if path.exists():
        if path.suffix in extensions:
            return path
        else:
            raise Exception(
                f"input file {filepath} does not have correct extension {extensions}")
    raise Exception(f"input file {filepath} was not found")
