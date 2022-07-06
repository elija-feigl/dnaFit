#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
""" utility module """
import shutil
import subprocess
from pathlib import Path
from typing import Dict
from typing import List
from typing import Tuple
from typing import Union

import numpy as np
import numpy.typing as npt

INS_OFFSET: float = 0.1  # offset base position in design to account for insertions
C1P_BASEDIST: float = 10.7
TOL: float = 10e-6
WC_DICT: Dict[str, str] = {
    "DC": "DG",
    "DG": "DC",
    "DT": "DA",
    "DA": "DT",
    "C": "G",
    "G": "C",
    "T": "A",
    "A": "T",
    "CYT": "GUA",
    "GUA": "CYT",
    "THY": "ADE",
    "ADE": "THY",
}
WC_HBONDS: Dict[str, Union[Tuple[str, str], Tuple[str, str, str]]] = {
    "DA": ("N6", "N1"),
    "A": ("N6", "N1"),
    "ADE": ("N6", "N1"),
    "DT": ("O4", "N3"),
    "T": ("O4", "N3"),
    "THY": ("O4", "N3"),
    "DG": ("O6", "N1", "N2"),
    "G": ("O6", "N1", "N2"),
    "GUA": ("O6", "N1", "N2"),
    "DC": ("N4", "N3", "O2"),
    "C": ("N4", "N3", "O2"),
    "CYT": ("N4", "N3", "O2"),
}
(N6O4, N1N3_AT) = (2.95, 2.88)
(O6N4, N1N3_GC, N2O2) = (2.85, 2.85, 2.85)
WC_HBONDS_DIST: Dict[str, Union[Tuple[float, float], Tuple[float, float, float]]] = {
    "DA": (N6O4, N1N3_AT),
    "A": (N6O4, N1N3_AT),
    "ADE": (N6O4, N1N3_AT),
    "DT": (N6O4, N1N3_AT),
    "T": (N6O4, N1N3_AT),
    "THY": (N6O4, N1N3_AT),
    "DG": (O6N4, N1N3_GC, N2O2),
    "G": (O6N4, N1N3_GC, N2O2),
    "GUA": (O6N4, N1N3_GC, N2O2),
    "DC": (O6N4, N1N3_GC, N2O2),
    "C": (O6N4, N1N3_GC, N2O2),
    "CYT": (O6N4, N1N3_GC, N2O2),
}
DH_ATOMS: Dict[str, Union[Tuple[str, str, str, str], Dict[str, Tuple[str, str, str, str]]]] = {
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
BB_ATOMS: List[str] = ["C1'", "O3'", "C3'", "C4'", "O5'", "C5'", "P"]
PYR_ATOMS: List[str] = ["N1", "C2"]
PUR_ATOMS: List[str] = ["N9", "C4"]

WC_PROPERTIES: List[str] = ["rise", "slide", "shift", "twist", "tilt", "roll"]


def _norm(vector: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """versor"""
    return vector / np.linalg.norm(vector)


def _proj(u: npt.NDArray[np.float64], v: npt.NDArray[np.float64]) -> float:
    """vector projection"""
    return np.inner(u, v) / np.linalg.norm(v)  # type: ignore


def _norm_proj(u: npt.NDArray[np.float64], v: npt.NDArray[np.float64]) -> float:
    """norm vector projection, cos(phi)"""
    return np.inner(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))  # type: ignore


def _v_proj(u: npt.NDArray[np.float64], v: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    """vector v projected on vector"""
    return _proj(u, v) * _norm(v)


def _project_v2plane(
    v: npt.NDArray[np.float64], n0: npt.NDArray[np.float64]
) -> npt.NDArray[np.float64]:
    v_projected = _v_proj(v, n0)
    return v - v_projected


def _proj2plane(
    vect: List[npt.NDArray[np.float64]], n0: npt.NDArray[np.float64]
) -> npt.NDArray[np.float64]:
    pro = []
    for x in vect:
        d = np.inner(x, n0)
        p = [d * n0[i] for i in range(len(n0))]
        pro.append([x[i] - p[i] for i in range(len(x))])
    return pro


def _save_arccos_deg(projection: float) -> float:
    if 1.0 < abs(projection) < 1.0 + TOL:
        projection = np.sign(projection)
    if projection > 0:
        phi = np.arccos(projection)
    else:
        phi = -np.arccos(abs(projection))
    return np.rad2deg(phi)  # type: ignore


# def _dh_angle(p: List[], as_rad=False):  # slow
#     """dihedral angle."""

#     v1 = p[1] - p[0]
#     v2 = p[2] - p[1]
#     v3 = p[3] - p[2]

#     n1 = _norm(np.cross(v1, v2))
#     n2 = _norm(np.cross(v2, v3))
#     m1 = np.cross(n1, _norm(v2))

#     x = np.dot(n1, n2)
#     y = np.dot(m1, n2)

#     angle = -np.arctan2(y, x)
#     return angle if as_rad else np.rad2deg(angle)


def _get_executable(name: str):
    exe = shutil.which(name)
    if exe is not None:
        return exe
    raise OSError(f"{name} was not found")


def _exec(cmd, logfile: Path):
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    with logfile.open(mode="w") as log:
        if process.stdout is not None:
            for line in process.stdout:
                log.write(line)
                log.flush()


def _check_path(filepath: Path, extensions: List[str]):
    if filepath.suffix not in extensions:
        raise Exception(f"input file {filepath} does not have correct extension {extensions}")
