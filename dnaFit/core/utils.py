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


def _norm(vector):
    """return versor"""
    return vector / np.linalg.norm(vector)


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
