#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
""" atomic coordinate file type"""
from dataclasses import dataclass
from typing import List
from typing import Tuple
from typing import Union

import numpy as np
import numpy.typing as npt

from .types import AtomName
from .types import ChainID
from .types import Number
from .types import ResName


@dataclass
class Atom:
    """pdb/cif atom class"""

    i_atom_coor: Union[npt.NDArray[np.float64], List[str], Tuple[float, float, float]]
    i_atom_number: Union[int, str]
    i_atom_name: str
    i_res_name: str
    i_chain_id: Union[int, str]
    i_res_number: Union[int, str]
    i_opacity: Union[float, str] = 0.0
    i_temperature: Union[float, str] = 0.0

    def __post_init__(self):
        self.atom_coor: npt.NDArray[np.float64] = self._convert_coor_input(self.i_atom_coor)
        self.atom_number: Number = Number(self.i_atom_number)
        self.atom_name: AtomName = AtomName(self.i_atom_name)
        self.res_name: ResName = ResName(self.i_res_name)
        self.chain_id: ChainID = ChainID(self.i_chain_id)
        self.res_number: Number = Number(self.i_res_number)
        self.opacity: float = float(self.i_opacity) if float(self.i_opacity) != 0.0 else 1.0
        self.temperature: float = float(self.i_temperature)
        self.element: str = self.atom_name.element_name()

    def _convert_coor_input(
        self,
        inpt: Union[npt.NDArray[np.float64], List[str], Tuple[float, float, float]],
    ) -> npt.NDArray[np.float64]:
        if isinstance(inpt, np.ndarray):
            return inpt
        elif isinstance(inpt, list):
            coor = [(float(c) if isinstance(c, str) else c) for c in inpt]
            return np.array(coor)
        else:
            raise NotImplementedError

    def as_cif(self) -> str:
        """return atom in mmCIF format"""
        return "".join(
            [
                "ATOM ",
                self.atom_number.as_str().ljust(12, " "),
                self.element.ljust(2, " "),
                self.atom_name.as_cif().ljust(7, " "),
                ". ",
                self.res_name.as_str().ljust(3, " "),
                self.chain_id.as_cif().ljust(3, " "),
                self.chain_id.as_str().ljust(4, " "),
                self.res_number.as_str().ljust(6, " "),
                "? ",
                f"{self.atom_coor[0]: .3f}".ljust(9, " "),
                f"{self.atom_coor[1]: .3f}".ljust(9, " "),
                f"{self.atom_coor[2]: .3f}".ljust(9, " "),
                f"{self.opacity: .2f}".ljust(6, " "),
                f"{self.temperature: .2f}".ljust(8, " "),
                "? ",
                self.res_number.as_str().ljust(6, " "),
                self.res_name.as_str().ljust(3, " "),
                self.chain_id.as_chimera().ljust(3, " "),
                self.atom_name.as_cif().ljust(7, " "),
                "1",
                "\n",
            ]
        )

    def as_pdb(self) -> str:
        """return atom in PDB format"""
        return "".join(
            [
                "ATOM  ",
                self.atom_number.as_pdb4namd(width=5),
                self.atom_name.as_pdb4namd(width=5),
                self.res_name.as_pdb4namd(width=4),
                self.chain_id.as_pdb4namd(width=2),
                self.res_number.as_pdb4namd(width=4),
                (4 * " "),
                f"{self.atom_coor[0]: .3f}".rjust(8, " "),
                f"{self.atom_coor[1]: .3f}".rjust(8, " "),
                f"{self.atom_coor[2]: .3f}".rjust(8, " "),
                f"{self.opacity: .2f}".rjust(6, " "),
                f"{self.temperature: .2f}".rjust(6, " "),
                (6 * " "),
                self.chain_id.as_segName4namd(width=4),
                self.element,
                (2 * " "),  # charge
                "\n",
            ]
        )
