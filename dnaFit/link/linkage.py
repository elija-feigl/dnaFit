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

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple

import MDAnalysis as mda
import pandas as pd

from ..data.basepair import BasePair
from ..data.crossover import Crossover


@dataclass
class Linkage(object):
    """ Linkage class stores the translation from cadnano base-indexing to
        namd-indexing of bases.

        COMMENTS:
        Insertions will be assigned a negative position value to retain unique
        keys. Multi-insertions (ins>1) receive an incremented positions.
        If multiple multi-insertions are close to each other, their positions
        might be ordered incorrectly.
    """
    Fbp: Dict[int, int] = field(default_factory=dict)
    DidFid: Dict[int, int] = field(default_factory=dict)
    DhpsDid: Dict[Tuple[int, int, bool], int] = field(default_factory=dict)
    Dcolor: Dict[int, int] = field(default_factory=dict)
    Fnicks: Dict[int, int] = field(default_factory=dict)
    FidSeq_local: Dict[int, str] = field(default_factory=dict)
    FidSeq_global: Dict[int, str] = field(default_factory=dict)
    FidHN: Dict[int, List[int]] = field(default_factory=dict)
    Fco: Dict[str, Crossover] = field(default_factory=dict)
    u: mda.Universe = None

    def __post_init__(self) -> None:
        self._reverse()

    def write_linkage(self, prefix: str, dest: Path) -> None:
        """ write human readable linkage information to dest folder
        """
        out = list()

        Fbp_file = dest / f"{prefix}_F-resID--basepairs.csv"
        Fbp_header = ["Atomic model resID_scaffold",
                      "Atomic model resID_staple"]
        Fbp_data = self.Fbp
        out.append((Fbp_file, Fbp_data, Fbp_header))

        FidDhps_file = dest / \
            f"{prefix}_F-resID--D-helix-position-strand.csv"
        FidDhps_header = ["Atomic model resID",
                          "cadnano helix", "cadnano position", "is_scaffold"]
        FidDhps_data = {Fid: self.DidDhps[Did]
                        for Fid, Did in self.FidDid.items()}
        out.append((FidDhps_file, FidDhps_data, FidDhps_header))
        # TODO: add more if needed

        try:
            for path, data, header in out:
                df = pd.DataFrame.from_dict(data=data, orient='index')
                df.reset_index(inplace=True)
                df.columns = header
                df.to_csv(path, index=False)
        except IOError:
            raise Exception("ERROR: write_linkage I/O error")

    def _reverse(self) -> None:
        def reverse_d(dict: dict) -> dict:
            return {v: k for k, v in iter(dict.items())}

        self.FidDid = reverse_d(self.DidFid)
        self.DidDhps = reverse_d(self.DhpsDid)
        self.Fbp_rev = reverse_d(self.Fbp)
        self.Fbp_full = {**self.Fbp, **self.Fbp_rev}

    # NOTE: why was this method needed by bdna?
    def relink_crossover_basepairs(self, bps: Dict[Tuple[int, int], BasePair]
                                   ) -> None:
        for co in self.Fco.values():
            Ps_relinked = list()
            for P in co.Ps:
                Ps_relinked.append(bps[P.hp] if P is not None else None)
            co.Ps = Ps_relinked

            Ls_relinked = list()
            for L in co.Ls:
                Ls_relinked.append(bps[L.hp] if L is not None else None)
            co.Ls = Ls_relinked
