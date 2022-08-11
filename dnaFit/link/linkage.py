#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
""" Linkage class module
    NOTE: (01.02.2021) does no longer support legacy pickle linkage
"""
from dataclasses import dataclass
from dataclasses import field
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Tuple

import MDAnalysis as mda
import pandas as pd

from ..data.crossover import Crossover


@dataclass
class Linkage:
    """Linkage class stores the translation from cadnano base-indexing to
    namd-indexing of bases.

    COMMENTS:
    Insertions will be assigned a negative position value to retain unique
    keys. Multi-insertions (ins>1) receive an incremented positions.
    If multiple multi-insertions are close to each other, their positions
    might be ordered incorrectly.
    """

    u: mda.Universe
    Fbp: Dict[int, int] = field(default_factory=dict)
    DidFid: Dict[int, int] = field(default_factory=dict)
    DhpsDid: Dict[Tuple[int, int, bool], int] = field(default_factory=dict)
    Dcolor: Dict[int, int] = field(default_factory=dict)
    Fnicks: Dict[int, int] = field(default_factory=dict)
    FidSeq_local: Dict[int, str] = field(default_factory=dict)
    FidSeq_global: Dict[int, str] = field(default_factory=dict)
    FidHN: Dict[int, List[int]] = field(default_factory=dict)
    Fco: Dict[str, Crossover] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self._reverse()

    def write_linkage(self, prefix: str, dest: Path) -> None:
        """write human readable linkage information to dest folder"""
        out = list()

        Fbp_file = dest / f"{prefix}_F-resID--basepairs.csv"
        Fbp_header = ["Atomic model resID_scaffold", "Atomic model resID_staple"]
        Fbp_data = self.Fbp
        out.append((Fbp_file, Fbp_data, Fbp_header))

        FidDhps_file = dest / f"{prefix}_F-resID--D-helix-position-strand.csv"
        FidDhps_header = [
            "Atomic model resID",
            "cadnano helix",
            "cadnano position",
            "is_scaffold",
        ]
        FidDhps_data = {Fid: self.DidDhps[Did] for Fid, Did in self.FidDid.items()}
        out.append((FidDhps_file, FidDhps_data, FidDhps_header))

        try:
            for path, data, header in out:
                dataframe = pd.DataFrame.from_dict(data=data, orient="index")
                dataframe.reset_index(inplace=True)
                dataframe.columns = header
                dataframe.to_csv(path, index=False)
        except OSError as exc:
            raise Exception("ERROR: write_linkage I/O error") from exc

    def _reverse(self) -> None:
        def reverse_d(dictionary: Dict[Any, Any]) -> Dict[Any, Any]:
            return {v: k for k, v in iter(dictionary.items())}

        self.FidDid = reverse_d(self.DidFid)
        self.DidDhps = reverse_d(self.DhpsDid)
        self.Fbp_rev = reverse_d(self.Fbp)
        self.Fbp_full = {**self.Fbp, **self.Fbp_rev}
