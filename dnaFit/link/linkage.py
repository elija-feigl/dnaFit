from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Set, Tuple

import MDAnalysis as mda
import pandas as pd

from ..data.basepair import BasePair
from ..data.crossover import Crossover

""" Linkage class stores the translation from cadnano base-indexing to
    namd-indexing of bases.

    COMMENTS:
"""


@dataclass
class Linkage(object):
    Fbp: Dict[int, int] = field(default_factory=dict)
    DidFid: Dict[int, int] = field(default_factory=dict)
    DhpsDid: Dict[Tuple[int, int, bool], int] = field(default_factory=dict)
    Dcolor: Dict[int, int] = field(default_factory=dict)
    Fnicks: Dict[int, int] = field(default_factory=dict)
    FidSeq_local: Dict[int, str] = field(default_factory=dict)
    FidSeq_global: Dict[int, str] = field(default_factory=dict)
    FidHN: Dict[int, List[int]] = field(default_factory=dict)
    Fco: Dict[str, Crossover] = field(default_factory=dict)
    Dhp_skips: Set[Tuple[int, int]] = field(default_factory=set)
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
        FidDhps_data = {Fid: self.DidDhps.get(Did, (-1, -1, -1))
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
