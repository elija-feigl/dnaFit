
import MDAnalysis as mda
import attr
from operator import attrgetter
from typing import List, Tuple
from pathlib import Path


""" DESCR:
    Fit Class manageing Mdanalyis structure for a given trajectory and
    configuration file of a namd simulation.

    COMMENTS:
    16.11.2020 not covering multi-scaffold structures
"""


@attr.s
class Fit(object):
    conf: Path = attr.ib()
    top: Path = attr.ib()

    def __attrs_post_init__(self):
        self.u: "mda.universe" = self._get_universe()
        self.scaffold, self.staples = self._split_strands()

    def _get_universe(self) -> "mda.universe":
        # TODO: use logger
        if self.top.exists() and self.conf.exists():
            u = mda.Universe(str(self.top), str(self.conf))
        else:
            raise FileNotFoundError
        return u

    def _split_strands(self) -> Tuple["mda.segment", List["mda.segment"]]:
        # TODO: -low- multiscaffold
        strands = self.u.segments
        scaffold = max(strands, key=attrgetter("residues.n_residues"))
        staples = [strand for strand in strands if len(
            strand.atoms) != len(scaffold.atoms)]
        return scaffold, staples
