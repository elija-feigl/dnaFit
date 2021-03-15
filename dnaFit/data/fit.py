
import MDAnalysis as mda
from MDAnalysis.core.groups import Segment
import attr
import logging
from operator import attrgetter
from typing import List, Tuple
from pathlib import Path


""" DESCR:
    Fit Class managing Mdanalyis structure for a given topology and
    configuration/trajectory file of a namd simulation.

    COMMENTS:
    16.11.2020 not covering multi-scaffold structures
"""


@attr.s
class Fit(object):
    conf: Path = attr.ib()
    top: Path = attr.ib()

    def __attrs_post_init__(self):
        self.u: mda.Universe = self._get_universe()
        self.scaffold, self.staples = self._split_strands()
        self.logger = logging.getLogger(__name__)

    def _get_universe(self) -> mda.Universe:
        if self.top.exists() and self.conf.exists():
            u = mda.Universe(str(self.top), str(self.conf))
        else:
            self.logger.fatal(
                f"Failed to initialize mda.Universe due to missing files: {self.top} {self.conf}")
            raise FileNotFoundError
        return u

    def _split_strands(self) -> Tuple[Segment, List[Segment]]:
        # TODO: -low- multiscaffold
        strands = self.u.segments
        scaffold = max(strands, key=attrgetter("residues.n_residues"))
        staples = [strand for strand in strands if len(
            strand.atoms) != len(scaffold.atoms)]
        return scaffold, staples
