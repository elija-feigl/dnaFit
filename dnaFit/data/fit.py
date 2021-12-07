#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
""" Fit Class managing MDAnalysis structure for a given topology and
    configuration/trajectory file of a namd simulation.

    COMMENTS:
    16.11.2020 not covering multi-scaffold structures
"""
import logging
from dataclasses import dataclass
from operator import attrgetter
from pathlib import Path
from typing import List
from typing import Tuple

import MDAnalysis as mda
from MDAnalysis.core.groups import Segment


@dataclass
class Fit:
    """atomic model class"""

    conf: Path
    top: Path

    def __post_init__(self):
        self.logger = logging.getLogger(__name__)
        self.u: mda.Universe = self._get_universe()
        self.scaffold, self.staples = self._split_strands()
        self.logger = logging.getLogger(__name__)

    def _get_universe(self) -> mda.Universe:
        if self.top.exists() and self.conf.exists():
            universe = mda.Universe(str(self.top), str(self.conf))
        else:
            self.logger.error(
                "Failed to initialize mda.Universe due to missing files: %s %s", self.top, self.conf
            )
            raise FileNotFoundError
        return universe

    def _split_strands(self) -> Tuple[Segment, List[Segment]]:
        # TODO: -low- multiscaffold
        strands = self.u.segments
        scaffold = max(strands, key=attrgetter("residues.n_residues"))
        staples = [strand for strand in strands if len(strand.atoms) != len(scaffold.atoms)]
        return scaffold, staples
