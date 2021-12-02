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

""" Fit Class managing MDAnalysis structure for a given topology and
    configuration/trajectory file of a namd simulation.

    COMMENTS:
    16.11.2020 not covering multi-scaffold structures
"""

import logging
from dataclasses import dataclass
from operator import attrgetter
from pathlib import Path
from typing import List, Tuple

import MDAnalysis as mda
from MDAnalysis.core.groups import Segment


@dataclass
class Fit:
    """ atomic model class"""
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
                "Failed to initialize mda.Universe due to missing files: %s %s",
                self.top, self.conf)
            raise FileNotFoundError
        return universe

    def _split_strands(self) -> Tuple[Segment, List[Segment]]:
        # TODO: -low- multiscaffold
        strands = self.u.segments
        scaffold = max(strands, key=attrgetter("residues.n_residues"))
        staples = [strand for strand in strands if len(
            strand.atoms) != len(scaffold.atoms)]
        return scaffold, staples
