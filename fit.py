#!/usr/bin/env python
# -*- coding: utf-8 -*-3#

import MDAnalysis as mda  # type:ignore
import attr

from operator import attrgetter
from typing import List, Tuple
from pathlib import Path

from project import Project

_author__ = "Elija Feigl"
__copyright__ = "Copyright 2019, Dietzlab (TUM)"
__credits__ = ["Autodesk: Nanodesign", "MDAnalysis", "mrcfile"]
__license__ = "None"
__version__ = "0.4"
__maintainer__ = "Elija Feigl"
__email__ = "elija.feigl@tum.de"
__status__ = "Development"


@attr.s
class Fit(object):
    project: Project = attr.ib()

    def __attrs_post_init__(self):
        self.infile: Path = self.project.input / self.project.name
        self.u: "mda.universe" = self._get_universe()
        self.scaffold, self.staples = self._split_strands()

    def _get_universe(self) -> "mda.universe":
        top = self.infile.with_suffix(".psf")
        trj = self.infile.with_suffix(".dcd")
        # TODO: -mid- if pdb, try invoke vmd animate write dcd
        if top.exists() and trj.exists():
            u = mda.Universe(str(top), str(trj))
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
