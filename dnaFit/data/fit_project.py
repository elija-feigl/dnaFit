
import MDAnalysis as mda
import attr
from operator import attrgetter
from typing import List, Tuple
from pathlib import Path

from ..core.project import Project

""" DESCR:
    Design Class manageing Mdanalyis structure for a given trajectory and
    configuration file of a namd simulation.

    COMMENTS:
    16.11.2020 not covering multi-scaffold structures
"""


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
        if not trj.exists():
            trj = self.infile.with_suffix(".coor")
            print("try to use  coor file instead")
        if not trj.exists():
            trj = self.infile.with_suffix(".pdb")
            print("try to use  pdb file instead")

        # TODO: use logger
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
