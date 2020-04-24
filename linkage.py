#!/usr/bin/env python
# -*- coding: utf-8 -*-3#
import pickle
import attr

from typing import Dict, Tuple, List, Set
import MDAnalysis as mda

from project import Project
from crossover import Crossover, CrossoverPicklable
from basepair import BasePair


""" DESCR:
    Linkage class stores the translation from cadnano base-indexing to
    namd-indexing of bases.

    COMMENTS:
"""


@attr.s(auto_attribs=True)
class Linkage(object):
    Fbp: Dict[int, int] = {}
    DidFid: Dict[int, int] = {}
    DhpsDid: Dict[Tuple[int, int, bool], int] = {}
    Dcolor: Dict[int, int] = {}
    Fnicks: Dict[int, int] = {}
    FidSeq_local: Dict[int, str] = dict()
    FidSeq_global: Dict[int, str] = dict()
    FidHN: Dict[int, List[int]] = {}
    Fco: Dict[str, Crossover] = {}
    Dhp_skips: Set[Tuple[int, int]] = set()
    u: "mda.universe" = None

    def __attrs_post_init__(self) -> None:
        self._reverse()

    def dump_linkage(self, project: Project) -> None:
        def pickle_universe(u: "mda.universe") -> Tuple[str, str]:
            top = project.input / u.filename
            suffix = u.filename.split(".")[-1]
            trj = project.input / "{}dcd".format(u.filename[:-len(suffix)])
            return (str(top.absolute()), str(trj.absolute()))

        def pickle_Fco(Fco: Dict[str, Crossover]
                       ) -> Dict[str, CrossoverPicklable]:
            return {key: co.transform2picklable() for key, co in Fco.items()}

        for name, link in vars(self).items():
            if name == "Fco":
                link = pickle_Fco(link)
            elif name == "u":
                link = pickle_universe(link)
            output = project.output / "{}__{}.p".format(project.name, name)
            pickle.dump(link, open(output, "wb"))

    def load_linkage(self, project: Project) -> None:
        def unpickle_universe(u: Tuple[str, str]) -> "mda.universe":
            return mda.Universe(*u)

        def unpickle_Fco(Fco: Dict[str, Crossover],
                         u: "mda.universe"
                         ) -> Dict[str, CrossoverPicklable]:
            return {key: co.transform(u=self.u) for key, co in Fco.items()}

        names = list(vars(self).keys())
        names.remove("Fco")
        names.append("Fco")
        for name in names:
            input = project.output / "{}__{}.p".format(project.name, name)
            value = pickle.load(open(input, "rb"))
            if name == "u":
                value = unpickle_universe(value)
            elif name == "Fco":
                value = unpickle_Fco(value, self.u)
            setattr(self, name, value)

    def _reverse(self) -> None:
        def reverse_d(dict: dict) -> dict:
            return {v: k for k, v in iter(dict.items())}

        self.FidDid = reverse_d(self.DidFid)
        self.DidDhps = reverse_d(self.DhpsDid)
        self.Fbp_rev = reverse_d(self.Fbp)
        self.Fbp_full = {**self.Fbp, **self.Fbp_rev}

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
