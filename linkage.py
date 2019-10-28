#!/usr/bin/env python3#
import pickle
import attr

from typing import Set, Dict, Tuple
import MDAnalysis as mda

from project import Project
from crossover import Crossover, CrossoverPicklable


@attr.s(auto_attribs=True)
class Linkage(object):
    Fbp: Dict[int, int] = {}
    DidFid: Dict[int, int] = {}
    DhpsDid: Dict[Tuple[int, int, bool], int] = {}
    Dcolor: Dict[int, int] = {}
    Dskips: Set[Tuple[int, int]] = set()
    Fnicks: Dict[int, int] = {}
    FidSeq: Dict[int, str] = {}
    Fco: Dict[str, Crossover] = {}
    u: "mda.universe" = None

    def dump_linkage(self, project: Project) -> None:
        def pickle_universe(u: "mda.universe") -> Tuple[str, str]:
            top = u.filename
            suffix = top.split(".")[-1]
            trj = top[:len(suffix)] + "dcd"
            return (top, trj)

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

        def unpickle_Fco(Fco: Dict[str, Crossover]
                         ) -> Dict[str, CrossoverPicklable]:
            return {key: co.transform(u=self.u) for key, co in Fco.items()}

        for name in vars(self):
            input = project.output / "{}__{}.p".format(project.name, name)
            value = pickle.load(open(input, "rb"))
            if name == "Fco":
                value = unpickle_Fco(value)
            elif name == "u":
                value = unpickle_universe(value)
            setattr(self, name, value)

    def reverse(self) -> None:
        def reverse_d(dict: dict) -> dict:
            return {v: k for k, v in iter(dict.items())}

        self.FidDid = reverse_d(self.DidFid)
        self.DidDhps = reverse_d(self.DhpsDid)
        self.Fbp_rev = reverse_d(self.Fbp)
        self.Fbp_full = {**self.Fbp, **self.Fbp_rev}

    def relink_crossover_basepairs(self, bps):
        for co in self.Fco.values():
            co.Ps = [None if P is None else bps[P.hp] for P in co.Ps]
            co.Ls = [None if L is None else bps[L.hp] for L in co.Ls]
