#!/usr/bin/env python3#
import attr

from typing import Tuple, List, Optional

import MDAnalysis as mda

from basepair import BasePair


@attr.s
class CrossoverPicklable(object):
    Ps: List[Optional[Tuple[Optional[int], Optional[int]]]] = attr.ib()
    Ls: List[Optional[Tuple[Optional[int], Optional[int]]]] = attr.ib()
    P_pos: List[Optional[Tuple[Optional[int], Optional[int]]]] = attr.ib()
    L_pos: List[Optional[Tuple[int, Optional[int]]]] = attr.ib()
    typ: str = attr.ib()
    is_scaf: bool = attr.ib()

    def transform(self, u: "mda.universe"):
        Ps: List[Optional[BasePair]] = list()
        Ls: List[Optional[BasePair]] = list()
        for Xs, hp, Xs_new, in [(self.Ps, self.P_pos, Ps),
                                (self.Ls, self.L_pos, Ls)]:

            for idx, P in enumerate(Xs):
                sc_index, st_index = P if P is not None else (None, None)
                sc = None if sc_index is None else u.residues[sc_index]
                st = None if st_index is None else u.residues[st_index]
                if sc is None and st is None:
                    Xs_new.append(None)
                else:
                    Xs_new.append(BasePair(sc=sc, st=st, hp=hp[idx]))
        return Crossover(Ps=Ps,
                         Ls=Ls,
                         typ=self.typ,
                         is_scaf=self.is_scaf,
                         )


@attr.s
class Crossover(object):
    Ps: List[Optional[BasePair]] = attr.ib()
    Ls: List[Optional[BasePair]] = attr.ib()
    typ: str = attr.ib()
    is_scaf: bool = attr.ib()

    def transform2picklable(self):
        Ps: List[Optional[Tuple[Optional[int], Optional[int]]]] = list()
        Ls: List[Optional[Tuple[Optional[int], Optional[int]]]] = list()
        P_pos: List[Optional[Tuple[Optional[int], Optional[int]]]] = list()
        L_pos: List[Optional[Tuple[Optional[int], Optional[int]]]] = list()
        for Xs, Xs_new, X_pos in [(self.Ps, Ps, P_pos), (self.Ls, Ls, L_pos)]:
            for P in Xs:
                if P is None:
                    Xs_new.append(None)
                    X_pos.append(None)
                else:
                    sc_index = P.sc.resindex if P.sc is not None else None
                    st_index = P.st.resindex if P.st is not None else None
                    Xs_new.append((sc_index, st_index))
                    X_pos.append(P.hp)

        return CrossoverPicklable(Ps=Ps,
                                  Ls=Ls,
                                  P_pos=P_pos,
                                  L_pos=L_pos,
                                  typ=self.typ,
                                  is_scaf=self.is_scaf,
                                  )
