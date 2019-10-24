#!/usr/bin/env python3#
import attr

from typing import Tuple, List, Optional


@attr.s
class Crossover(object):
    Ps: List[Optional[Tuple[int, Optional[int]]]] = attr.ib()
    Ls: List[Optional[Tuple[int, Optional[int]]]] = attr.ib()
    P_pos: List[Optional[Tuple[int, Optional[int]]]] = attr.ib()
    L_pos: List[Optional[Tuple[int, Optional[int]]]] = attr.ib()
    typ: str = attr.ib()
    is_scaf: bool = attr.ib()

    def transform2bp(self, bps: dict):
        self.Ps, self.Ls = list(), list()
        for hp, hp_ in zip(self.P_pos, self.L_pos):
            self.Ps.append(bps.get(hp, None))
            self.Ls.append(bps.get(hp_, None))
