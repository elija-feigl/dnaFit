
import attr

from typing import List, Optional

from .basepair import BasePair


""" Crossover Classes represent holiday-junction in DNA-Origami
"""


@attr.s
class Crossover(object):
    Ps: List[Optional[BasePair]] = attr.ib()
    Ls: List[Optional[BasePair]] = attr.ib()
    typ: str = attr.ib()  # [full, half, end]
    is_scaf: bool = attr.ib()
