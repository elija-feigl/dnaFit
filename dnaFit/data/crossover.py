
from dataclasses import dataclass
from typing import List, Optional

from .basepair import BasePair

""" Crossover Classes represent holiday-junction in DNA-Origami
"""


@dataclass
class Crossover(object):
    Ps: List[Optional[BasePair]]
    Ls: List[Optional[BasePair]]
    typ: str  # [full, half, end]
    is_scaf: bool
