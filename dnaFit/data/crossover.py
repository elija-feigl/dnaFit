#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
""" Crossovers are defined via a holiday-junction in DNA-Origami"""
from dataclasses import dataclass
from typing import List
from typing import Optional
from typing import Tuple


@dataclass
class Crossover:
    """crossover class
    to define angles, each crossover arm is composed of
        * positional base (part of the holiday junction)
        * leg base (n bases away)
    there are three types [full, half, end]
    """

    positional: List[Optional[Tuple[int, int]]]
    legs: List[Optional[Tuple[int, int]]]
    typ: str
    is_scaf: bool
