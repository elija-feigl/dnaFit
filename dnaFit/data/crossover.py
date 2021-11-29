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

""" Crossovers are defined via a holiday-junction in DNA-Origami"""

from dataclasses import dataclass
from typing import List, Optional

from .basepair import BasePair


@dataclass
class Crossover:
    """ crossover class
        to define angles, each crossover arm is composed of
            * positional base (part of the holiday junction)
            * leg base (n bases away)
        there are three types [full, half, end]
    """
    positionals: List[Optional[BasePair]]
    legs: List[Optional[BasePair]]
    typ: str
    is_scaf: bool
