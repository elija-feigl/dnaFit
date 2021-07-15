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

from typing import List

__version__ = "0.6"
__authors__ = ["Elija Feigl"]
__copyright__ = "Copyright (C) 2021  Elija Feigl"
__credits__ = ["Autodesk: Nanodesign", "MDAnalysis", "mrcfile"]
__license__ = "GPL-3.0"
__maintainer__ = "Elija Feigl"
__email__ = "elija.feigl@tum.de"
__status__ = "Development"


def get_version() -> str:
    return __version__


def get_authors() -> List[str]:
    return __authors__
