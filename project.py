#!/usr/bin/env python
# -*- coding: utf-8 -*-3

import attr
from pathlib import Path

_author__ = "Elija Feigl"
__copyright__ = "Copyright 2019, Dietzlab (TUM)"
__credits__ = ["Autodesk: Nanodesign", "MDAnalysis", "mrcfile"]
__license__ = "None"
__version__ = "0.4"
__maintainer__ = "Elija Feigl"
__email__ = "elija.feigl@tum.de"
__status__ = "Development"


@attr.s(slots=True)
class Project(object):
    input: Path = attr.ib()
    output: Path = attr.ib()
    name: str = attr.ib()
    # specific FitAnalysis
    frames: int = attr.ib(default=1)
    dev: float = attr.ib(default=0.)
    relink: bool = attr.ib(default=False)
    localres: bool = attr.ib(default=False)
    pdb: bool = attr.ib(default=False)
    # specific FitLinker
    ENmodify: bool = attr.ib(default=False)
    EN: str = attr.ib(default="11111110")
    # specific segmentation
    context: int = attr.ib(default=5)
    range: int = attr.ib(default=10)
    halfmap: bool = attr.ib(default=False)
    star: bool = attr.ib(default=False)
