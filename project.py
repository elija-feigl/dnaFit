#!/usr/bin/env python3

import attr
from pathlib import Path


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
