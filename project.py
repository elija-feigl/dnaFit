#!/usr/bin/env python3

import attr
from pathlib import Path


@attr.s(slots=True)
class Project(object):
    input: Path = attr.ib()
    output: Path = attr.ib()
    name: str = attr.ib()
    # specific
    frames: int = attr.ib(default=1)
    dev: float = attr.ib(default=0.)
    ENmodify: bool = attr.ib(default=False)
    EN: str = attr.ib(default="11111110")
    context: int = attr.ib(default=5)
    range: int = attr.ib(default=10)
    halfmap: bool = attr.ib(default=False)
    localres: bool = attr.ib(default=False)
    star: bool = attr.ib(default=False)
