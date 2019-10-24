#!/usr/bin/env python3

import attr
from pathlib import Path


@attr.s(slots=True)
class ProjectAnalysis(object):
    input: Path = attr.ib()
    output: Path = attr.ib()
    name: str = attr.ib()
    # specific
    frames: int = attr.ib()
    dev: float = attr.ib()


@attr.s(slots=True)
class ProjectLink(object):
    input: Path = attr.ib()
    output: Path = attr.ib()
    name: str = attr.ib()
    # specific
    frames: int = attr.ib()
    dev: float = attr.ib()
    ENmodify: bool = attr.ib()
    EN: str = attr.ib()
