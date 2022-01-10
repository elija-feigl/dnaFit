#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
import logging
from pathlib import Path
from typing import Literal
from typing import NamedTuple
from typing import Optional


def get_resource(resources: str) -> Path:
    return Path(__file__).parent / "resources" / resources


def _init_logging():
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s | [%(name)s] %(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)


ReleaseType = Optional[Literal["alpha", "beta", "candidate", "final", "dev"]]


class VersionInfo(NamedTuple):
    """Version class of dnaFit."""

    major: int
    minor: int
    micro: int

    release_level: ReleaseType = None
    serial: int = 1

    def __repr__(self) -> str:
        return f"{self.major}.{self.minor}.{self.micro}" + (
            f".{self.release_level}{self.serial}" * (self.release_level is not None)
        )


version_info = VersionInfo(0, 8, 5, "dev")
__version__ = repr(version_info)

_init_logging()
