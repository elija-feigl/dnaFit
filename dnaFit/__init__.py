#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
import logging
from pathlib import Path


def get_version() -> str:
    return __version__


def get_resource(resources: str) -> Path:
    return Path(__file__).parent / "resources" / resources


def _init_logging():
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s | [%(name)s] %(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)


_init_logging()

version_info = [0, 10, 1, "alpha"]

__version__ = ".".join([str(sub) for sub in version_info])
__all__ = ["__version__"]
