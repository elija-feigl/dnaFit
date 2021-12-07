#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
from dnaFit import __version__


def get_dependencies(path: str) -> str:
    with open(path) as f:
        return "\n\t".join(f.read().strip().splitlines())


def main():
    with open("VERSION", "w") as f:
        f.write(__version__)

    with open("./setup_config/base.cfg") as f:
        base = f.read()

    dependencies = {
        "requires": get_dependencies("requirements.txt"),
        "dev_requires": get_dependencies("requirements_dev.txt"),
    }

    with open("setup.cfg", "w") as f:
        f.write(base.format(version=__version__, **dependencies))


if __name__ == "__main__":
    main()
