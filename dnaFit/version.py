from typing import List

__version__ = "0.5"
__authors__ = ["Elija Feigl"]
__copyright__ = "Copyright 2021, Dietzlab (TUM)"
__credits__ = ["Autodesk: Nanodesign", "MDAnalysis", "mrcfile"]
__license__ = "GPL-3.0"
__maintainer__ = "Elija Feigl"
__email__ = "elija.feigl@tum.de"
__status__ = "Development"


def get_version() -> str:
    return __version__


def get_authors() -> List[str]:
    return __authors__
