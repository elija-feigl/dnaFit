[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Python-version:](https://img.shields.io/badge/python-v3.8-green)]() | [Usage](#usage) | [Dependencies](#dependencies) | [Installation](#installation) | [References](#references)

# dnaFit: cascading mrDNA-driven flexibel fitting

A Python3 package that contains a collection of scripts related to the fitting of pseudoatomic models into cryo-EM maps of DNA-Origami structures.


## Usage
 * `dnaFit`: create an atomic model for a cryo-EM perform mrDNA simulation followed by cascading mrDNA-driven flexibel fitting of a selected cadnano design
    design and model are linked to allow analysis of basepairs and comparison between map and design
 link design to existing atomic model to allow analysis of basepairs and comparison between map and design
...


## Dependencies

* Python >= 3.8
  * numpy >= 1.20
  * click >= 8.0
  * mdanalysis >= 1.0
  * mrcfile >= 1.3
  * mrdna >= 1.0a0.dev51
  * nanodesign >= 1.0 (python3 fork)

### without mrDNA:

* Operating system independet
* [NAMD](https://www.ks.uiuc.edu/Research/namd/) molecular dynamics

### with mrDNA:

* Linux operating system
* g++ >= 4.8
* [CUDA toolkit](https://developer.nvidia.com/cuda-toolkit) >= 6
* [ARBD](http://bionano.physics.illinois.edu/arbd) simulation engine
* [NAMD](https://www.ks.uiuc.edu/Research/namd/) molecular dynamics

## Installation

install all required Python Packages. the following are not available via the PyPI
### mrDNA install
link: https://gitlab.engr.illinois.edu/tbgl/tools/mrdna
NOTE: requires install from git repository (git clone ...) to verify version number

### nanodesign3
python3 fork of autodesk/nanodesign
https://github.com/elija-feigl/nanodesign


    git clone https://github.com/elija-feigl/DNA_Fit
    cd dnaFit
    python setup.py install


## References

When using dnaFit Python Package in published work, please cite the following paper:

Kube, M., Kohler, F., Feigl, E. et al. Revealing the structures of megadalton-scale DNA complexes with nucleotide resolution. Nat Commun 11, 6229 (2020). https://doi.org/10.1038/s41467-020-20020-7
