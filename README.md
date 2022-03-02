[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![Python-version:](https://img.shields.io/badge/python-v3.8-green)]() | [Usage](#usage) | [Dependencies](#dependencies) | [Installation](#installation) | [References](#references)

# dnaFit: cascading mrDNA-driven flexible fitting

A Python3 package that contains a collection of scripts related to the fitting of pseudoatomic models into cryo-EM maps of DNA-Origami structures.


## Usage
```
Usage: dnaFit [OPTIONS] COMMAND [ARGS]...

  Cascaded mrDNA-driven MD flexible fitting script module:

  dnaFit commands:
      main pipeline:
      1. mrdna:   mrDNA structure prediction with custom settings for cadnano file.
                  (includes creation of prep folder for fitting)
      2. vmd_info:    print info for rigid body docking with VMD
      3. fit:     shrink wrap fitting if of rigid body docked mrDNA prediction
                  includes mask, pdb2cif, and link

Options:
  -v, --version  Show __version__ and exit.
  -h, --help     Show this message and exit.

Commands:
  center-on-map  Recenter atomic model center-of-mass on mrc cryo map...
  fit            Cascaded mrDNA-driven MD flexible fitting (shrink wrap...
  link           Links structural information of the CADNANO design file...
  mask           Mask mrc map to fitted atomic model.
  mrdna          MrDNA simulation of CADNANO design file with custom...
  pdb2cif        Generate atomic model in mmCIF format from namd PDB...
  prep           Prepare mrDNA results for fitting "dnaFit fit".
  vmd-info       Print VMD command for rotation of pdb around center of...
```

## Dependencies

* see setup.cfg

### without mrDNA:

* Operating system independent
* [NAMD](https://www.ks.uiuc.edu/Research/namd/) molecular dynamics
* [VMD](https://www.ks.uiuc.edu/Research/vmd/) visual molecular dynamics
### with mrDNA:

* Linux operating system
* g++ >= 4.8
* [mrDNA](https://gitlab.engr.illinois.edu/tbgl/tools/mrdna) Multi-Resolution DNA simulations
* [CUDA toolkit](https://developer.nvidia.com/cuda-toolkit) >= 6
* [ARBD](http://bionano.physics.illinois.edu/arbd) simulation engine
* [NAMD](https://www.ks.uiuc.edu/Research/namd/) molecular dynamics
* [VMD](https://www.ks.uiuc.edu/Research/vmd/) visual molecular dynamics

## Installation
  ```
    pip install git+https://github.com/elija-feigl/DNA_Fit#egg=dnaFit
  ```
or
  ```
    git clone https://github.com/elija-feigl/DNA_Fit
    cd dnaFit
    pip install .
  ```

The simulation frameworks mrDNA + ARBD, NAMD + VMD have to be installed separately.


## References

When using dnaFit Python Package in published work, please cite the following paper:

Kube, M., Kohler, F., Feigl, E. et al. Revealing the structures of megadalton-scale DNA complexes with nucleotide resolution. Nat Commun 11, 6229 (2020). https://doi.org/10.1038/s41467-020-20020-7
