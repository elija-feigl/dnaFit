# DNA_Fit: cascading mrDNA-driven flexibel fitting

A Python3 package that contains a collection of scripts related to the fitting of pseudoatomic models into cryo-EM maps of DNA-Origami structures.


# Usage
 * `dnaFit`: create an atomic model for a cryo-EM perform mrDNA simulation followed by cascading mrDNA-driven flexibel fitting of a selected cadnano design
    design and model are linked to allow analysis of basepairs and comparison between map and design
 * `dnaLink`: link design to existing atomic model to allow analysis of basepairs and comparison between map and design
...


# Dependencies

* Python >= 3.7
  * numpy >= 1.14
  * scipy >= 1.1
  * mdanalysis >= 1.0 
  * mrcfile >= ??
  * attrs >= ??
  * nanodesign3  >= ?? # not publised yet

## without  mrDNA:

* Operating system independet

## with mrDNA:

* Linux operating system
* g++ >= 4.8
* [CUDA toolkit](https://developer.nvidia.com/cuda-toolkit) >= 6
* [ARBD](http://bionano.physics.illinois.edu/arbd) simulation engine


# Installation

install all required Python Packages. the following are not available via the PyPI
## mrDNA install
link: https://gitlab.engr.illinois.edu/tbgl/tools/mrdna
## nanodesign3
https://github.com/elija-feigl/nanodesign_dietz


git clone https://github.com/elija-feigl/DNA_Fit
cd dnaFit
python setup.py install
cd ..


# Citation

When using dnaFit Python Package in published work, please cite the following paper:

Kube, M., Kohler, F., Feigl, E. et al. Revealing the structures of megadalton-scale DNA complexes with nucleotide resolution. Nat Commun 11, 6229 (2020). https://doi.org/10.1038/s41467-020-20020-7
