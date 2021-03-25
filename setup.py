from dnaFit.version import get_version
from setuptools import setup, find_packages

"""Atomic model construction and analysis for lattice based DNAOrigami."""

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("LICENSE", "r") as fh:
    license = fh.read()

setup(
    name="dnaFit",
    version=get_version(),
    author="Elija Feigl",
    author_email="elija.feigl@tum.de",
    description=__doc__,
    license=license,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/elija-feigl/DNA_Fit",
    packages=find_packages(),
    include_package_data=True,
    install_requires=(
        'numpy',
        'attrs',
        'click',
        'mdanalysis>=1.0',
        'mrcfile>=1.0',
        # 'mrdna>=0.0'  # TODO: add link in readme
        # 'nanodesign>=???'  # TODO: add version
    ),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: GNU General Public License Version 3",
        "Operating System :: OS Independent",
    ),
    entry_points='''
        [console_scripts]
        dnaFit=dnaFit.scripts.dnafit:cli
    ''',
)
