from dnaFit.version import get_version
import setuptools

description = """\
Atomic model construction and analysis for lattice based DNAOrigami\
"""

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("LICENSE", "r") as fh:
    license = fh.read()


setuptools.setup(
    name="dnaFit",
    version=get_version(),
    author="Elija Feigl",
    author_email="elija.feigl@tum.de",
    description=description,
    license=license,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/elija-feigl/DNA_Fit",
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=(
        'numpy>=1.14',
        'attrs>=19.3',
        'mdanalysis>=1.0',
        'mrcfile>=1.1.1',
        'mrdna>=1.0'  # TODO: add link in readme
        'nanodesign>=???'  # TODO: add version
    ),
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: GNU General Public License Version 3",
        "Operating System :: OS Independent",
    ),
)
