#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2021  Elija Feigl
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see https://www.gnu.org/licenses/gpl-3.0.html.

""" Atomic model construction and analysis for lattice based DNAOrigami."""

from setuptools import find_packages, setup

from dnaFit.version import get_version

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
        'numpy>=1.20',
        'click>=8.0.0',
        'mdanalysis>=1.0',
        'mrcfile>=1.3',
        #'mrdna>=1.0',
        'nanodesign>=1.0',
    ),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: GNU General Public License Version 3",
        "Operating System :: OS Independent",
    ],
    entry_points='''
        [console_scripts]
        dnaFit=dnaFit.scripts.dnafit:cli
    ''',
)
