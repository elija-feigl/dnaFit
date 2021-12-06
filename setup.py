#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.

""" Atomic model construction and analysis for lattice based DNAOrigami."""

from setuptools import find_packages, setup


setup(
    packages=find_packages(),
    include_package_data=True,
    entry_points='''
        [console_scripts]
        dnaFit=dnaFit.scripts.dnafit:cli
    ''',
)
