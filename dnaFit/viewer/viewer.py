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

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List

import ipywidgets as widgets
import MDAnalysis as mda
import mrcfile
from IPython.display import display
from MDAnalysis.core.groups import AtomGroup

from ..data.mrc import write_mrc_from_atoms
from ..link.linkage import Linkage
from ..link.linker import Linker
from ..pdb.structure import Structure

"""
    VIEWERTOOL:
    This module defines the classes used to define the connectivity and
    geometry of a DNA structure.

    A DNA structure consists of a number of scaffold and staple strands
    (DNA origami), or oligo strands alone, bound together to form a designed
    geometric shape.

    NOTE: requires ipywidget and which is not listed in dnaFit requirements
"""


@dataclass
class Viewer(object):
    conf: Path
    top: Path
    mrc: Path
    json: Path
    seq: Path
    is_mrdna: bool = True

    def __post_init__(self):
        self.logger = self._setup_logger()
        self.linker = Linker(conf=self.conf, top=self.top, json=self.json, seq=self.seq,
                             generated_with_mrdna=self.is_mrdna)
        try:
            self.link: Linkage = self.linker.create_linkage()
        except:
            self.logger.error(
                "ERROR: The provided design is not compatible with the atomic.")

        self.u = self.linker.fit.u
        self.Hid2H = self.linker.design.design.structure_helices_map
        self.Hcoor2H = self.linker.design.design.structure_helices_coord_map

    def _setup_logger(self):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            '%(asctime)s | [%(name)s] %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        return logger

    def select_by_helixandbase(self, helices: List, bases: List):
        def _DhpsFid(h, p, s) -> int:
            return self.link.DidFid[
                self.link.DhpsDid[(h, p, s)]
            ]
        u = self.link.u
        stsc_bases = self._parse_selection(bases, helices)
        atoms = mda.AtomGroup([], u)
        for idx, bases in enumerate(stsc_bases):
            for base in bases:
                Fid = _DhpsFid(base.h, base.p, bool(idx))
                atoms += u.residues[Fid].atoms
        return atoms, self.link.Dcolor

    def select_dsDNA(self, atoms=None):
        if atoms is None:
            residues = self.u.residues
        else:
            residues = atoms.residues

        atoms_ds = self.empty_atomGroup()
        for residue in residues:
            if residue.resindex in self.link.Fbp_full.keys():
                atoms_ds += residue.atoms
        return atoms_ds

    def select_scaffold(self, atoms=None):
        if atoms is None:
            residues = self.u.residues
        else:
            residues = atoms.residues

        atoms_sc = self.empty_atomGroup()
        for residue in residues:
            if residue.resindex in self.link.Fbp.keys():
                atoms_sc += residue.atoms
        return atoms_sc

    def select_withoutH(self, atoms=None):
        if atoms is None:
            atoms = self.u.atoms

        atoms_noH = self.empty_atomGroup()
        for atom in atoms:
            if "H" not in atom.name:
                atoms_noH += atom
        return atoms_noH

    def write_custom_gridpdb(self, atoms_selected, name=None, destination=None):
        self.u.add_TopologyAttr('tempfactor')
        self.u.add_TopologyAttr('occupancy')

        for atom in self.link.u.atoms:
            if "H" not in atom.name and atom in atoms_selected:
                atom.tempfactor = atom.mass
                atom.occupancy = 1.
        if name is None:
            name = self.conf.name + "_grid"
        self.writepdb(atoms=self.u.atoms, name=name,
                      destination=destination, as_mmCif=False)
        self.u.atoms.write("/Users/elija/Desktop/grid_s.pdb", bonds=None)
        return

    def _parse_selection(self, base_pos: List, helix_ids: List) -> List:
        helices = [self.Hid2H[idx] for idx in helix_ids]
        scaffold = [
            b for h in helices for b in h.scaffold_bases if b.p in base_pos
        ]
        staples = [
            b for h in helices for b in h.staple_bases if b.p in base_pos
        ]
        return [staples, scaffold]

    def select_widget(self):
        def _button(r, c, lattice):
            h = self.Hcoor2H.get((r, c), None)
            h_id = "" if h is None else str(h.id)
            if lattice in ["square", "honeycomb"]:
                return widgets.ToggleButton(
                    description=h_id,
                    layout=widgets.Layout(width='30px', height='30px'),
                    disabled=(h is None),
                )
            else:
                raise TypeError

        def _minmax(alist: List, asrange=False):
            if asrange:
                return range(min(alist), max(alist) + 1)
            else:
                return min(alist), max(alist)

        minb, maxb = _minmax([base.p for base in self.linker.design.allbases])
        r_range = _minmax([coor[0] for coor in self.Hcoor2H.keys()], True)
        c_range = _minmax([coor[1] for coor in self.Hcoor2H.keys()], True)

        lattice = ("square" if self.linker.design.design.lattice_type == 0
                   else "honeycomb")
        helixButtons = [
            [_button(r, c, lattice) for c in c_range] for r in r_range
        ]

        buttonsBox = widgets.VBox([widgets.HBox(row) for row in helixButtons])
        baseSlider = widgets.IntRangeSlider(
            value=[minb, maxb], min=minb, max=maxb, step=1,
            description='base:',
            disabled=False, continuous_update=True, orientation='horizontal',
            readout=True, readout_format='d',
            layout=widgets.Layout(width='900px', height='70px')
        )
        contextSlider = widgets.IntSlider(
            value=4, min=1, max=10, step=1,
            description='context [Angstr.]:',
            disabled=False, continuous_update=True, orientation='vertical',
            readout=True, readout_format='d',
            layout=widgets.Layout(width='100px')
        )

        vBox = widgets.VBox([buttonsBox, baseSlider])
        mainWidget = widgets.Box(
            [vBox, contextSlider],
            layout=widgets.Layout(width='1000px', border='1px solid black')
        )
        display(mainWidget)
        return (helixButtons, baseSlider, contextSlider)

    def eval_sliders(self, helix_buttons, slider_b, slider_c):

        selection_bases = range(slider_b.lower, slider_b.upper + 1)

        flat_list = [item for sublist in helix_buttons for item in sublist]
        selection_helices = [int(b.description) for b in flat_list
                             if (not b.disabled and b.value)
                             ]

        with mrcfile.open(self.mrc, mode='r+') as mrc:
            voxel_size = mrc.voxel_size.x
        context = int(slider_c.value / voxel_size) + 1

        return (selection_helices, selection_bases), context

    def writepdb(self, atoms, name=None, destination=None, singleframe=True, frame=-1, as_mmCif=True):
        import warnings
        warnings.filterwarnings('ignore')

        if not len(atoms):
            self.logger.warning("Empty atom selection. No file written.")
            return

        if destination is None:
            self.logger.info(
                "No target destination provided. Using configuration file folder.")
            destination = self.conf.parent
        elif isinstance(destination, str):
            destination = Path(destination)
        else:
            self.logger.info(
                "Invalid target destination path. Using configuration file folder.")
            destination = self.conf.parent
        if name is None:
            self.logger.info(
                "No target name provided. Using configuration file name.")
            name = self.conf.name

        path = destination / f"{name}.pdb"

        with mda.Writer(path, multiframe=True,
                        n_atoms=atoms.n_atoms) as W:
            if singleframe:
                self.link.u.trajectory[frame]
                W.write(atoms)
            else:
                for _ in self.link.u.trajectory:
                    W.write(atoms)

        if as_mmCif:
            structure = Structure(path=path, remove_H=False)
            structure.parse_pdb()
            mmcif = self.conf.with_suffix(".cif")
            structure.write_cif(mmcif)

    def writedcd(self, atoms, name=None, destination=None):
        if not len(atoms):
            self.logger.warning("Empty atom selection. No file written.")
            return

        if destination is None:
            self.logger.info(
                "No target destination provided. Using configuration file folder.")
            destination = self.conf.parent
        elif isinstance(destination, str):
            destination = Path(destination)
        else:
            self.logger.info(
                "Invalid target destination path. Using configuration file folder.")
            destination = self.conf.parent
        if name is None:
            self.logger.info(
                "No target name provided. Using configuration file name.")
            name = self.conf.name

        path = destination / f"{name}.dcd"

        with mda.Writer(path.as_posix(), n_atoms=atoms.n_atoms) as W:
            for _ in self.link.u.trajectory:
                W.write(atoms)

    def writemrc(self, atomsXX, name=None, destination=None, context=4, cut_box=True):
        if destination is None:
            self.logger.info(
                "No target destination provided. Using configuration file folder.")
            destination = self.conf.parent
        elif isinstance(destination, str):
            destination = Path(destination)
        else:
            self.logger.info(
                "Invalid target destination path. Using configuration file folder.")
            destination = self.conf.parent
        if name is None:
            self.logger.info(
                "No target name provided. Using configuration file name.")
            name = self.conf.name

        path = destination / f"{name}.mrc"
        write_mrc_from_atoms(path=self.mrc, atoms=atomsXX,
                             path_out=path, context=context, cut_box=cut_box)

    def empty_atomGroup(self) -> AtomGroup:
        return AtomGroup([], self.u)
