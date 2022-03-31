#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
"""
    VIEWERTOOL:
    This module defines the classes used to define the connectivity and
    geometry of a DNA structure.

    A DNA structure consists of a number of scaffold and staple strands
    (DNA origami), or oligo strands alone, bound together to form a designed
    geometric shape.

    NOTE: requires ipywidgets and which is not listed in dnaFit requirements but with the ViewerApp
"""
import logging
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from typing import List

import ipywidgets as widgets
import MDAnalysis as mda
import mrcfile
from IPython.display import display
from MDAnalysis.core.groups import AtomGroup
from pdb2cif.pdb.structure import Structure

from ..data.mrc import write_binary_mrc_from_atoms
from ..data.mrc import write_mrc_from_atoms
from ..link.linkage import Linkage
from ..link.linker import Linker


@dataclass
class Viewer:
    """FitViewer main app class"""

    conf: Path
    top: Path
    mrc: Path
    json: Path
    seq: Path
    is_mrdna: bool = True

    def __post_init__(self):
        self.logger = logging.getLogger(__name__)
        warnings.filterwarnings("ignore")
        self.linker = Linker(
            conf=self.conf,
            top=self.top,
            json=self.json,
            seq=self.seq,
            generated_with_mrdna=self.is_mrdna,
        )
        try:
            self.link: Linkage = self.linker.create_linkage()
        except Exception:
            self.logger.error("ERROR: The provided design is not compatible with the atomic.")

        self.u = self.linker.fit.u
        self.Hid2H = self.linker.design.design.structure_helices_map
        self.Hcoor2H = self.linker.design.design.structure_helices_coord_map

    def select_by_helixandbase(self, helix_range: List[int], base_range: List[int]):
        """select all atoms of the cadnano subset of helices and baseposition-range"""

        def DhpsFid(h, p, s) -> int:
            return self.link.DidFid[self.link.DhpsDid[(h, p, s)]]

        u = self.link.u
        stsc_bases = self._parse_selection(base_range, helix_range)
        atoms = mda.AtomGroup([], u)
        for idx, bases in enumerate(stsc_bases):
            for base in bases:
                Fid = DhpsFid(base.h, base.p, bool(idx))
                atoms += u.residues[Fid].atoms
        return atoms, self.link.Dcolor

    def select_ds_dna(self, atoms=None):
        """select all atoms that are part of double stranded DNA
        from provided atoms or the whole universe if None
        """
        if atoms is None:
            residues = self.u.residues
        else:
            residues = atoms.residues

        atoms_ds = self.empty_atom_group()
        for residue in residues:
            if residue.resindex in self.link.Fbp_full.keys():
                atoms_ds += residue.atoms
        return atoms_ds

    def select_scaffold(self, atoms=None):
        """select all scaffold atoms
        from provided atoms or the whole universe if None
        """
        if atoms is None:
            residues = self.u.residues
        else:
            residues = atoms.residues

        atoms_sc = self.empty_atom_group()
        for residue in residues:
            if residue.resindex in self.link.Fbp.keys():
                atoms_sc += residue.atoms
        return atoms_sc

    def select_without_H(self, atoms=None):
        """select all scaffold not containing H (hydrogen) in their name
        from provided atoms or the whole universe if None
        """
        if atoms is None:
            atoms = self.u.atoms

        atoms_no_hydrogen = self.empty_atom_group()
        for atom in atoms:
            if "H" not in atom.name:
                atoms_no_hydrogen += atom
        return atoms_no_hydrogen

    def _parse_selection(self, base_pos: List[int], helix_ids: List[int]) -> List[Any]:
        helices = [self.Hid2H[idx] for idx in helix_ids]
        scaffold = [b for h in helices for b in h.scaffold_bases if b.p in base_pos]
        staples = [b for h in helices for b in h.staple_bases if b.p in base_pos]
        return [staples, scaffold]

    def select_widget(self, lattice=None):
        """create selection widget
        helix selection: click helices to activate
        base_position: move sliders to specify base-position range
        context: distance for zoning surrounding selected atoms [in Angstrom]
        """

        def _button(row, clmn, lattice):
            h = self.Hcoor2H.get((row, clmn), None)
            h_id = "" if h is None else str(h.id)
            if lattice in ["square", "honeycomb"]:
                description = h_id
                button = widgets.ToggleButton(
                    description=description,
                    layout=widgets.Layout(width="30px", height="30px"),
                    disabled=(h is None),
                )
                if lattice == "square":
                    return button
                else:
                    spacer = widgets.ToggleButton(
                        description="",
                        layout=widgets.Layout(width="30px", height="15px"),
                        disabled=True,
                    )
                    odd_clmn = not clmn % 2
                    even_row = row % 2
                    if (odd_clmn + even_row) % 2:
                        return widgets.VBox([button, spacer])
                    else:
                        return widgets.VBox([spacer, button])
            else:
                raise TypeError

        def _minmax(alist: List[Any], asrange=False):
            if asrange:
                return range(min(alist), max(alist) + 1)
            else:
                return min(alist), max(alist)

        minb, maxb = _minmax([base.p for base in self.linker.design.allbases])
        row_range = _minmax([coor[0] for coor in self.Hcoor2H], True)
        clmn_range = _minmax([coor[1] for coor in self.Hcoor2H], True)

        if lattice is None or lattice not in ["square", "honeycomb"]:
            lattice = "honeycomb" if self.linker.design.design.lattice_type == 0 else "square"
            self.logger.debug("Using %s-lattice from design file.", lattice)
        self.logger.info("Using %s-lattice widget.", lattice)

        helix_buttons = [[_button(row, clmn, lattice) for clmn in clmn_range] for row in row_range]

        button_box = widgets.VBox([widgets.HBox(row) for row in helix_buttons])
        base_slider = widgets.IntRangeSlider(
            value=[minb, maxb],
            min=minb,
            max=maxb,
            step=1,
            description="base:",
            disabled=False,
            continuous_update=True,
            orientation="horizontal",
            readout=True,
            readout_format="d",
            layout=widgets.Layout(width="900px", height="70px"),
        )
        context_slider = widgets.IntSlider(
            value=4,
            min=1,
            max=10,
            step=1,
            description="context [Angstrom]:",
            disabled=False,
            continuous_update=True,
            orientation="vertical",
            readout=True,
            readout_format="d",
            layout=widgets.Layout(width="100px"),
        )

        v_box = widgets.VBox([button_box, base_slider])
        main_widget = widgets.Box(
            [v_box, context_slider],
            layout=widgets.Layout(width="1000px", border="1px solid black"),
        )
        display(main_widget)
        return (helix_buttons, base_slider, context_slider)

    def eval_sliders(self, helix_selection, slider_b, slider_c):
        """evaluate selection widget"""

        selection_bases = range(slider_b.lower, slider_b.upper + 1)

        selection_helices = []
        for helix in helix_selection:
            for button in helix:
                if isinstance(button, widgets.VBox):
                    try:
                        button = next(b for b in button.children if (not b.disabled))
                    except StopIteration:
                        continue
                if not button.disabled and button.value:
                    selection_helices.append(int(button.description))

        with mrcfile.open(self.mrc, mode="r+") as mrc:
            voxel_size = mrc.voxel_size.x
        context = int(slider_c.value / voxel_size) + 1

        return (selection_helices, selection_bases), context

    def _check_destination(self, destination):
        if destination is None:
            self.logger.info("No target destination provided. Using configuration file folder.")
            destination = self.conf.parent
        elif isinstance(destination, str):
            destination = Path(destination)
        else:
            self.logger.info("Invalid target destination path. Using configuration file folder.")
            destination = self.conf.parent
        return destination

    def _check_name(self, name):
        if name is None:
            self.logger.info("No target name provided. Using configuration file name.")
            name = self.conf.name
        return name

    def write_pdb(
        self,
        atoms,
        name=None,
        destination=None,
        single_frame=True,
        frame=-1,
        as_cif=True,
    ):
        """write an atomic coordinate file for an atom group"""

        name = self._check_name(name)
        destination = self._check_destination(destination)
        path = destination / f"{name}.pdb"

        with mda.Writer(path, multiframe=True, n_atoms=atoms.n_atoms) as mda_writer:
            if single_frame:
                _ = self.link.u.trajectory[frame]
                mda_writer.write(atoms)
            else:
                for _ in self.link.u.trajectory:
                    mda_writer.write(atoms)

        if as_cif:
            structure = Structure(path=path, remove_H=False)
            structure.parse_pdb()
            mmcif = path.with_suffix(".cif")
            structure.write_cif(mmcif)

    def write_custom_gridpdb(self, atoms_selected, name=None, destination=None):
        """create a custom grid.pdb file based on an group of atoms"""
        self.u.add_TopologyAttr("tempfactor")
        self.u.add_TopologyAttr("occupancy")

        for atom in self.link.u.atoms:
            if "H" not in atom.name and atom in atoms_selected:
                atom.tempfactor = atom.mass
                atom.occupancy = 1.0
        if name is None:
            name = self.conf.name + "_grid"
        self.write_pdb(atoms=self.u.atoms, name=name, destination=destination, as_cif=False)

    def write_dcd(self, atoms, name=None, destination=None):
        """write an atomic trajectory file for an atom group [only if dcd supplied to viewer]"""
        if not atoms:
            self.logger.warning("Empty atom selection. No file written.")
            return

        name = self._check_name(name)
        destination = self._check_destination(destination)
        path = destination / f"{name}.dcd"

        with mda.Writer(path.as_posix(), n_atoms=atoms.n_atoms) as mda_writer:
            for _ in self.link.u.trajectory:
                mda_writer.write(atoms)

    def write_mrc(self, atoms, name=None, destination=None, context=4, cut_box=True):
        """write an mrc file zoned for the area within context of an atom group"""
        name = self._check_name(name)
        destination = self._check_destination(destination)
        path = destination / f"{name}.mrc"
        write_mrc_from_atoms(
            path=self.mrc, atoms=atoms, path_out=path, context=context, cut_box=cut_box
        )

    def write_binary_mask(self, atoms, name=None, destination=None, context=4):
        """write a binary mrc mask file for the area within context of an atom group"""
        name = self._check_name(name)
        destination = self._check_destination(destination)
        path = destination / f"{name}.mrc"
        write_binary_mrc_from_atoms(path=self.mrc, atoms=atoms, path_out=path, context=context)

    def empty_atom_group(self) -> AtomGroup:
        """return an empty MDanalysis atom group"""
        return AtomGroup([], self.u)
