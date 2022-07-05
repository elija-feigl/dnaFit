#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
""" Design Class manageing cadnano design. The autodesk/nanodesign package is
    used to read the json file. Some modifications and additions are necessary
    to the nanodesign functionality.

    COMMENTS:
    20.11.2020 A modified version of nanodesign is used (available on github)
    03.12.2020 anticipates singlescaffold structure with circular scaffold
"""
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict
from typing import List
from typing import Optional

from nanodesign.converters.converter import Converter
from nanodesign.data.base import DnaBase
from nanodesign.data.dna_structure import DnaStructure


@dataclass
class Design:
    """DNA Origami design class, based on nanodesign"""

    json: Path
    seq: Optional[Path] = None
    reorder_helices: bool = True
    ss_reduction: bool = False  # NOTE: SNUPI2.0

    def __post_init__(self):
        self.logger = logging.getLogger(__name__)
        logging.getLogger("nanodesign").setLevel(logging.ERROR)
        self.design = self._get_design()
        self.strands = self.design.strands
        if self.ss_reduction:
            self._remove_ss_extension()
        self.scaffold = self._scaffold()
        self.staples = self._staple()

        self.helixorder = self._create_helix_order()
        self.stapleorder = self._create_staple_order()
        self.allbases = [b for s in self.design.strands for b in s.tour]
        self.Dhps_base = self._init_hps_base()

    def _init_hps_base(self):
        hps_base = dict()
        for base in self.allbases:
            position = (base.h, base.p, base.is_scaf)
            hps_base[position] = base
        return hps_base

    def _close_strand(self, tour: List[DnaBase]):
        start = tour[0]
        end = tour[-1]
        start.up, end.down = end, start
        self.design.strands[start.strand].is_circular = True
        return tour

    def _scaffold(self) -> List[DnaBase]:
        # TODO: -low- multiscaffold
        scaffolds = [s.tour for s in self.design.strands if s.is_scaffold]
        if len(scaffolds) > 1:
            self.logger.error("Design contains multiple scaffold strands")
            sys.exit(1)
        scaffold: List[DnaBase] = scaffolds.pop()
        self._close_strand(tour=scaffold)
        return scaffold

    def _staple(self) -> List[List[DnaBase]]:
        return [s.tour for s in self.design.strands if not s.is_scaffold]

    def _remove_ss_extension(self, terminal_threshold=10) -> None:
        def delete_base(base, strand):
            strand.tour.remove(base)
            helix = self.design.structure_helices_map[base.h]
            if base.is_scaf:
                helix.scaffold_bases.remove(base)
            else:
                helix.staple_bases.remove(base)
            if base.up is not None:
                base.up.down = None
            if base.down is not None:
                base.down.up = None
            if base in self.design.base_connectivity:
                self.design.base_connectivity.remove(base)
            del base
            strand.base_id_list.clear()

        for strand in self.strands:
            removable_bases = []
            for idx, base in enumerate(strand.tour):
                is_ss_staple = base.across is None and not base.is_scaf
                is_terminal = not (terminal_threshold < idx < len(strand.tour) - terminal_threshold)
                if is_ss_staple and is_terminal:
                    removable_bases.append(base)

            for base in removable_bases:
                delete_base(base, strand)

    def _get_design(self) -> DnaStructure:
        seq = str(self.seq) if self.seq is not None else None
        converter = Converter()
        converter.modify = True
        if self.json.exists():
            converter.read_cadnano_file(
                file_name=str(self.json),
                seq_file_name=seq,
                seq_name=None,
            )
        else:
            self.logger.error(
                "Failed to initialize nanodesign due to missing files: %s %s",
                self.json,
                self.seq,
            )
            raise FileNotFoundError
        converter.dna_structure.compute_aux_data()
        return converter.dna_structure

    def _create_helix_order(self) -> Dict[int, int]:
        """helices are not listed in the json in same order as they are listed
            by helix-index. NOTE: enrgMD is helix reorder sensitive.
        -------
            Returns
            -------
            helixorder
        """
        helices_dict = self.design.structure_helices_map
        helixorder = {i: h.load_order for (i, h) in iter(helices_dict.items())}
        return helixorder

    def _create_staple_order(self) -> Dict[int, int]:
        """enrgMD and nanodesign number staples differently.
            enrgMD: first occurrence of staple sorted by h, p
            nanodesign: 5' end of staple sorted by h, p
        -------
            Returns
            -------
            stapleorder
                nanodesign -> enrgMD
        """
        if self.reorder_helices:
            Dhps = [(s[0].h, s[0].p) for s in self.staples]
        else:
            Dhps = [(self.helixorder[s[0].h], s[0].p) for s in self.staples]
        Dhps_sorted = sorted(Dhps, key=lambda x: (x[0], x[1]))
        order_nanodesign = [Dhps.index(Dhps_sorted[i]) for i, _ in enumerate(Dhps)]
        stapleorder = {nd: idx for (idx, nd) in enumerate(order_nanodesign)}
        return stapleorder
