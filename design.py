#!/usr/bin/env python
# -*- coding: utf-8 -*-3#
import attr

from typing import List, Dict, Set, Tuple
from nanodesign.converters import Converter
from nanodesign.data.base import DnaBase
from nanodesign.data.dna_structure import DnaStructure

from project import Project

""" DESCR:
    Design Class manageing cadnano design. The autodesk/nanodesign package is
    used to read the json file. Some modifiactions and additions are necessary
    to the nanodesign functionaly.

    COMMENTS:
    20.11.2020 A modified version of nanodesign is used (available on github)
    03.12.2020 anticipates circular scaffold structures
"""


@attr.s
class Design(object):
    project: Project = attr.ib()

    def __attrs_post_init__(self):
        self.infile = self.project.input / self.project.name
        self.design = self._get_design()
        self.scaffold = self._scaffold()
        self.staples = self._staple()
        self.helixorder = self._create_helix_order()
        self.stapleorder = self._create_staple_order()
        self.allbases = [b for s in self.design.strands for b in s.tour]
        self.Dhps_base: dict = self._init_hps_base()
        self.Dhp_skips: Set[Tuple[int, int]] = set()

    def _init_hps_base(self):
        hps_base = dict()
        for base in self.allbases:
            position = (base.h, base.p, base.is_scaf)
            hps_base[position] = base
        return hps_base

    def _close_strand(self, strand):
        start = strand[0]
        end = strand[-1]
        start.up, end.down = end, start
        self.design.strands[start.strand].is_circular = True
        return strand

    def _scaffold(self) -> List[DnaBase]:
        # TODO: -low- multiscaffold
        scaffold = [s.tour for s in self.design.strands if s.is_scaffold][0]
        self._close_strand(strand=scaffold)
        return scaffold

    def _staple(self) -> List[List[DnaBase]]:
        return [s.tour for s in self.design.strands if not s.is_scaffold]

    def _get_design(self) -> DnaStructure:
        fil = self.infile.with_suffix(".json")
        seq = self.infile.with_suffix(".seq")
        converter = Converter(modify=True)
        if fil.exists() and seq.exists():
            converter.read_cadnano_file(
                file_name=str(fil),
                seq_file_name=str(seq),
                seq_name=None,
            )
        else:
            raise FileNotFoundError
        converter.dna_structure.compute_aux_data()
        self.Dhp_skips = converter.dna_structure.Dhp_skips
        return converter.dna_structure

    def _create_helix_order(self) -> Dict[int, int]:
        """ helices are not listed in the json in same order as they are listed
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
        """ enrgMD and nanodesign number staples differently.
            enrgMD: first occurence of staple sorted by h, p
            nanodesign: 5' end of staple sorted by h, p
        -------
            Returns
            -------
            stapleorder
                nanodesign -> enrgMD
        """
        Dhps = [(self.helixorder[s[0].h], s[0].p)
                for s in self.staples
                ]
        Dhps_sorted = sorted(Dhps, key=lambda x: (x[0], x[1]))
        order_ND = [
            Dhps.index(Dhps_sorted[i])
            for i, _ in enumerate(Dhps)
        ]
        stapleorder = {nd: idx for (idx, nd) in enumerate(order_ND)}
        return stapleorder
