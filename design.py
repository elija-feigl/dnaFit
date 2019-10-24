#!/usr/bin/env python3#
import attr
import nanodesign as nd


from typing import List, Dict, Any
from nanodesign.converters import Converter

from project import ProjectLink as Project


@attr.s
class Design(object):
    project: Project = attr.ib()

    def __attrs_post_init__(self):
        self.infile = self.project.input / self.project.name
        self.design = self._get_design()
        self.strands = self.design.strands
        self.scaffold = self._clean_scaffold(self.strands)
        self.excl = self.scaffold[0].strand
        self.staples = self._clean_staple(self.strands)
        self.helixorder = self._create_helix_order()
        self.stapleorder = self._create_staple_order()
        self.allbases = [b for s in self.strands for b in s.tour]
        self.Dhps_base: dict = self._init_hps_base()

    def _init_hps_base(self):
        hps_base = dict()
        for base in self.allbases:
            position = (base.h, base.p, base.is_scaf)
            hps_base[position] = base
        return hps_base

    def _is_del(self, base: "nd.base") -> bool:
        return base.num_deletions != 0

    def _clean_scaffold(self, strands: List["nd.base"]) -> List["nd.base"]:
        # TODO: -low- multiscaffold
        # TODO: -low- insertions
        scaffold = [s.tour for s in strands if s.is_scaffold][0]
        scaffold_clean = [b for b in scaffold if not self._is_del(b)]
        return scaffold_clean

    def _clean_staple(self, strands: List["nd.base"]) -> List[List["nd.base"]]:
        # TODO: -low- insertions
        staples = [s.tour for s in strands if not s.is_scaffold]
        staples_clean = [[b for b in s if not self._is_del(b)]
                         for s in staples
                         ]
        return staples_clean

    def _get_design(self) -> Any:
        fil = self.infile.with_suffix(".json")
        seq = self.infile.with_suffix(".seq")
        converter = Converter()
        if fil.exists() and seq.exists():
            converter.read_cadnano_file(file_name=str(fil),
                                        seq_file_name=str(seq),
                                        seq_name=None,
                                        )
        else:
            raise FileNotFoundError
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
        order_ND = [Dhps.index(Dhps_sorted[i])
                    for i, _ in enumerate(Dhps)
                    ]
        stapleorder = {nd: idx for (idx, nd) in enumerate(order_ND)}
        return stapleorder
