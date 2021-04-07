import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import attr
from nanodesign.converters import Converter
from nanodesign.data.base import DnaBase
from nanodesign.data.dna_structure import DnaStructure

""" Design Class manageing cadnano design. The autodesk/nanodesign package is
    used to read the json file. Some modifications and additions are necessary
    to the nanodesign functionality.

    COMMENTS:
    20.11.2020 A modified version of nanodesign is used (available on github)
    03.12.2020 anticipates singlescaffold structure with circular scaffold
"""


@attr.s
class Design(object):
    json: Path = attr.ib()
    seq: Optional[Path] = attr.ib()
    generated_with_mrdna: bool = attr.ib(default=True)

    def __attrs_post_init__(self):
        self.design = self._get_design()
        self.scaffold = self._scaffold()
        self.staples = self._staple()
        self.helixorder = self._create_helix_order()
        self.stapleorder = self._create_staple_order()
        self.allbases = [b for s in self.design.strands for b in s.tour]
        self.Dhps_base: dict = self._init_hps_base()
        self.Dhp_skips: Set[Tuple[int, int]] = set()
        self.logger = logging.getLogger(__name__)
        logging.getLogger("nanodesign").setLevel(logging.ERROR)

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
        try:
            scaffold, = [s.tour for s in self.design.strands if s.is_scaffold]
        except ValueError:
            self.logger.exception("Design contains multiple scaffold strands")
        self._close_strand(strand=scaffold)
        return scaffold

    def _staple(self) -> List[List[DnaBase]]:
        return [s.tour for s in self.design.strands if not s.is_scaffold]

    def _get_design(self) -> DnaStructure:
        converter = Converter()
        converter.modify = True
        if self.json.exists() and self.seq is None:
            converter.read_cadnano_file(
                file_name=str(self.json),
                seq_file_name=None,
                seq_name=None,
            )
        elif self.json.exists() and self.seq.exists():
            converter.read_cadnano_file(
                file_name=str(self.json),
                seq_file_name=str(self.seq),
                seq_name=None,
            )
        else:
            self.logger.error(
                f"Failed to initialize nanodesign due to missing files: {self.json} {self.seq}")
            raise FileNotFoundError
        converter.dna_structure.compute_aux_data()
        dnaStructure = converter.dna_structure

        # TODO: implement without reusing converter
        converter = Converter()
        converter.modify = False
        converter.read_cadnano_file(self.json, None, self.seq)
        converter.dna_structure.compute_aux_data()
        dna_structure_del_ins = converter.dna_structure

        hps_deletions = set()
        for strand in dna_structure_del_ins.strands:
            for base in strand.tour:
                if base.num_deletions != 0:
                    hps_deletions.add((base.h, base.p))

        self.Dhp_skips = hps_deletions
        return dnaStructure

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
            enrgMD: first occurrence of staple sorted by h, p
            nanodesign: 5' end of staple sorted by h, p
        -------
            Returns
            -------
            stapleorder
                nanodesign -> enrgMD
        """
        if self.generated_with_mrdna:
            Dhps = [(s[0].h, s[0].p)
                    for s in self.staples
                    ]
        else:
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
