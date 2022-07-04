#!/usr/bin/env python
# Copyright (C) 2021-Present  Elija Feigl
# Full GPL-3 License can be found in `LICENSE` at the project root.
""" module for linking cadnano design files with atomic model"""
import logging
from dataclasses import dataclass
from dataclasses import field
from pathlib import Path
from typing import Any
from typing import Dict
from typing import List
from typing import Optional
from typing import Tuple

import MDAnalysis as mda
from MDAnalysis.core.groups import AtomGroup
from nanodesign.data.base import DnaBase
from nanodesign.data.dna_structure_helix import DnaStructureHelix

from ..data.bdna import BDna
from ..data.crossover import Crossover
from ..data.design import Design
from ..data.fit import Fit
from .linkage import Linkage


@dataclass
class Linker:
    """Linker class

    COMMENTS:
    Insertions will be assigned a negative position value to retain unique
    keys. Multi-insertions (ins>1) receive an incremented positions.
    If multiple multi-insertions are close to each other, their positions
    might be ordered incorrectly.
    """

    conf: Path
    json: Path
    top: Path
    seq: Optional[Path] = None
    reorder_helices: bool = True
    link: Linkage = field(init=False)

    Fbp: Dict[int, int] = field(default_factory=dict)
    DidFid: Dict[int, int] = field(default_factory=dict)
    DhpsDid: Dict[Tuple[int, int, bool], int] = field(default_factory=dict)
    Fnicks: Dict[int, int] = field(default_factory=dict)
    FidSeq_local: Dict[int, str] = field(default_factory=dict)
    FidSeq_global: Dict[int, str] = field(default_factory=dict)
    Fco: Dict[str, Crossover] = field(default_factory=dict)

    FidSeq: Dict[int, str] = field(default_factory=dict)
    FidHN: Dict[int, List[int]] = field(default_factory=dict)
    Dcolor: Dict[int, int] = field(default_factory=dict)

    def __post_init__(self) -> None:

        self.fit: Fit = Fit(top=self.top, conf=self.conf)
        self.design: Design = Design(
            json=self.json, seq=self.seq, reorder_helices=self.reorder_helices
        )
        self.logger = logging.getLogger(__name__)

    def _eval_sequence(self, steps=5) -> None:
        """Affects
        -------
            self.FidSeq_local
            self.FidSeq_global

        """
        for base in self.design.allbases:
            # local: scaffold 5'->3'
            sequence = ""
            for stp in range(steps + 1):
                neighbor = self._get_n_strand(base=base, direct="down", steps=stp, local=True)
                if neighbor is None:
                    sequence += "N"
                elif neighbor.h != base.h:
                    sequence += "X"
                else:
                    n_resindex = self.DidFid[neighbor.id]
                    sequence += self.fit.u.residues[n_resindex].resname[0]
            resindex = self.DidFid[base.id]
            self.FidSeq_local[resindex] = sequence

            # global: even helix scaffold 5'->3', odd helix scaffold 3'->5'
            sequence = ""
            for stp in range(0, -(steps + 1), -1):
                neighbor = self._get_n_strand(
                    base=base,
                    direct="down",
                    steps=stp,
                    local=False,
                )
                if neighbor is None:
                    sequence += "N"
                elif neighbor.h != base.h:
                    sequence += "X"
                else:
                    n_resindex = self.DidFid[neighbor.id]
                    sequence += self.fit.u.residues[n_resindex].resname[0]
            resindex = self.DidFid[base.id]
            self.FidSeq_global[resindex] = sequence

    def _eval_FidHelixneighbors(self, steps=5) -> None:
        """Affects
        -------
            self.FidHN
        """

        def is_occupied_helix(
            H: Optional[DnaStructureHelix],
            p: int,
        ) -> bool:
            if H is None:
                return False
            elif (H.id, p, True) in self.design.Dhps_base:
                return True
            elif (H.id, p, False) in self.design.Dhps_base:
                return True
            else:
                return False

        for base in self.design.allbases:
            HidH = self.design.design.structure_helices_map
            HrcH = self.design.design.structure_helices_coord_map
            h_col = HidH[base.h].lattice_col
            h_row = HidH[base.h].lattice_row

            nhelices = list()
            for stp in range(1, steps + 1):
                number_nhelices = 0
                for rc in [
                    (h_row - stp, h_col),
                    (h_row + stp, h_col),
                    (h_row, h_col - stp),
                    (h_row, h_col + stp),
                ]:
                    n_H = HrcH.get(rc, None)
                    is_nh_occupied = is_occupied_helix(H=n_H, p=base.p)
                    if is_nh_occupied:
                        number_nhelices += 1
                nhelices.append(number_nhelices)
            resindex = self.DidFid[base.id]
            self.FidHN[resindex] = nhelices

    def create_linkage(self) -> Linkage:
        """invoke _link_scaffold, _link_staples, _link_bp to compute mapping
        of every base design-id to fit-id as well as the basepair mapping.
        basepairs are mapped from scaffold to staple, unique (invertable).
        updates linker attributes corresponding to the respective mapping
        and returns them.
        """
        try:
            self._sequence_match_link()
        except Exception:
            self.logger.warning(
                "Could not match design and model via strand sequences. Attempting strand-order specific matching."
            )
            self._link()

        self._identify_bp()
        self._identify_crossover()
        self._identify_nicks()
        self._eval_sequence()
        self._eval_FidHelixneighbors()
        self.link = Linkage(
            Fbp=self.Fbp,
            DidFid=self.DidFid,
            DhpsDid=self.DhpsDid,
            Dcolor=self.Dcolor,
            Fco=self.Fco,
            Fnicks=self.Fnicks,
            FidSeq_local=self.FidSeq_local,
            FidSeq_global=self.FidSeq_global,
            FidHN=self.FidHN,
            u=self.fit.u,
        )
        return self.link

    def _get_hp_list(self, strand: List[DnaBase], is_scaffold: bool) -> List[Tuple[int, int, bool]]:
        Dhp = list()
        for base in strand:
            hps = (base.h, base.p, is_scaffold)
            if hps in Dhp:  # insertion
                hps = (base.h, -base.p, is_scaffold)
            if hps in Dhp:  # multi-insertion
                pos = -base.p
                while hps in Dhp:
                    pos -= 1
                    hps = (base.h, pos, is_scaffold)
            Dhp.append(hps)
        return Dhp

    # NOTE: deprecated, only used as fallback fot _sequence_match_link() failure
    def _link(
        self,
    ) -> Tuple[Dict[int, int], Dict[Tuple[int, int, bool], int]]:
        def link_scaffold() -> Tuple[
            Dict[int, int],
            Dict[Tuple[int, int, bool], int],
        ]:
            """collect position in scaffold (0-x) by comparing index in list
                of scaffold_design positions
            -------
                Returns
                -------
                DidFid
                    design-id -> fit-id
                DhpsDid
                    helix-number, base-position, is_scaffold -> design-id
            """
            Dscaffold = self.design.scaffold
            Did = [base.id for base in Dscaffold]
            Dhp = self._get_hp_list(strand=Dscaffold, is_scaffold=True)
            Fid_local = [Did.index(base.id) for base in Dscaffold]
            Fid_global = self.fit.scaffold.residues[Fid_local].resindices

            DidFid = dict(zip(Did, Fid_global))
            DhpsDid = dict(zip(Dhp, Did))
            return (DidFid, DhpsDid)

        def link_staples() -> Tuple[
            Dict[int, int],
            Dict[Tuple[int, int, bool], int],
            Dict[int, int],
        ]:
            """same procedure as scaffold for each
            -------
            Returns
                -------
                DidFid
                    design-id -> fit-id
                DhpsDid
                    helix-number, base-position, is_scaffold -> design-id
                dict color
                    fit-segment-id -> color
            """

            def get_resid(segindex: int, resindex_local: int) -> Any:
                segment = self.fit.staples[segindex]
                return segment.residues[resindex_local].resindex

            DidFid: Dict[int, int] = {}
            DhpsDid: Dict[Tuple[int, int, bool], int] = {}
            color: Dict[int, int] = {}

            for i, staple in enumerate(self.design.staples):
                seg_id = self.design.stapleorder[i]

                Did = [base.id for base in staple]
                Dhp = self._get_hp_list(strand=staple, is_scaffold=False)

                Fid_local = [Did.index(base.id) for base in staple]
                Fid_global = [get_resid(seg_id, resid) for resid in Fid_local]

                icolor = self.design.design.strands[staple[0].strand].icolor
                segidxforcolor = self.fit.staples[seg_id].segindex
                color[segidxforcolor] = icolor

                DidFid_add = dict(zip(Did, Fid_global))
                DhpsDid_add = dict(zip(Dhp, Did))
                DidFid = {**DidFid, **DidFid_add}
                DhpsDid = {**DhpsDid, **DhpsDid_add}

            return (DidFid, DhpsDid, color)

        try:
            DidFid_sc, DhpsDid_sc = link_scaffold()
            DidFid_st, DhpsDid_st, self.Dcolor = link_staples()
        except IndexError:
            self.logger.exception("Linker failed: Design and Atomic Model incompatible")
            raise

        self.DidFid = {**DidFid_sc, **DidFid_st}
        self.DhpsDid = {**DhpsDid_sc, **DhpsDid_st}
        return (self.DidFid, self.DhpsDid)

    def _sequence_match_link(
        self,
    ) -> Tuple[Dict[int, int], Dict[Tuple[int, int, bool], int]]:
        design_sequences = [
            "".join([base.seq for base in seg.tour]).upper().replace("N", "T")
            for seg in self.design.strands
        ]

        sequence_pos = (
            1 if self.fit.u.residues[0].resname[0] == "D" else 0
        )  # NOTE: DT,DA,DG,DC vs THY,ADE,GUA,CYT nomenclature
        fit_sequences = [
            "".join([res.resname[sequence_pos] for res in seg.residues])
            for seg in self.fit.u.segments
        ]
        strand_DidFid = [fit_sequences.index(seq) for seq in design_sequences]

        for strand in self.design.design.strands:
            for base_tour_idx, base in enumerate(strand.tour):
                strand_Fid = strand_DidFid[strand.id]
                self.DidFid[base.id] = (
                    self.fit.u.segments[strand_Fid].residues[base_tour_idx].resindex
                )
                self.DhpsDid[(base.h, base.p, base.is_scaf)] = base.id

                self.Dcolor[strand_Fid] = strand.icolor  # TODO: might be wrong

        return (self.DidFid, self.DhpsDid)

    def _identify_bp(self) -> Dict[int, int]:
        """link basepairs by mapping indices according to json (cadnano).
            basepairs are mapped from scaffold to staple, unique (invertable).
        -------
         Returns
            -------
            self.Fbp
                fit-id -> fit-id
        """
        self.Fbp = {
            self.DidFid[base.id]: self.DidFid[base.across.id]
            for base in self.design.scaffold
            if base.across is not None
        }
        return self.Fbp

    def _get_n_strand(self, base: DnaBase, direct: str, steps=1, local=True) -> Optional[DnaBase]:
        """direct = ["up","down"]"""
        if steps == 0:
            return base
        if steps < 0:
            direct = "up" if direct == "down" else "down"
        if (base.h % 2) == 1 and not local:
            direct = "up" if direct == "down" else "down"

        for _ in range(abs(steps)):
            base = base.up if direct == "up" else base.down
            if base is None:
                return None
        return base

    def _get_n_helix(self, base: DnaBase, direct: int, steps=1) -> Optional[DnaBase]:
        """direct = [1,-1]"""
        if steps == 0:
            return base
        helix, position, is_scaf = base.h, base.p, base.is_scaf
        if steps < 0:
            steps = abs(steps)
            direct = -direct
        # first check the number of skips passed
        n_skips = 0
        for n_steps in range(direct, direct * (steps + 1), direct):
            n_position = position + n_steps
            if (helix, n_position, True) not in self.DhpsDid:
                n_skips += 1
        # move one position further if on skip
        n_position = position + direct * (steps + n_skips)
        if (helix, n_position, True) not in self.DhpsDid:
            n_position += direct
        return self.design.Dhps_base.get((helix, n_position, is_scaf), None)

    def is_co(self, base: DnaBase, neighbor: Optional[DnaBase]) -> bool:
        """determine if a base is part of a crossover."""
        if neighbor is None:
            return False
        else:
            return bool(neighbor.h != base.h)

    def _identify_crossover(self) -> None:
        """Affects
        -------
            self.Fco
        """

        def get_crossover_leg(base: Optional[DnaBase], direct: int) -> Optional[DnaBase]:
            if base is None:
                return None
            else:
                return self._get_n_helix(base=base, direct=direct, steps=2)

        def get_crossover(
            bA: DnaBase,
            bC: DnaBase,
            bB: Optional[DnaBase],
            bD: Optional[DnaBase],
            direct: int,
            typ: str,
        ) -> Tuple[str, Crossover]:

            bA_ = get_crossover_leg(base=bA, direct=(-1 * direct))
            bB_ = get_crossover_leg(base=bB, direct=direct)
            bC_ = get_crossover_leg(base=bC, direct=(-1 * direct))
            bD_ = get_crossover_leg(base=bD, direct=direct)

            positional, legs = [], []
            co_pos = []
            for bP, bL in zip([bA, bB, bC, bD], [bA_, bB_, bC_, bD_]):
                hp = None if bP is None else (bP.h, bP.p)
                positional.append(hp)
                if bP is not None:
                    co_pos.append((bP.h, bP.p))
                hp = None if bL is None else (bL.h, bL.p)
                legs.append(bL)
            co = Crossover(
                positional=positional,
                legs=legs,
                typ=typ,
                is_scaf=bA.is_scaf,
            )
            key = str(sorted(co_pos))
            return key, co

        co_subparts = set()
        for base in self.design.allbases:
            for direct in ["up", "down"]:
                neighbor = self._get_n_strand(base, direct)
                if self.is_co(base=base, neighbor=neighbor):
                    co_subparts.add(frozenset([base, neighbor]))
                    break

        while co_subparts:
            bA, bC = co_subparts.pop()
            co_direct = "up" if self._get_n_strand(bA, "up") == bC else "down"
            bN = bA.up if co_direct == "down" else bA.down
            direct_int = bA.p - bN.p
            bB = self._get_n_helix(base=bA, direct=direct_int)
            bD = self._get_n_helix(base=bC, direct=direct_int)
            bBbD = frozenset([bB, bD])

            if bBbD in co_subparts:
                co_subparts.remove(bBbD)
                typ = "full"
            elif (bB is not None) and (bD is not None):
                typ = "half"
            else:
                typ = "end"

            key, crossover = get_crossover(bA=bA, bB=bB, bC=bC, bD=bD, direct=direct_int, typ=typ)
            self.Fco[key] = crossover

    def _identify_nicks(self) -> None:
        """Affects
        -------
            self.Fnicks
        """

        def is_nick(candidate: DnaBase, base: DnaBase) -> bool:
            is_onhelix = candidate.h == base.h
            is_neighbor = abs(base.p - candidate.p) <= 2  # skip = 2
            is_base = candidate is base
            b_Fid = self.DidFid[base.id]
            c_Fid = self.DidFid[candidate.id]
            is_ds = all([(x in self.Fbp.values()) for x in [b_Fid, c_Fid]])
            return all([is_onhelix, is_neighbor, not is_base, is_ds])

        def Fid(Did: int) -> int:
            return self.DidFid[Did]

        start_bases = [s[0] for s in self.design.staples]
        end_bases = [s[-1] for s in self.design.staples]

        self.Fnicks = {
            Fid(start.id): Fid(candi.id)
            for start in start_bases
            for candi in end_bases
            if is_nick(candidate=candi, base=start)
        }

    def write_internal_gridpdb(
        self, dest: Path, exclude_ss=True, exclude_id=None, exclude_resid=None
    ):
        """write custom grid.pdb file for a given set of indices"""
        if not hasattr(self, "link"):
            self.logger.debug("Creating Linkage")
            _ = self.create_linkage()  # initializes and returns link
        u: mda.Universe = self.link.u
        u.add_TopologyAttr("tempfactor")  # init as 0.
        u.add_TopologyAttr("occupancy")  # init as 0.

        atoms_exclude = AtomGroup([], u)
        if exclude_ss:
            for residue in u.residues:
                res_id = residue.resindex
                is_bp_scaffold = res_id in self.Fbp.keys()
                is_bp_staple = res_id in self.Fbp.values()
                if not is_bp_scaffold and not is_bp_staple:
                    atoms_exclude += residue.atoms
        if exclude_id is not None:
            for atom_id in exclude_id:
                atoms_exclude += u.atoms[atom_id]
        if exclude_resid is not None:
            for res_id in exclude_resid:
                atoms_exclude += u.residues[res_id].atoms

        for atom in u.residues.atoms:
            is_H = "H" in atom.name
            is_excluded = atom in atoms_exclude
            if not is_H and not is_excluded:
                atom.tempfactor = atom.mass
                atom.occupancy = 1.0
        u.atoms.write(str(dest), bonds=None)

    def compute_bdna(self) -> BDna:
        return BDna(self.link)
