#!/usr/bin/env python3#
import attr
import nanodesign as nd

from typing import Set, Dict, Tuple, Optional

from project import Project
from fit import Fit
from design import Design
from crossover import Crossover
from linkage import Linkage
from basepair import BasePair


@attr.s
class Linker(object):
    """ Linker class
    """
    # TODO: move categorize to linker?
    project: Project = attr.ib()

    Fbp: Dict[int, int] = dict()
    DidFid: Dict[int, int] = dict()
    DhpsDid: Dict[Tuple[int, int, bool], int] = dict()
    Fnicks: Dict[int, int] = dict()
    FidSeq: Dict[int, str] = dict()
    Dskips: Set[Tuple[int, int]] = set()
    Fco: Dict[str, Crossover] = dict()

    def __attrs_post_init__(self) -> None:
        self.fit: Fit = Fit(self.project)
        self.design: Design = Design(self.project)

    def _eval_sequence(self) -> Dict[int, str]:
        FidSeq = {res.resindex: res.resname[0]
                  for res in self.fit.u.residues
                  }
        return FidSeq

    def _eval_skips(self) -> Set[Tuple[int, int]]:
        Dskips = set((base.h, base.p)
                     for base in self.design.allbases
                     if self._is_del(base)
                     )
        return Dskips

    def _is_del(self, base: "nd.base") -> bool:
        return base.num_deletions != 0

    def create_linkage(self) -> Linkage:
        """ invoke _link_scaffold, _link_staples, _link_bp to compute mapping
            of every base design-id to fit-id as well as the basepair mapping.
            basepairs are mapped from scaffold to staple, unique (invertable).
            updates linker attributes corresponding to the respective mapping
            and returns them.
        """
        self.Dskips = self._eval_skips()
        self.FidSeq = self._eval_sequence()
        self._link()
        self._identify_bp()
        self._identify_crossover()
        self._identify_nicks()
        self.link = Linkage(Fbp=self.Fbp,
                            DidFid=self.DidFid,
                            DhpsDid=self.DhpsDid,
                            Dcolor=self.Dcolor,
                            Dskips=self.Dskips,
                            Fco=self.Fco,
                            Fnicks=self.Fnicks,
                            FidSeq=self.FidSeq,
                            u=self.fit.u,
                            )
        return self.link

    def _link(self) -> Tuple[Dict[int, int],
                             Dict[Tuple[int, int, bool], int],
                             ]:
        def link_scaffold() -> Tuple[Dict[int, int],
                                     Dict[Tuple[int, int, bool], int],
                                     ]:
            """ collect position in scaffold (0-x) by comparing index in list
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
            Dhp = [(base.h, base.p, True) for base in Dscaffold]
            Fid_local = [Did.index(base.id) for base in Dscaffold]
            Fid_global = self.fit.scaffold.residues[Fid_local].resindices

            DidFid = dict(zip(Did, Fid_global))
            DhpsDid = dict(zip(Dhp, Did))
            return (DidFid, DhpsDid)

        def link_staples() -> Tuple[Dict[int, int],
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
            def get_resid(segindex: int, resindex_local: int) -> int:
                segment = self.fit.staples[segindex]
                return segment.residues[resindex_local].resindex

            DidFid: Dict[int, int] = {}
            DhpsDid: Dict[Tuple[int, int, bool], int] = {}
            color: Dict[int, int] = {}

            for i, staple in enumerate(self.design.staples):
                seg_id = self.design.stapleorder[i]

                Did = [base.id for base in staple]
                Dhp = [(base.h, base.p, False) for base in staple]

                Fid_local = [Did.index(base.id)for base in staple]
                Fid_global = [get_resid(seg_id, resid) for resid in Fid_local]

                icolor = self.design.design.strands[staple[0].strand].icolor
                segidxforcolor = self.fit.staples[seg_id].segindex
                color[segidxforcolor] = icolor

                DidFid_add = dict(zip(Did, Fid_global))
                DhpsDid_add = dict(zip(Dhp, Did))
                DidFid = {**DidFid, **DidFid_add}
                DhpsDid = {**DhpsDid, **DhpsDid_add}

            return (DidFid, DhpsDid, color)

        DidFid_sc, DhpsDid_sc = link_scaffold()
        DidFid_st, DhpsDid_st, self.Dcolor = link_staples()

        self.DidFid = {**DidFid_sc, **DidFid_st}
        self.DhpsDid = {**DhpsDid_sc, **DhpsDid_st}
        return (self.DidFid, self.DhpsDid)

    def _identify_bp(self) -> Dict[int, int]:
        """ link basepairs by mapping indices according to json (cadnano).
            basepairs are mapped from scaffold to staple, unique (invertable).
        -------
         Returns
            -------
            self.Fbp
                fit-id -> fit-id
        """
        self.Fbp = {self.DidFid[base.id]: self.DidFid[base.across.id]
                    for base in self.design.scaffold
                    if base.across is not None
                    }
        return self.Fbp

    def _get_n_strand(self, base: "nd.base", direct: str, steps=1
                      ) -> Optional["nd.base"]:
        """direct = ["up","down"]"""
        if steps == 0:
            return base
        if steps < 0:
            direct = "up" if direct == "down" else "down"

        for _ in range(steps):
            base = (base.up if direct == "up" else base.down)
            if base is None:
                return None
            else:
                n_position = (base.h, base.p)
                while n_position in self.Dskips:
                    base = (base.up if direct == "up" else base.down)
                    if base is None:
                        return None
                    n_position = (base.h, base.p)
        return base

    def _get_n_helix(self, base: "nd.base", direct: int, steps=1
                     ) -> Optional["nd.base"]:
        """direct = [1,-1]"""
        if steps == 0:
            return base
        helix, position, is_scaf = base.h, base.p, base.is_scaf

        if steps < 0:
            steps = abs(steps)
            direct = -direct

        n_skips = 0
        for n in range(direct, direct * (steps + 1), direct):
            n_position = position + n
            if (helix, n_position) in self.Dskips:
                n_skips += 1

        n_position = position + direct * (steps + n_skips)
        while (helix, n_position) in self.Dskips:
            n_position += direct

        return self.design.Dhps_base.get((helix, n_position, is_scaf), None)

    def _get_bp(self, base: "nd.residue") -> Optional[BasePair]:
        if base is None:
            return None
        else:
            resindex = self.DidFid[base.id]
            Fbp_all = {**self.Fbp, **{v: k for k, v in iter(self.Fbp.items())}}
            wcindex = Fbp_all.get(resindex, None)
            sc_index = resindex if base.is_scaf else wcindex
            st_index = wcindex if base.is_scaf else resindex

            sc = None if sc_index is None else self.fit.u.residues[sc_index]
            st = None if st_index is None else self.fit.u.residues[st_index]
            hp = (base.h, base.p)
            return BasePair(sc=sc, st=st, hp=hp)

    def _identify_crossover(self) -> None:
        """ Affects
            -------
                self.Fco
        """
        def get_co_leg(base: Optional["nd.base"], direct: int
                       ) -> Optional["nd.base"]:
            if base is None:
                return None
            else:
                return self._get_n_helix(base=base, direct=direct, steps=2)

        def is_co(base: "nd.base", neighbor: Optional["nd.base"], direct: str
                  ) -> bool:
            if neighbor is None:
                return False
            else:
                while self._is_del(neighbor):
                    neighbor = self._get_n_strand(neighbor, direct)
                return neighbor.h != base.h

        def get_co(bA: "nd.base",
                   bC: "nd.base",
                   bB: Optional["nd.base"],
                   bD: Optional["nd.base"],
                   direct: int,
                   typ: str,
                   ) -> Tuple[str, Crossover]:

            bA_ = get_co_leg(base=bA, direct=(-1 * direct))
            bB_ = get_co_leg(base=bB, direct=direct)
            bC_ = get_co_leg(base=bC, direct=(-1 * direct))
            bD_ = get_co_leg(base=bD, direct=direct)

            Ps, Ls = list(), list()
            co_pos = list()
            for bP, bL in zip([bA, bB, bC, bD], [bA_, bB_, bC_, bD_]):
                Ps.append(self._get_bp(base=bP))
                if bP is not None:
                    co_pos.append((bP.h, bP.p))
                Ls.append(self._get_bp(base=bL))
            co = Crossover(Ps=Ps,
                           Ls=Ls,
                           typ=typ,
                           is_scaf=bA.is_scaf,
                           )
            key = str(sorted(co_pos))
            return key, co

        co_subparts = set()
        for base in self.design.allbases:
            if self._is_del(base):
                continue
            for direct in ["up", "down"]:
                neighbor = self._get_n_strand(base, direct)
                if is_co(base, neighbor, direct):
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
                typ = "full"
            elif (bB is not None) and (bD is not None):
                typ = "half"
            else:
                typ = "end"

            key, co = get_co(bA=bA, bB=bB, bC=bC, bD=bD,
                             direct=direct_int, typ=typ)
            self.Fco[key] = co

    def _identify_nicks(self) -> None:
        """ Affects
            -------
                self.Fnicks
        """
        def is_nick(candidate: "nd.base", base: "nd.base") -> bool:
            is_onhelix = (candidate.h == base.h)
            is_neighbor = (abs(base.p - candidate.p) <= 2)  # skip = 2
            is_base = (candidate is base)
            b_Fid = self.DidFid[base.id]
            c_Fid = self.DidFid[candidate.id]
            is_ds = all([(x in self.Fbp.values()) for x in [b_Fid, c_Fid]])
            return all([is_onhelix, is_neighbor, not is_base, is_ds])

        def Fid(Did: int) -> int:
            return self.DidFid[Did]

        start_bases = [s[0] for s in self.design.staples]
        end_bases = [s[-1] for s in self.design.staples]

        self.Fnicks = {Fid(start.id): Fid(candi.id)
                       for start in start_bases
                       for candi in end_bases
                       if is_nick(candidate=candi, base=start)
                       }
