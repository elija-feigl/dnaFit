#!/usr/bin/env python3#
import pickle
import attr
import nanodesign as nd

from itertools import chain
from typing import Set, Dict, Tuple, Any, Optional, List

from project import Project
from fit import Fit
from design import Design


@attr.s
class Crossover(object):
    Ps: List[Tuple[int, int]] = attr.ib()
    Ls: List[Tuple[int, int]] = attr.ib()
    P_pos: List[Tuple[int, int]] = attr.ib()
    L_pos: List[Tuple[int, int]] = attr.ib()
    typ: str = attr.ib()
    is_scaf: bool = attr.ib()

    def transform2bp(self, bps: dict):
        self.Ps, self.Ls = list(), list()
        for hp, hp_ in zip(self.P_pos, self.L_pos):
            self.Ps.append(bps[hp])
            self.Ls.append(bps[hp_])


@attr.s(auto_attribs=True)
class Linkage(object):
    Fbp: Dict[int, int] = {}
    DidFid: Dict[int, int] = {}
    DhpsDid: Dict[Tuple[int, int, bool], int] = {}
    Dcolor: Dict[int, int] = {}
    Dskips: Set[Tuple[int, int, bool]] = set()
    Fnicks: Dict[int, int] = {}
    FidSeq: Dict[int, str] = {}
    Fco: Dict[int, Any] = {}
    universe: Tuple[str, str] = ("", "")

    def dump_linkage(self, project: Project) -> None:
        for name, link in vars(self).items():
            output = project.output / "{}__{}.p".format(project.name, name)
            pickle.dump(link, open(output, "wb"))
        return

    def load_linkage(self, project: Project) -> None:
        for name in vars(self).keys():
            input = project.output / "{}__{}.p".format(project.name, name)
            value = pickle.load(open(input, "rb"))
            setattr(self, name, value)
        return

    def reverse(self) -> None:
        def reverse_d(dict: dict) -> dict:
            return {v: k for k, v in iter(dict.items())}

        self.FidDid = reverse_d(self.DidFid)
        self.DidDhps = reverse_d(self.DhpsDid)
        self.Fbp_rev = reverse_d(self.Fbp)
        self.Fbp_full = {**self.Fbp, **self.Fbp_rev}
        return


@attr.s
class Linker(object):
    """ Linker class
    """
    # TODO: move categorize to linker?
    project: Project = attr.ib()
    Fbp: Dict[int, int] = {}
    DidFid: Dict[int, int] = {}
    DhpsDid: Dict[Tuple[int, int, bool], int] = {}
    Fnicks: Dict[int, int] = {}
    FidSeq: Dict[int, str] = {}
    Dskips: Set[Tuple[int, int, bool]] = set()
    Fco: Dict[int, Any] = {}

    def __attrs_post_init__(self):
        self.fit: Fit = Fit(self.project)
        self.design: Design = Design(self.project)
        self.Dskips = self._eval_skips()
        self.FidSeq = self._eval_sequence()

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

        self._link()
        self._identify_bp()
        self._identify_crossover()
        self._identify_nicks()
        universe = self._get_universe_tuple()

        return Linkage(Fbp=self.Fbp,
                       DidFid=self.DidFid,
                       DhpsDid=self.DhpsDid,
                       Dcolor=self.Dcolor,
                       Dskips=self.Dskips,
                       Fco=self.Fco,
                       Fnicks=self.Fnicks,
                       FidSeq=self.FidSeq,
                       universe=universe,
                       )

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
        def Fid(Did: int) -> int:
            return self.DidFid[Did]

        self.Fbp = {Fid(base.id): Fid(base.across.id)
                    for base in self.design.scaffold
                    if base.across is not None
                    }
        return self.Fbp

    def _get_n_strand(self, base: "nd.base", direct: str, steps=1
                      ) -> Optional["nd.base"]:
        """direct = ["up","down"]"""
        if steps == 0:
            return base

        # TODO: steps missing
        n = (base.up if direct == "up" else base.down)
        if n is None:
            return None
        else:
            n_position = (n.h, n.p)
            while n_position in self.Dskips:
                n = (n.up if direct == "up" else n.down)
                if n is None:
                    return None
                n_position = (n.h, n.p)
            return n

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
        # check if we passed skips
        for n in range(direct, direct * (steps + 1), direct):
            n_position = position + n
            if (helix, n_position) in self.Dskips:
                n_skips += 1

        # check if land on skip
        n_position = position + direct * (steps + n_skips)
        while (helix, n_position) in self.Dskips:
            n_position += direct

        if (helix, n_position, is_scaf) in self.design.Dhps_base.keys():
            return self.design.Dhps_base[(helix, n_position, is_scaf)]
        else:
            return None

    def _get_bp_tuple(self, base: Optional["nd.residue"]
                      ) -> Optional[Tuple[Tuple[int, int], Tuple[int, int]]]:
        if base is None:
            return None, None

        h, p, is_scaf = base.h, base.p, base.is_scaf
        resindex = self.DidFid[base.id]
        Fbp_full = {**self.Fbp, **{v: k for k, v in iter(self.Fbp.items())}}
        if resindex in Fbp_full:
            wcindex = Fbp_full[resindex]
        else:
            wcindex = None

        if is_scaf:
            sc_index, st_index = resindex, wcindex
        else:
            sc_index, st_index = wcindex, resindex

        bp_resindices = (sc_index, st_index)
        pos = (h, p)

        return pos, bp_resindices

    def _identify_crossover(self) -> Dict[int, Any]:
        """ for every base id that is involved in a crossover
            updates linker attribute of crossovers and returns it
        """
        def get_co_leg(base: Optional["nd.base"], direct: str
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
                   ) -> Crossover:

            bA_ = get_co_leg(base=bA, direct=(-1 * direct))
            bB_ = get_co_leg(base=bB, direct=direct)
            bC_ = get_co_leg(base=bC, direct=(-1 * direct))
            bD_ = get_co_leg(base=bD, direct=direct)

            Ps, Ls, P_pos, L_pos = [], [], [], []
            for bP, bL in zip([bA, bB, bC, bD], [bA_, bB_, bC_, bD_]):
                posP, P_index = self._get_bp_tuple(base=bP)
                Ps.append(P_index)
                P_pos.append(posP)
                posL, L_index = self._get_bp_tuple(base=bL)
                Ls.append(L_index)
                L_pos.append(posL)

            co = Crossover(Ps=Ps,
                           Ls=Ls,
                           P_pos=P_pos,
                           L_pos=L_pos,
                           typ=typ,
                           is_scaf=bA.is_scaf,
                           )
            key = str(sorted(filter(None, co.P_pos)))
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
            co_direct = "up" if bA.up == bC else "down"
            bN = bA.up if co_direct == "down" else bA.down
            direct = bA.p - bN.p
            if direct not in [-1, 1]:
                import ipdb; ipdb.set_trace()
            bB = self._get_n_helix(base=bA, direct=direct)
            bD = self._get_n_helix(base=bC, direct=direct)
            bBbD = frozenset([bB, bD])

            if bBbD in co_subparts:
                typ = "full"
            elif (bB is not None) and (bD is not None):
                typ = "half"
            else:
                typ = "end"

            key, co = get_co(bA=bA, bB=bB, bC=bC, bD=bD,
                             direct=direct, typ=typ)
            self.Fco[key] = co
        return self.Fco

    def _identify_nicks(self) -> Dict[int, int]:  # slow!
        """ for every nick, provide id of base accross nick, bidirectional
        -------
         Returns
            -------
            self.Fnicks
                fit-id -> fit_id
        """
        # called 145924 times
        def is_nick(candidate: "nd.base", base: "nd.base") -> bool:
            is_onhelix = (candidate.h == base.h)
            is_neighbor = (abs(base.p - candidate.p) <= 2)  # skip = 2
            is_base = (candidate is base)
            b_Fid = self.DidFid[base.id]
            c_Fid = self.DidFid[candidate.id]
            is_ds_staple = (b_Fid in self.Fbp.values() and
                            c_Fid in self.Fbp.values()
                            )
            return (is_onhelix and
                    is_neighbor and
                    not is_base and
                    is_ds_staple
                    )

        def Fid(Did: int) -> int:
            return self.DidFid[Did]

        staple_end_bases = list(chain.from_iterable((s[0], s[-1])
                                for s in self.design.staples)
                                )
        self.Fnicks = {Fid(base.id): Fid(candidate.id)
                       for base in staple_end_bases
                       for candidate in staple_end_bases
                       if is_nick(candidate, base)
                       }
        return self.Fnicks

    def _get_universe_tuple(self) -> Tuple[str, str]:
        infile = self.project.input / self.project.name
        top = infile.with_suffix(".psf")
        trj = infile.with_suffix(".dcd")
        return (str(top), str(trj))
