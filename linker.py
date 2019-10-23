#!/usr/bin/env python3#
import pickle
import attr
import nanodesign as nd

from itertools import chain
from typing import Set, Dict, Tuple, Any, Optional, List

from project import Project
from fit import Fit
from design import Design
from basepair import BasePair


@attr.s
class Crossover(object):
    A: BasePair = attr.ib()
    A_: BasePair = attr.ib()
    B: BasePair = attr.ib()
    B_: BasePair = attr.ib()
    C: BasePair = attr.ib()
    C_: BasePair = attr.ib()
    D: BasePair = attr.ib()
    D_: BasePair = attr.ib()
    typ: str = attr.ib()
    is_scaf: bool = attr.ib()

    def __attrs_post_init__(self):
        self.Ps: List[BasePair] = list()
        self.Ls: List[BasePair] = list()
        for P in [self.A, self.B, self.C, self.D]:
            self.Ps.append(P)
        for L in [self.A_, self.B_, self.C_, self.D_]:
            self.Ls.append(L)

        self.position = sorted([P.hp for P in self.Ps if P is not None])


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
        """direct = [1,-1]"""
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
        """direct = ["up","down"]"""
        if steps == 0:
            return base
        helix, position, is_scaf = base.h, base.p, base.is_scaf
        if ((helix % 2) == 1) and is_scaf:
            direct = -direct
        elif ((helix % 2) == 0) and not is_scaf:
            direct = -direct

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

    def _get_bp(self, base: Optional["nd.residue"]
                ) -> Optional[BasePair]:
        if base is None:
            return None

        h, p, is_scaf = base.h, base.p, base.is_scaf
        resindex = self.DidFid[base.id]
        res = self.fit.u.residues[resindex]
        Fbp_full = {**self.Fbp, **{v: k for k, v in iter(self.Fbp.items())}}
        if resindex in Fbp_full:
            wcindex = Fbp_full[resindex]
            wc = self.fit.u.residues[wcindex]
        else:
            wcindex = None
            wc = None

        if is_scaf:
            sc, st = res, wc
        else:
            sc, st = wc, res

        bp = BasePair(sc=sc, st=st, hp=(h, p))
        return bp

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
            bC_ = get_co_leg(base=bC, direct=direct)
            bD_ = get_co_leg(base=bD, direct=(-1 * direct))

            co = Crossover(A=self._get_bp(base=bA),
                           A_=self._get_bp(base=bA_),
                           B=self._get_bp(base=bB),
                           B_=self._get_bp(base=bB_),
                           C=self._get_bp(base=bC),
                           C_=self._get_bp(base=bC_),
                           D=self._get_bp(base=bD),
                           D_=self._get_bp(base=bD_),
                           typ=typ,
                           is_scaf=bA.is_scaf,
                           )
            key = str(co.position)
            return key, co

        dir_str2int = {"up": -1, "down": 1}
        co_subparts = set()
        co_cubpart_dir = dict()
        for base in self.design.allbases:
            if self._is_del(base):
                continue
            for direct in ["up", "down"]:
                neighbor = self._get_n_strand(base, direct)
                if is_co(base, neighbor, direct):
                    AC = frozenset([base, neighbor])
                    co_subparts.add(AC)
                    co_cubpart_dir[AC] = direct
                    break

        while co_subparts:
            bAbC = co_subparts.pop()
            dir_int = dir_str2int[co_cubpart_dir[bAbC]]
            bA, bC = bAbC
            bB = self._get_n_helix(base=bA, direct=dir_int)
            bD = self._get_n_helix(base=bC, direct=(-1 * dir_int))
            bBbD = frozenset([bB, bD])

            if bBbD in co_subparts:
                typ = "full"
            elif (bB is not None) and (bD is not None):
                typ = "half"
            else:
                typ = "end"

            key, co = get_co(bA=bA, bB=bB, bC=bC, bD=bD,
                             direct=dir_int, typ=typ)
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
