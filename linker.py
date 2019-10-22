#!/usr/bin/env python3#
import numpy as np
import pickle
import attr
import nanodesign as nd

from itertools import chain
from typing import Set, Dict, Tuple, Any

from project import Project
from fit import Fit
from design import Design


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

    def _eval_skips(self) -> Set[Tuple[int, int, bool]]:
        """ Returns
            -------
                (base.h, base.p, base.is_scaf) of all skips
        """
        design_allbases = [base
                           for strand in self.design.strands
                           for base in strand.tour
                           ]
        Dskips = set((base.h, base.p, base.is_scaf)
                     for base in design_allbases
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
        -------
         Returns
            -------
                Linkage
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

    def _get_nextInHelix(self, h: int, p: int, is_scaf: bool, i: int
                         ) -> Tuple[int, int, bool]:
        Dhps = (h, p + i, is_scaf)
        while Dhps in self.Dskips:
            i += np.sign(i)
            Dhps = (h, p + i, is_scaf)
        return Dhps

    def _identify_crossover(self) -> Dict[int, Any]:
        # TODO: make crossover object 
        """ for every base id that is involved in a crossover
            updates linker attribute of crossovers and returns it
        -------
         Returns
            -------
            self.Fco
                fit-id ->
                "co_index": co_index, "co": co_id,
                "leg": leg_id,
                "is_scaffold": is_scaf,
                "position": (h, p, is_scaf)
                if end, single, dpouble
                    "type": ("end", None)
                    "type": ("single", single, single_leg)
                    "type": ("double", double)
        """
        # TODO: double, single -> full, half
        def add_co_type() -> None:

            def get_co_leg_p(n_Dhps: Tuple[int, int, bool], i: int) -> int:
                """determine leg base by Dhps and +1/-1 (def: 2 bases away)"""
                h, p, is_scaf = n_Dhps
                Dhps = self._get_nextInHelix(h, p, is_scaf, i)
                leg_Did = self.DhpsDid[Dhps]
                return self.DidFid[leg_Did]

            def is_nextInStrand(b_Dhps: Tuple[int, int, bool],
                                a_Dhps: Tuple[int, int, bool]
                                ) -> bool:
                a_Fid = self.DidFid[self.DhpsDid[a_Dhps]]
                b_Fid = self.DidFid[self.DhpsDid[b_Dhps]]
                return True if abs(a_Fid - b_Fid) == 1 else False

            for co in iter(self.Fco.values()):
                h, p, is_scaf = co["position"]
                for i in [-1, 1]:
                    n_Dhps = self._get_nextInHelix(h, p, is_scaf, i)
                    n_Did = self.DhpsDid.get(n_Dhps, None)
                    n_Fid = self.DidFid.get(n_Did, None)

                    if n_Fid is None:
                        co["type"] = ("end", None, None)
                    elif is_nextInStrand(co["position"], n_Dhps):
                        continue
                    else:
                        is_double = (n_Fid in iter(self.Fco.keys()))
                        if is_double:
                            co["type"] = ("double", n_Fid, None)
                        else:
                            leg = get_co_leg_p(n_Dhps, i)
                            co["type"] = ("single", n_Fid, leg)
                co.pop("base")
            return

        def get_co_leg_id(base: "nd.base", direct: str) -> int:
            """determine leg base base and up/down (def: 2 bases away)"""
            m = base.down.p if direct == "up" else base.up.p
            i = (m - base.p) * 2.
            Dhps = self._get_nextInHelix(base.h, base.p, base.is_scaf, i)
            leg_Did = self.DhpsDid[Dhps]
            return self.DidFid[leg_Did]

        def get_next(base: "nd.base", direct: str) -> "nd.base":
            n = (base.up if direct == "up" else base.down)
            if n is None:
                return None
            else:
                n_position = (n.h, n.p, n.is_scaf)
                while n_position in self.Dskips:
                    n = (n.up if direct == "up" else n.down)
                    if n is None:
                        return None
                    n_position = (n.h, n.p, n.is_scaf)
                return n

        def is_co(base: "nd.base", neighbor: "nd.base", direct: str) -> bool:
            if neighbor is None:
                return False
            else:
                while self._is_del(neighbor):
                    neighbor = get_next(neighbor, direct)
                return neighbor.h != base.h

        def get_co(base: "nd.base",
                   neighbor: "nd.base",
                   direct: str,
                   run_id: int
                   ) -> Tuple[dict, int]:
            co_id = self.DidFid[neighbor.id]
            leg_id = get_co_leg_id(base, direct)
            Dhps = (base.h, base.p, base.is_scaf)

            try:  # TODO -high cleanup co_index
                index = self.Fco[co_id]["co_index"]
            except KeyError:
                run_id += 1
                index = run_id
            data = {"co_index": index,
                    "co": co_id,
                    "leg": leg_id,
                    "is_scaffold": base.is_scaf,
                    "position": Dhps,
                    "base": base,
                    }
            return data, run_id

        allbases = self.design.scaffold.copy()
        staples = [base for staple in self.design.staples for base in staple]
        allbases.extend(staples)

        run_id = 0
        for base in allbases:
            base_Fid = self.DidFid[base.id]
            for direct in ["up", "down"]:
                neighbor = get_next(base, direct)
                if is_co(base, neighbor, direct):
                    co_data, run_id = get_co(base, neighbor, direct, run_id)
                    self.Fco[base_Fid] = co_data
                    break
        add_co_type()
        return self.Fco

    def _identify_nicks(self) -> Dict[int, int]:
        """ for every nick, provide id of base accross nick, bidirectional
        -------
         Returns
            -------
            self.Fnicks
                fit-id -> fit_id
        """
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
