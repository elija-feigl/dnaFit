#!/usr/bin/env python3#
import sys
import os

import MDAnalysis as mda
import numpy as np
import pickle
import contextlib
import argparse

from pathlib import Path
from itertools import chain
from operator import attrgetter
from typing import List, Set, Dict, Tuple, Optional
from collections import namedtuple
from nanodesign.converters import Converter

DICTS = ["bp", "idid", "hpid", "color",
         "coid", "nicks", "skips", "universe"]
# TODO: -low move or rather make use obsolete


@contextlib.contextmanager
def ignored(*exceptions):
    try:
        yield
    except exceptions:
        pass


class Linker(object):
    """ Linker class
    """

    def __init__(self, project):
        self.project = project
        self.fit = Fit(project)
        self.design = Design(project)
        self.Fbp = None
        self.DidFid = None
        self.DhpsDid = None
        self.Fco = None
        self.Dskips = self._get_skips()
        self.Fnicks = None
        self.FidSeq = self._get_sequence()

    def _get_sequence(self) -> Dict[int, str]:
        FidSeq = {res.resindex: res.resname[0]
                  for res in self.fit.u.residues
                  }
        return FidSeq

    def _is_del(self, base) -> bool:
        return base.num_deletions != 0

    def _get_skips(self) -> List[Tuple[int, int, bool]]:
        """ collect position of skips in design
        -------
            Returns
            -------
            list Dskips
                (base.h, base.p, base.is_scaf) of all skips
        """
        design_allbases = [base
                           for strand in self.design.strands
                           for base in strand.tour
                           ]
        Dskips = [(base.h, base.p, base.is_scaf)
                  for base in design_allbases
                  if self._is_del(base)
                  ]
        return Dskips

    def create_linkage(self):
        Linkage = namedtuple("Linkage", ["Fbp",
                                         "DidFid",
                                         "DhpsDid",
                                         "Dcolor",
                                         "Dskips",
                                         "Fnicks",
                                         "universe"
                                         "Fco",
                                         ]
                             )
        Link = self.link()
        Fco = self.identify_crossover()
        Fnicks = self.identify_nicks()
        universe = self.get_universe_tuple()
        return Linkage(Fbp=Link.Fbp,
                       DidFid=Link.DidFid,
                       DhpsFid=Link.DhpsDid,
                       Dcolor=Link.Dcolor,
                       Dskips=Link.Dskips,
                       Fco=Fco,
                       Fnicks=Fnicks,
                       universe=universe,
                       )

    def _link_scaffold(self) -> (Dict[int, int],
                                 Dict[Tuple[int, int, bool], int],
                                 ):
        """ collect position in scaffold (0-x) by comparing to index in list
            of scaffold_design positions
        -------
            Returns
            -------
            dict DidFid
                design-id (int) -> fit-id (int)
            dict DhpsDid
                helix-number (int), base-position (int), is_scaffold (bool)
                -> design-id (int)
        """
        Dscaffold = self.design.scaffold
        Did = [base.id for base in Dscaffold]
        Dhp = [(base.h, base.p, True) for base in Dscaffold]
        Fid_local = [Did.index(base.id) for base in Dscaffold]
        Fid_global = self.fit.scaffold.residues[Fid_local].resindices

        DidFid = dict(zip(Did, Fid_global))
        DhpsDid = dict(zip(Dhp, Did))
        return DidFid, DhpsDid

    def _link_staples(self) -> (Dict[int, int],
                                Dict[Tuple[int, int, bool], int],
                                Dict[int, int],
                                ):
        """same procedure as scaffold for each
        -------
         Returns
            -------
            dict DidFid
                design-id (int) -> fit-id (int)
            dict DhpsDid
                helix-number (int), base-position (int), is_scaffold (bool)
                -> design-id (int)
            dict color
                fit-segment-id (int) -> color (hex?)
        """
        def get_resid(segindex: int, resindex_local: int) -> int:
            return self.fit.staples[segindex].residues[resindex_local].resindex

        DidFid = {}
        DhpsDid = {}
        color = {}

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

        return DidFid, DhpsDid, color

    def _link_bp(self) -> Dict[int, int]:
        """ link basepairs by mapping indices according to json (cadnano).
            basepairs are mapped from scaffold to staple, unique (invertable).
        -------
         Returns
            -------
            dict Fbp
                fit-id (int) -> fit-id (int)
        """
        def Fid(Did):
            return self.DidFid[Did]

        return {Fid(base.id): Fid(base.across.id)
                for base in self.design.scaffold
                if base.across is not None
                }

    def link(self):
        """ invoke _link_scaffold, _link_staples, _link_bp to compute mapping
            of every base design-id to fit-id as well as the basepair mapping.
            basepairs are mapped from scaffold to staple, unique (invertable).
            updates linker attributes corresponding to the respective mapping
            and returns them.
        -------
         Returns
            -------
            dict self.Fbp
                fit-id (int) -> fit-id (int)
            dict self.DidFid
                design-id (int) -> fit-id (int)
            dict self.DhpsDid
                helix-number (int), base-position (int), is_scaffold (bool)
                -> design-id (int)
            dict self.color
                fit-segment-id (int) -> color (hex?)
        """
        DidFid_sc, DhpsDid_sc = self._link_scaffold()
        DidFid_st, DhpsDid_st, self.Dcolor = self._link_staples()

        self.DidFid = {**DidFid_sc, **DidFid_st}
        self.DhpsDid = {**DhpsDid_sc, **DhpsDid_st}
        self.Fbp = self._link_bp()
        Link = namedtuple("Link", ["Fbp",
                                   "DidFid",
                                   "DhpsDid",
                                   "Dcolor",
                                   "Dskips",
                                   ]
                          )
        return Link(Fbp=self.Fbp,
                    DidFid=self.DidFid,
                    DhpsDid=self.DhpsDid,
                    Dcolor=self.Dcolor,
                    Dskips=self.Dskips,
                    )

    def _get_nextInHelix(self, h, p, is_scaf, i):
        Dhps = (h, p+i, is_scaf)
        while Dhps in self.Dskips:
            i += np.sign(i)
            Dhps = (h, p+i, is_scaf)
        return Dhps

    def identify_crossover(self) -> Dict:
        """ for every base id that is involved in a crossover
            updates linker attribute of crossovers and returns it
        -------
         Returns
            -------
            dict self.Fco
                fit-id (int) ->
                "co_index": co_index (int), "co": co_id (int),
                "leg": leg_id (int),
                "is_scaffold": is_scaf (bool),
                "position": (h, p, is_scaf) /(int,int, bool),
                if end, single, dpouble
                    "type": ("end", None)
                    "type": ("single", single, single_leg) /(str, int, int)
                    "type": ("double", double) /(str, int)
        """
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

        def get_co_leg_id(base, direct: str) -> int:
            """determine leg base base and up/down (def: 2 bases away)"""
            l = base.down.p if direct == "up" else base.up.p
            i = (l - base.p) * 2.
            Dhps = self._get_nextInHelix(base.h, base.p, base.is_scaf, i)
            leg_Did = self.DhpsDid[Dhps]
            return self.DidFid[leg_Did]

        def get_next(base, direct: str):
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

        def is_co(base, neighbor, direct: str) -> bool:
            if neighbor is None:
                return False
            else:
                while self._is_del(neighbor):
                    neighbor = get_next(neighbor, direct)
                return neighbor.h != base.h

        def get_co(base, neighbor, direct: str, run_id: int) -> (Dict, int):
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

        self.Fco = {}
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

    def identify_nicks(self) -> Dict[int, int]:
        """ for every nick, provide id of base accross nick
        -------
         Returns
            -------
            dict self.Fnicks
                fit-id (int) -> fit_id (int)
        """
        def is_nick(candidate, base) -> bool:
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

    def get_universe_tuple(self) -> Tuple[str, str]:
        infile = self.project.input / self.project.name
        top = infile.with_suffix(".psf")
        trj = infile.with_suffix(".dcd")
        return (str(top), str(trj))


class Fit(object):

    def __init__(self, project):
        self.infile = project.input / project.name
        self.u = self._get_universe()
        self.scaffold, self.staples = self._split_strands()

    def _get_universe(self):
        top = self.infile.with_suffix(".psf")
        trj = self.infile.with_suffix(".dcd")
        # TODO: -mid- if pdb, try invoke vmd animate write dcd
        if top.exists() and trj.exists():
            u = mda.Universe(str(top), str(trj))
        else:
            raise FileNotFoundError
        return u

    def _split_strands(self):
        # TODO: -low- multiscaffold
        strands = self.u.segments
        scaffold = max(strands, key=attrgetter("residues.n_residues"))
        staples = [strand for strand in strands if len(
            strand.atoms) != len(scaffold.atoms)]
        return scaffold, staples


class Design(object):

    def __init__(self, project):
        self.infile = project.input / project.name
        self.design = self._get_design()
        self.strands = self.design.strands
        self.scaffold = self._clean_scaffold(self.strands)
        self.excl = self.scaffold[0].strand
        self.staples = self._clean_staple(self.strands)
        self.helixorder = self._create_helix_order()
        self.stapleorder = self._create_staple_order()

    def _is_del(self, base) -> bool:
        return base.num_deletions != 0

    def _clean_scaffold(self, strands):
        # TODO: -low- multiscaffold
        # TODO: -low- insertions
        scaffold = [s.tour for s in strands if s.is_scaffold][0]
        scaffold_clean = [b for b in scaffold if not self._is_del(b)]
        return scaffold_clean

    def _clean_staple(self, strands):
        # TODO: -low- insertions
        staples = [s.tour for s in strands if not s.is_scaffold]
        staples_clean = [[b for b in s if not self._is_del(b)]
                         for s in staples
                         ]
        return staples_clean

    def _get_design(self):
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
            dict helixorder
                (int) -> (int)
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
            dict stapleorder
                (int) -> (int)
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


def get_description():
    return """links structural information of the cadnano designfile
              [design.json] to fitted atomic model [design.psf, design.dcd].
              stores dictionaries as pickles containg mapping for motifs,
              residue-id, lattice position and base-pairing."""


def proc_input():
    parser = argparse.ArgumentParser(
        description=get_description(),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument("--folder",
                        help="input folder",
                        type=str,
                        default="./",
                        )
    parser.add_argument("--name",
                        help="name of design and files",
                        type=str,
                        required="True",
                        default=argparse.SUPPRESS,
                        )
    args = parser.parse_args()
    Project = namedtuple("Project", ["input", "output", "name"])
    project = Project(input=Path(args.folder),
                      output=Path(args.folder) / "analysis",
                      name=args.name,
                      )
    with ignored(FileExistsError):
        os.mkdir(project.output)
    return project


def main():

    # process input
    project = proc_input()
    print("output to ", project.output)
    linker = Linker(project)
    linkage = linker.create_linkage()

    for name in linkage:
        pickle.dump(eval(name), open(
            str(project.output) + "__" + name + ".p", "wb"))


if __name__ == "__main__":
    main()
