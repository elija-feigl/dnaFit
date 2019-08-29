#!/usr/bin/env python3#
import sys
import os

import MDAnalysis as mda
import numpy as np
import pickle

from nanodesign.converters import Converter
from operator import attrgetter
from typing import List, Set, Dict, Tuple, Optional

DICTS = ["dict_bp", "dict_idid", "dict_hpid", "dict_color",
         "dict_coid", "dict_nicks", "list_skips", "universe"]
# TODO: -low move or rather make use obsolete


class Linker(object):
    """ Linker class
    """

    def __init__(self, path: str):
        self.path = path
        self.fit = Fit(self.path)
        self.design = Design(self.path)
        self.dict_Fbp = None
        self.dict_DidFid = None
        self.dict_DhpsDid = None
        self.dict_Fco = None
        self.list_Dskips = self._get_skips()
        self.dict_Fnicks = None
        self.dict_FidSeq = self._get_sequence()

    def _get_sequence(self) -> Dict[int, str]:
        dict_FidSeq = {}
        for res in self.fit.u.residues:
            dict_FidSeq[res.resindex] = res.resname[0]
        return dict_FidSeq

    def _is_del(self, base) -> bool:
        return base.num_deletions != 0

    def _get_skips(self) -> List[Tuple[int, int, bool]]:
        """ collect position of skips in design
        -------
            Returns
            -------
            list list_Dskips
                (base.h, base.p, base.is_scaf) of all skips
        """
        design_allbases = [
            base for strand in self.design.strands for base in strand.tour]
        list_Dskips = []
        for base in design_allbases:
            if self._is_del(base):
                Dhps = (base.h, base.p, base.is_scaf)
                list_Dskips.append(Dhps)

        return list_Dskips

    def _link_scaffold(self) -> (Dict[int, int], 
                                 Dict[Tuple[int, int, bool], int],
                                 ):
        """ collect position in scaffold (0-x) by comparing to index in list
            of scaffold_design positions
        -------
            Returns
            -------
            dict dict_DidFid
                design-id (int) -> fit-id (int)
            dict d_DhpsDid
                helix-number (int), base-position (int), is_scaffold (bool)
                -> design-id (int)
        """

        D_ids = [base.id for base in self.design.scaffold]
        D_hp = [(base.h, base.p, True) for base in self.design.scaffold]

        F_id_local = []
        for base in self.design.scaffold:
            indx = D_ids.index(base.id)
            F_id_local.append(indx)

        F_id_global = self.fit.scaffold.residues[F_id_local].resindices

        dict_DidFid = dict(zip(D_ids, F_id_global))
        dict_DhpsDid = dict(zip(D_hp, D_ids))
        return dict_DidFid, dict_DhpsDid

    def _link_staples(self) -> (Dict[int, int],
                                Dict[Tuple[int, int, bool], int],
                                Dict[int, int],
                                ):
        """same procedure as scaffold for each
        -------
         Returns
            -------
            dict dict_DidFid
                design-id (int) -> fit-id (int)
            dict dict_DhpsDid
                helix-number (int), base-position (int), is_scaffold (bool)
                -> design-id (int)
            dict dict_color
                fit-segment-id (int) -> color (hex?)
        """
        def get_resid(segindex: int, resindex_local: int) -> int:
            return self.fit.staples[segindex].residues[resindex_local].resindex

        dict_DidFid = {}
        dict_DhpsDid = {}
        dict_color = {}

        for i, staple in enumerate(self.design.staples):
            seg_id = self.design.dict_stapleorder[i]

            D_ids = [base.id for base in staple]
            D_hp = [(base.h, base.p, False) for base in staple]

            F_id_local = [D_ids.index(base.id)for base in staple]
            F_id_global = [get_resid(seg_id, resid) for resid in F_id_local]

            color = self.design.design.strands[staple[0].strand].icolor
            segidxforcolor = self.fit.staples[seg_id].segindex
            dict_color[segidxforcolor] = color

            dict_DidFid_add = dict(zip(D_ids, F_id_global))
            dict_DhpsDid_add = dict(zip(D_hp, D_ids))
            dict_DidFid = {**dict_DidFid, **dict_DidFid_add}
            dict_DhpsDid = {**dict_DhpsDid, **dict_DhpsDid_add}

        return dict_DidFid, dict_DhpsDid, dict_color

    def _link_bp(self) -> Dict[int, int]:
        """ link basepairs by mapping indices according to json (cadnano).
            basepairs are mapped from scaffold to staple, unique (invertable).
        -------
         Returns
            -------
            dict dict_Fbp
                fit-id (int) -> fit-id (int)
        """
        dict_Fbp = {}
        for base in self.design.scaffold:
            if base.across is not None:
                base_Fid = self.dict_DidFid[base.id]
                across_Fid = self.dict_DidFid[base.across.id]
                dict_Fbp[base_Fid] = across_Fid

        return dict_Fbp

    def link(self):
        """ invoke _link_scaffold, _link_staples, _link_bp to compute mapping
            of every base design-id to fit-id as well as the basepair mapping.
            basepairs are mapped from scaffold to staple, unique (invertable).
            updates linker attributes corresponding to the respective mapping
            and returns them.
        -------
         Returns
            -------
            dict self.dict_Fbp
                fit-id (int) -> fit-id (int)
            dict self.dict_DidFid
                design-id (int) -> fit-id (int)
            dict self.dict_DhpsDid
                helix-number (int), base-position (int), is_scaffold (bool)
                -> design-id (int)
            dict self.dict_color
                fit-segment-id (int) -> color (hex?)
        """
        dict_DidFid_sc, dict_DhpsDid_sc = self._link_scaffold()
        dict_DidFid_st, dict_DhpsDid_st, self.dict_color = self._link_staples()

        self.dict_DidFid = {**dict_DidFid_sc, **dict_DidFid_st}
        self.dict_DhpsDid = {**dict_DhpsDid_sc, **dict_DhpsDid_st}
        self.dict_Fbp = self._link_bp()

        return (self.dict_Fbp,
                self.dict_DidFid,
                self.dict_DhpsDid,
                self.dict_color,
                self.list_Dskips,
                )  # TODO: named tuple?

    def identify_crossover(self) -> Dict:
        """ for every base id that is involved in a crossover
            updates linker attribute of crossovers and returns it
        -------
         Returns
            -------
            dict self.dict_Fco
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
                Dhps = (h, p+i, is_scaf)
                while Dhps in self.list_Dskips:
                    i += np.sign(i)
                    Dhps = (h, p+i, is_scaf)
                leg_Did = self.dict_DhpsDid[Dhps]
                return self.dict_DidFid[leg_Did]

            def is_nextInStrand(b_Dhps: Tuple[int, int, bool],
                                a_Dhps: Tuple[int, int, bool]
                                ) -> bool:
                a_Fid = self.dict_DidFid[self.dict_DhpsDid[a_Dhps]]
                b_Fid = self.dict_DidFid[self.dict_DhpsDid[b_Dhps]]
                return True if abs(a_Fid - b_Fid) == 1 else False

            for co in iter(self.dict_Fco.values()):
                h, p, is_scaf = co["position"]
                for i in [-1, 1]:
                    n_Dhps = (h, p+i, is_scaf)
                    while n_Dhps in self.list_Dskips:
                        i += np.sign(i)
                        n_Dhps = (h, p+i, is_scaf)

                    n_Did = self.dict_DhpsDid.get(n_Dhps, None)
                    n_Fid = self.dict_DidFid.get(n_Did, None)
                    if n_Fid is None:
                        co["type"] = ("end", None, None)
                    elif is_nextInStrand(co["position"], n_Dhps):
                        continue
                    else:
                        is_double = (n_Fid in iter(self.dict_Fco.keys()))

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
            Dhps = (base.h, base.p+i, base.is_scaf)
            while Dhps in self.list_Dskips:
                i += np.sign(i)
                Dhps = (base.h, base.p+i, base.is_scaf)
            leg_Did = self.dict_DhpsDid[Dhps]
            return self.dict_DidFid[leg_Did]

        def get_next(base, direct: str):
            n = (base.up if direct == "up" else base.down)
            if n is None:
                return None
            else:
                n_position = (n.h, n.p, n.is_scaf)
                while n_position in self.list_Dskips:
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
            co_id = self.dict_DidFid[neighbor.id]
            leg_id = get_co_leg_id(base, direct)
            Dhps = (base.h, base.p, base.is_scaf)

            try:  # TODO -high cleanup co_index
                index = self.dict_Fco[co_id]["co_index"]
            except KeyError:
                run_id += 1
                index = run_id

            dict_data = {"co_index": index,
                         "co": co_id,
                         "leg": leg_id,
                         "is_scaffold": base.is_scaf,
                         "position": Dhps,
                         "base": base,
                         }
            return dict_data, run_id

        allbases = self.design.scaffold.copy()
        staples = [base for staple in self.design.staples for base in staple]
        allbases.extend(staples)

        self.dict_Fco = {}
        run_id = 0
        for base in allbases:
            base_Fid = self.dict_DidFid[base.id]

            for direct in ["up", "down"]:
                neighbor = get_next(base, direct)
                if is_co(base, neighbor, direct):
                    co_data, run_id = get_co(base, neighbor, direct, run_id)
                    self.dict_Fco[base_Fid] = co_data
                    break

        add_co_type()
        return self.dict_Fco

    def identify_nicks(self) -> Dict[int, int]:
        """ for every nick, provide id of base accross nick
        -------
         Returns
            -------
            dict self.dict_Fnicks
                fit-id (int) -> fit_id (int)
        """
        def is_nick(candidate, base) -> bool:
            is_onhelix = (candidate.h == base.h)
            is_neighbor = (abs(base.p - candidate.p) <= 2)  # skip = 2
            is_base = (candidate is base)
            b_Fid = self.dict_DidFid[base.id]
            c_Fid = self.dict_DidFid[candidate.id]
            is_ds_staple = (b_Fid in self.dict_Fbp.values() and c_Fid in self.dict_Fbp.values())
            if is_onhelix and is_neighbor and not is_base and not is_ds_staple:
                import ipdb
                ipdb.set_trace()
            return is_onhelix and is_neighbor and not is_base and is_ds_staple

        dict_Fnicks = {}
        staple_end_bases = []
        for staple in self.design.staples:
            staple_end_bases.extend((staple[0], staple[-1]))

        for base in staple_end_bases:
            for candidate in staple_end_bases:
                if is_nick(candidate, base):
                    base_Fid = self.dict_DidFid[base.id]
                    nick_Fid = self.dict_DidFid[candidate.id]
                    dict_Fnicks[base_Fid] = nick_Fid

        self.dict_Fnicks = dict_Fnicks
        return self.dict_Fnicks


class Fit(object):

    def __init__(self, path: str):
        self.path = path
        self.u = self._get_universe()
        self.scaffold, self.staples = self._split_strands()

    def _get_universe(self):
        # TODO: -mid- if no dcd, try invoke vmd animate write dcd from pdb
        top = self.path + ".psf"
        trj = self.path + ".dcd" #TODO with ignored
        u = mda.Universe(top, trj)
        return u

    def _split_strands(self):
        # TODO: -low- multiscaffold
        strands = self.u.segments
        scaffold = max(strands, key=attrgetter('residues.n_residues'))
        staples = [strand for strand in strands if len(
            strand.atoms) != len(scaffold.atoms)]
        return scaffold, staples


class Design(object):

    def __init__(self, path: str):
        self.path = path
        self.design = self._get_design()
        self.strands = self.design.strands
        self.scaffold = self._clean_scaffold(self.strands)
        self.excl = self.scaffold[0].strand
        self.staples = self._clean_staple(self.strands)
        self.dict_helixorder = self._create_helix_order()
        self.dict_stapleorder = self._create_staple_order()

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
        staples_clean = []
        for s in staples:
            tour = [b for b in s if not self._is_del(b)]
            staples_clean.append(tour)
        return staples_clean

    def _get_design(self):
        file_name = self.path + ".json"
        try:  # TODO with ignored
            seq_file = self.path + ".seq"
        except FileNotFoundError:
            seq_file = None
        seq_name = None

        converter = Converter()
        converter.read_cadnano_file(file_name, seq_file, seq_name)
        return converter.dna_structure

    def _create_helix_order(self) -> Dict[int, int]:
        """ helices are not listed in the json in same order as they are listed
            by helix-index. NOTE: enrgMD is helix reorder sensitive.
        -------
            Returns
            -------
            dict dict_helixorder
                (int) -> (int)
        """
        helices_dict = self.design.structure_helices_map
        dict_helixorder = {i: h.load_order for (i, h) in iter(helices_dict.items())}
        return dict_helixorder

    def _create_staple_order(self) -> Dict[int, int]:
        """ enrgMD and nanodesign number staples differently.
            enrgMD: first occurence of staple sorted by h, p
            nanodesign: 5' end of staple sorted by h, p
        -------
            Returns
            -------
            dict dict_stapleorder
                (int) -> (int)
        """
        list_Dhps = [(self.dict_helixorder[s[0].h], s[0].p)
                     for s in self.staples]
        list_Dhps_sorted = sorted(list_Dhps, key=lambda x: (x[0], x[1]))

        order_EMD = range(len(list_Dhps)) #TODO  use dict(enumerate()))
        order_ND = [list_Dhps.index(list_Dhps_sorted[i]) for i in order_EMD]
        dict_stapleorder = dict(zip(order_ND, order_EMD))

        return dict_stapleorder


def print_usage():  # TODO:cleanup
    print("""
    initializes MDAnalysis universe and creates dictionaries that link bse id
    in design and fit to design position, wc-partner, nicks, crossovers and
    staple color.

    usage: projectname designname

    return: creates the following pickles:
        ["dict_bp", "dict_idid", "dict_hpid", "dict_color", "dict_coid",
        "dict_nicks", "list_skips"]
        NOTE: (mda-universe cannot be pickled)
        """)


def proc_input():
    if len(sys.argv) < 2:
        print_usage()
        sys.exit(0)

    project = sys.argv[1]
    name = sys.argv[2]
    cwd = os.getcwd()
    in_put = cwd + "/" + project + "/" + name
    out_put = cwd + "/" + project + "/analysis/"
    try:
        os.mkdir(out_put)
    except FileExistsError:
        pass
    return in_put, out_put + name


def main():

    # process input
    in_put, out_put = proc_input()
    print("output to ", out_put)
    linker = Linker(in_put)
    universe = (in_put + ".psf", in_put + ".dcd")
    dict_bp, dict_idid, dict_hpid, dict_color, list_skips = linker.link()
    dict_coid = linker.identify_crossover()
    dict_nicks = linker.identify_nicks()
    for dict_name in DICTS:
        pickle.dump(eval(dict_name), open(
            out_put + "__" + dict_name + ".p", "wb"))


if __name__ == "__main__":
    main()
