#!/usr/bin/env python3#
import sys
import os

import MDAnalysis as mda
import numpy as np
import pickle

from nanodesign.converters import Converter
from operator import attrgetter

DICTS = ["dict_bp", "dict_idid", "dict_hpid", "dict_color",
         "dict_coid", "dict_nicks", "list_skips"]
# TODO: -low move


class Linker(object):
    """ Linker class
    """

    def __init__(self, path):
        self.path = path
        self.fit = Fit(self.path,)
        self.design = Design(self.path)
        self.d_Fbp = None
        self.d_DidFid = None
        self.d_DhpsDid = None
        self.d_Fco = None
        self.l_Dskips = self._get_skips()
        self.d_Fnicks = None
        self.d_FidSeq = self._get_sequence()
        # TODO: d_DhpsFid

    def _get_sequence(self):
        d_FidSeq = {}
        for res in self.fit.u.residues:
            d_FidSeq[res.resindex] = res.resname[0]
        return d_FidSeq

    def _is_del(self, base):
        return base.num_deletions != 0

    def _get_skips(self):
        """ collect position of skips in design
        -------
            Returns
            -------
            list l_Dskips
                (base.h, base.p, base.is_scaf) of all skips
        """
        design_allbases = [
            base for strand in self.design.strands for base in strand.tour]
        l_Dskips = []
        for base in design_allbases:
            if self._is_del(base):
                hp_s = (base.h, base.p, base.is_scaf)
                l_Dskips.append(hp_s)

        return l_Dskips

    def _link_scaffold(self):
        """ collect position in scaffold (0-x) by comparing to index in list
            of scaffold_design positions
        -------
            Returns
            -------
            dict d_DidFid
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

        d_DidFid = dict(zip(D_ids, F_id_global))
        d_DhpsDid = dict(zip(D_hp, D_ids))
        return d_DidFid, d_DhpsDid

    def _link_staples(self):
        """same procedure as scaffold for each
        -------
         Returns
            -------
            dict d_DidFid
                design-id (int) -> fit-id (int)
            dict d_DhpsDid
                helix-number (int), base-position (int), is_scaffold (bool)
                -> design-id (int)
            dict d_color
                fit-segment-id (int) -> color (hex?)
        """
        def get_resid(segindex, resindex_local):
            return self.fit.staples[segindex].residues[resindex_local].resindex

        d_DidFid = {}
        d_DhpsDid = {}
        d_color = {}

        for i, staple in enumerate(self.design.staples):
            seg_id = self.design.s_dict[i]

            D_ids = [base.id for base in staple]
            D_hp = [(base.h, base.p, False) for base in staple]

            F_id_local = [D_ids.index(base.id)for base in staple]
            F_id_global = [get_resid(seg_id, resid) for resid in F_id_local]

            color = self.design.design.strands[staple[0].strand].icolor
            segidxforcolor = self.fit.staples[seg_id].segindex
            d_color[segidxforcolor] = color

            d_DidFid_add = dict(zip(D_ids, F_id_global))
            d_DhpsDid_add = dict(zip(D_hp, D_ids))
            d_DidFid = {**d_DidFid, **d_DidFid_add}
            d_DhpsDid = {**d_DhpsDid, **d_DhpsDid_add}

        return d_DidFid, d_DhpsDid, d_color

    def _link_bp(self):
        """ link basepairs by mapping indices according to json (cadnano).
            basepairs are mapped from scaffold to staple, unique (invertable).
        -------
         Returns
            -------
            dict d_Fbp
                fit-id (int) -> fit-id (int)
        """
        d_Fbp = {}
        for base in self.design.scaffold:
            if base.across is not None:
                base_Fid = self.d_DidFid[base.id]
                across_Fid = self.d_DidFid[base.across.id]
                d_Fbp[base_Fid] = across_Fid

        return d_Fbp

    def link(self):
        """ invoke _link_scaffold, _link_staples, _link_bp to compute mapping
            of every base design-id to fit-id as well as the basepair mapping.
            basepairs are mapped from scaffold to staple, unique (invertable).
            updates linker attributes corresponding to the respective mapping
            and returns them.
        -------
         Returns
            -------
            dict self.d_Fbp
                fit-id (int) -> fit-id (int)
            dict self.d_DidFid
                design-id (int) -> fit-id (int)
            dict self.d_DhpsDid
                helix-number (int), base-position (int), is_scaffold (bool)
                -> design-id (int)
            dict self.d_color
                fit-segment-id (int) -> color (hex?)
        """
        d_DidFid_sc, d_DhpsDid_sc = self._link_scaffold()
        d_DidFid_st, d_DhpsDid_st, self.d_color = self._link_staples()

        self.d_DidFid = {**d_DidFid_sc, **d_DidFid_st}
        self.d_DhpsDid = {**d_DhpsDid_sc, **d_DhpsDid_st}
        self.d_Fbp = self._link_bp()

        return (self.d_Fbp,
                self.d_DidFid,
                self.d_DhpsDid,
                self.d_color,
                self.l_Dskips,
                )

    def identify_crossover(self):
        """ for every base id that is involved in a crossover
            updates linker attribute of crossovers and returns it
        -------
         Returns
            -------
            dict self.d_Fco
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
        def add_co_type(dict_co):  # TODO: -midclean-
            for value in dict_co.values():
                h, p, is_scaf = value["position"]
                is_double, is_single, is_end = False, False, False
                n_Fid = None
                for i in [-1, 1]:
                    try:
                        if (h, p+i, is_scaf) in self.l_Dskips:
                            i *= 2.
                        n_Fid = self.d_DidFid[self.d_DhpsDid[(
                            h, p+i, is_scaf)]]
                        try:
                            leg = self.d_DidFid[self.d_DhpsDid[(
                                h, p+2*i, is_scaf)]]
                        except KeyError:
                            leg = self.d_DidFid[self.d_DhpsDid[(
                                h, p+3*i, is_scaf)]]
                    except KeyError:
                        is_end = True

                    if n_Fid in dict_co.keys():
                        double = n_Fid
                        is_double = True
                    elif not is_end:
                        neighbors = []
                        for direct in ["up", "down"]:
                            neigh = get_next(value["base"], direct)
                            if neigh is not None:
                                while self._is_del(neigh):
                                    n_new = get_next(neigh, direct)
                                    if n_new is None:
                                        break
                                    neigh = n_new
                            neighbors.append(self.d_DidFid[neigh.id])
                        if n_Fid not in neighbors:
                            single = n_Fid
                            single_leg = leg
                            is_single = True

                if is_end:
                    value["type"] = ("end", None)
                elif is_single:
                    value["type"] = ("single", single, single_leg)
                elif is_double:
                    value["type"] = ("double", double)
                else:
                    exit(0)

                value.pop("base")
            return dict_co

        def get_co_leg_id(base, direct):
            """determine leg base (def: 2 bases away)"""
            l = base.down.p if direct == "up" else base.up.p
            i = (l - base.p) * 2.
            hp_s = (base.h, base.p+i, base.is_scaf)
            while hp_s in self.l_Dskips:
                i = i + np.sign(i)
                hp_s = (base.h, base.p+i, base.is_scaf)
            leg_Did = self.d_DhpsDid[hp_s]
            return self.d_DidFid[leg_Did]

        def get_next(base, direct):
            return (base.up if direct == "up" else base.down)

        design_allbases = self.design.scaffold.copy()
        staples = [base for staple in self.design.staples for base in staple]
        design_allbases.extend(staples)

        dict_co = {}
        co_running_index = 0
        for design_base in design_allbases:
            base_Fid = self.d_DidFid[design_base.id]

            for direct in ["up", "down"]:
                neighbor = get_next(design_base, direct)
                if neighbor is not None:
                    while self._is_del(neighbor):
                        n_new = get_next(neighbor, direct)
                        if n_new is None:
                            break
                        neighbor = n_new
                    if neighbor.h != design_base.h:  # crossover condition
                        leg_id = get_co_leg_id(design_base, direct)
                        try:
                            co_id = self.d_DidFid[neighbor.id]
                        except KeyError:
                            exit(0)
                        position = (design_base.h, design_base.p,
                                    design_base.is_scaf)

                        try:
                            co_index = dict_co[co_id]["co_index"]
                        except KeyError:
                            co_running_index += 1
                            co_index = co_running_index

                        dict_co[base_Fid] = {
                            "co_index": co_index,
                            "co": co_id,
                            "leg": leg_id,
                            "is_scaffold": design_base.is_scaf,
                            "position": position,
                            "base": design_base,
                        }

        dict_co = add_co_type(dict_co)

        self.d_Fco = dict_co
        return self.d_Fco

    def identify_nicks(self):
        """ for every nick, provide id of base accross nick
        -------
         Returns
            -------
            dict self.d_Fnicks
                fit-id (int) -> fit_id (int)
        """
        def is_nick(candidate, base):
            is_onhelix = (candidate.h == base.h)
            is_neighbor = (abs(base.p - candidate.p) <= 2)  # skip = 2
            is_base = (candidate is not base)
            return is_onhelix and is_neighbor and not is_base

        d_Fnicks = {}
        staple_end_bases = []
        for staple in self.design.staples:
            staple_end_bases.extend((staple[0], staple[-1]))

        for base in staple_end_bases:
            for candidate in staple_end_bases:
                if is_nick(candidate, base):
                    base_Fid = self.d_DidFid[base.id]
                    nick_Fid = self.d_DidFid[candidate.id]
                    d_Fnicks[base_Fid] = nick_Fid

        self.d_Fnicks = d_Fnicks
        return self.d_Fnicks


class Fit(object):

    def __init__(self, path):
        self.path = path
        self.u = self._get_universe()
        self.scaffold, self.staples = self._split_strands()

    def _get_universe(self):
        # TODO: -mid- if no dcd, try invoke vmd animate write dcd from pdb
        top = self.path + ".psf"
        trj = self.path + ".dcd"
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

    def __init__(self, path):
        self.path = path
        self.design = self._get_design()
        self.strands = self.design.strands
        self.scaffold = self._clean_scaffold(self.strands)
        self.excl = self.scaffold[0].strand
        self.staples = self._clean_staple(self.strands)
        self.h_dict = self._create_helix_order()
        self.s_dict = self._create_staple_order()

    def _is_del(self, base):
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
        try:
            seq_file = self.path + ".seq"
        except FileNotFoundError:
            seq_file = None
        seq_name = None

        converter = Converter()
        converter.read_cadnano_file(file_name, seq_file, seq_name)
        return converter.dna_structure

    def _create_helix_order(self):
        """ helices are not listed in the json in same order as they are listed
            by helix-index. NOTE: enrgMD is helix reorder sensitive.
        -------
            Returns
            -------
            dict h_dict
                (int) -> (int)
        """
        helices_dict = self.design.structure_helices_map
        h_dict = {i: h.load_order for (i, h) in helices_dict.items()}
        return h_dict

    def _create_staple_order(self):
        """ enrgMD and nanodesign number staples differently.
            enrgMD: first occurence of staple sorted by h, p
            nanodesign: 5' end of staple sorted by h, p
        -------
            Returns
            -------
            dict s_dict
                (int) -> (int)
        """
        list_hp = [(self.h_dict[s[0].h], s[0].p) for s in self.staples]
        list_hp_s = sorted(list_hp, key=lambda x: (x[0], x[1]))

        idx_list = [list_hp.index(list_hp_s[i]) for i in range(len(list_hp))]
        s_dict = dict(zip(idx_list, range(len(idx_list))))

        return s_dict


def print_usage():
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

    output = cwd + "/" + project + "/" + name

    return output


def main():

    # process input
    path = proc_input()
    print("output to ", path)
    linker = Linker(path)

    dict_bp, dict_idid, dict_hpid, dict_color, list_skips = linker.link()
    dict_coid = linker.identify_crossover()
    dict_nicks = linker.identify_nicks()
    for dict_name in DICTS:
        pickle.dump(eval(dict_name), open(
            path + "__" + dict_name + ".p", "wb"))


if __name__ == "__main__":
    main()
