#!/usr/bin/env python3#
import sys
import os

import MDAnalysis as mda
import numpy as np
import pickle
import ipdb

from nanodesign.converters import Converter
from operator import attrgetter


class Linker(object):
    def __init__(self, path):
        self.path = path
        self.fit = Fit(self.path,)
        self.design = Design(self.path)
        self.d_Fbp, self.d_DidFid, self.d_DhpsDid = None, None, None
        self.d_Fco = None

    def _link_scaffold(self):
        """ My numpydoc description of a kind
            of very exhautive numpydoc format docstring.

            Parameters
            ----------
            first : array_like
                the 1st param name `first`
            second :
                the 2nd param
            third : {'value', 'other'}, optional
                the 3rd param, by default 'value'

            Returns
            -------
            string
                a value in a string

            Raises
            ------
            KeyError
                when a key error
            OtherError
                when an other error
            """
        """ idea: collect position in scaffold (0-x) by comparing to index in list of scaffold_design positions
            hp-id is design to design
            id-id is design to fit
        """

        idx_bases = []
        # get all possible positions for a scaffold base
        design_idx = [base.id for base in self.design.scaffold]
        design_hp = [(base.h, base.p, True) for base in self.design.scaffold]

        for base in self.design.scaffold:
            idx = design_idx.index(base.id)
            idx_bases.append(idx)

        res_ids = self.fit.scaffold.residues[idx_bases].resindices

        d_DidFid = dict(zip(design_idx, res_ids))
        d_hp = dict(zip(design_hp, design_idx))
        return d_DidFid, d_hp

    def _link_staples(self):
        """ loop over all staples, then perform the same procesdure as for scaffold
        """
        d_DidFid = {}
        d_hp = {}
        d_color = {}

        for i, staple in enumerate(self.design.staples):
            idx_segment = self.design.s_dict[i]

            # get all possible positions for a  base in this specific staple
            design_idx = [base.id for base in staple]
            design_hp = [(base.h, base.p, False) for base in staple]

            idx_bases = [design_idx.index(base.id) for base in staple]
            res_ids = [
                self.fit.staples[idx_segment].residues[j].resindex for j in idx_bases]

            # get color
            color = self.design.design.strands[staple[0].strand].icolor
            segidxforcolor = self.fit.staples[idx_segment].segindex
            d_color[segidxforcolor] = color

            idid_add = dict(zip(design_idx, res_ids))
            hp_add = dict(zip(design_hp, res_ids))
            d_DidFid = {**d_DidFid, **idid_add}
            d_hp = {**d_hp, **hp_add}

        return d_DidFid, d_hp, d_color

    def _link_bp(self, d_scaffold, d_staple):
        """ returns dict: key = id fit scaffold  value = id fit
            bp fit only
        """
        bp_u = {}
        for base in self.design.scaffold:
            if base.across is not None:
                bp_u[d_scaffold[base.id]] = d_staple[base.across.id]

        return bp_u

    def link(self):

        d_DidFid_sc, d_DhpsDid_sc = self._link_scaffold()
        d_DidFid_st, d_DhpsDid_st, self.d_color = self._link_staples()

        self.d_Fbp = self._link_bp(d_DidFid_sc, d_DidFid_st)
        self.d_DidFid = {**d_DidFid_sc, **d_DidFid_st}
        self.d_DhpsDid = {**d_DhpsDid_sc, **d_DhpsDid_st}

        return self.d_Fbp, self.d_DidFid, self.d_DhpsDid, self.d_color

    def identify_crossover(self):
        """ for every base id that is involved in a crossover, #ALL DESIGN
            list: co-partner id "co", "is_scaffold", "type" (single/double, id) 
        """
        def add_co_type(dict_co):
            for value in dict_co.values():
                h, p, is_scaf = value["position"]
                is_single = True
                for i in [-1, 1]:
                    try:
                        if self.d_DhpsDid[(h, p+i, is_scaf)] in dict_co.keys():
                            is_single = False
                            double = self.d_DhpsDid[(h, p+i, is_scaf)]
                    except KeyError:
                        pass #helix end

                value["type"] = ("single", None) if is_single else (
                    "double", double)
            return dict_co

        def get_co_leg_id(base, direct):
            # determine leg base (def: 2 bases away)
            l = base.down.p if direct == "up" else base.up.p
            i = (l - base.p) * 2.
            return self.d_DhpsDid[(base.h, base.p+i, base.is_scaf)]

        design_allbases = self.design.scaffold.copy()
        design_allbases.extend(
            [base for staple in self.design.staples for base in staple])

        dict_co = {}
        for design_base in design_allbases:
            for direct in ["up", "down"]:
                neighbor = design_base.up if direct == "up" else design_base.down

                if neighbor is not None:
                    if neighbor.h != design_base.h:  # crossover condition
                        leg_id = get_co_leg_id(design_base, direct)
                        co_id = self.d_DidFid[neighbor.id]
                        position = (design_base.h, design_base.p,
                                    design_base.is_scaf)
                        dict_co[self.d_DidFid[design_base.id]] = {
                            "co": co_id, "leg": leg_id, "is_scaffold": design_base.is_scaf, "position": position}

        dict_co = add_co_type(dict_co)

        self.d_Fco = dict_co
        return self.d_Fco

    def identify_nicks(self):
        # TODO: -low-
        raise NotImplementedError

    def _idx_incl(self, idx):
        """ include scaffold from idcount
        """
        return idx + 1 if idx >= self.design.excl else idx


class Fit(object):

    def __init__(self, path):
        self.path = path
        self.u = self._get_universe()
        self.scaffold, self.staples = self._split_strands()

    def _get_universe(self):
        top = self.path + ".psf"
        trj = self.path + ".dcd"
        u = mda.Universe(top, trj)
        return u

    def _split_strands(self):
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
        self.scaffold = self._clean_scaffold(
            self.strands)  # TODO: -low- multiscaffold
        self.excl = self.scaffold[0].strand
        self.staples = self._clean_staple(self.strands)
        self.h_dict = self._create_helix_order()
        self.s_dict = self._create_staple_order()

    def _clean_scaffold(self, strands):
        scaffold = [strand.tour for strand in strands if strand.is_scaffold][0]
        scaffold_clean = [
            base for base in scaffold if base.num_deletions == 0 and base.num_insertions == 0]
        return scaffold_clean

    def _clean_staple(self, strands):
        staples = [strand.tour for strand in strands if not strand.is_scaffold]
        staples_clean = []
        for stra in staples:
            staples_clean.append(
                [base for base in stra if base.num_deletions == 0 and base.num_insertions == 0])
        return staples_clean

    def _get_design(self):
        file_name = self.path + ".json"
        seq_file = self.path + ".seq"
        seq_name = None

        converter = Converter()
        converter.read_cadnano_file(file_name, seq_file, seq_name)
        return converter.dna_structure

    def _create_helix_order(self):
        """ helices are not prcessed in the same order as they are listed by idx
        """
        helices_dict = self.design.structure_helices_map
        h_dict = {i: h.load_order for (i, h) in helices_dict.items()}

        return h_dict

    def _create_staple_order(self):
        """ exchange staple id with load_id 
            map design-staple-order to universe-staple-order
            staple order is not the same for ergy-server and nd
            enrg/MDanalysis: sort left to right top to bottom; only counting 3' ends
            nanodesigns/auto: sort left to right top to bottom; counting every staple piece
        """
        list_hp = [(self.h_dict[s[0].h], s[0].p) for s in self.staples]
        list_hp_s = sorted(list_hp, key=lambda x: (x[0], x[1]))

        idx_list = [list_hp.index(list_hp_s[i]) for i in range(len(list_hp))]
        s_dict = dict(zip(idx_list, range(len(idx_list))))

        return s_dict


def print_usage():
    print("""
    initializes MDAnalysis universe and compoutes watson scric base pairs. they are returned as to dictionaries. this process is repeated for each Hbond-deviation criterion
    subsequently universe and dicts are stored into a pickle. each deviation criterion is stored in one pickle

    usage: projectname designname  

    return: creates a pickle for each deviation: the pickle contains: (top, trj), wc_pairs, wc_index_pairs
            the tuple contains the absolute path of the files md-files (universe cannot be pickled), second and third are the two dictionaries
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

    dict_bp, dict_idid, dict_hpid, dict_color = linker.link()
    dict_coid = linker.identify_crossover()
    for dict_name in ["dict_bp", "dict_idid", "dict__hpid", "dict_color", "dict_coid"]:
        pickle.dump(eval(dict_name), open(
            path + "__" + dict_name + ".p", "wb"))


if __name__ == "__main__":
    main()
