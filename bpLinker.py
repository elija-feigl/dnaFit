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
        self.d_bp, self.d_idid, self.d_hpid = None, None, None
        self.d_co_id = None

    def _link_scaffold(self):
        """idea: collect position in scaffold (0-x) by comparing to index in list of scaffold_design positions"""

        idx_bases = []
        # get all possible positions for a scaffold base
        design_idx = [base.id for base in self.design.scaffold]
        design_hp = [(base.h, base.p, 0) for base in self.design.scaffold]

        for base in self.design.scaffold:
            idx = design_idx.index(base.id)
            idx_bases.append(idx)

        res_ids = self.fit.scaffold.residues[idx_bases].resindices

        d_idid = dict(zip(design_idx, res_ids))
        d_hp = dict(zip(design_hp, res_ids))
        return d_idid, d_hp
  

    def _link_staples(self):
        """ loop over all staples, then perform the same procesdure as for scaffold
        """
        d_idid = {}
        d_hp = {}
        d_color = {}

        for i, staple in enumerate(self.design.staples):
            idx_segment = self.design.s_dict[i]
            #print("seg", idx_segment, "staple runs", (staple[0].h,staple[0].p), ":", (staple[-1].h,staple[-1].p) )

            idx_bases = []
            # get all possible positions for a  base in this specific staple
            design_idx = [base.id for base in staple]
            design_hp = [(base.h, base.p, 1)for base in staple]
            for base in staple:
                idx = design_idx.index(base.id)
                idx_bases.append(idx)

            # TODO: cleanup after final testing
            res_ids = []
            for j in idx_bases:
                try:
                    res_ids.append(
                        self.fit.staples[idx_segment].residues[j].resindex)

                except:
                    print("something went wrong - excluded base ",
                          j, "on helix ", staple[j].h, )
                    res_ids.append(-1)
                    pass
            # get color
            color = self.design.design.strands[staple[0].strand].icolor
            segidxforcolor = self.fit.staples[idx_segment].segindex
            d_color[segidxforcolor] = color
            

            idid_add = dict(zip(design_idx, res_ids))
            hp_add = dict(zip(design_hp, res_ids))
            d_idid = {**d_idid, **idid_add}
            d_hp = {**d_hp, **hp_add}

        return d_idid, d_hp, d_color

    def _link_bp(self, d_scaffold, d_staple):
        """ returns id id in fit?
        """
        bp_u = {}

        for base in self.design.scaffold:
            if base.across is not None:
                bp_u[d_scaffold[base.id]] = d_staple[base.across.id]

        return bp_u

    def link(self):

        d_idid_sc, d_hpid_sc = self._link_scaffold()
        d_idid_st, d_hpid_st, self.d_color = self._link_staples()

        self.d_bp = self._link_bp(d_idid_sc, d_idid_st)
        self.d_idid = {**d_idid_sc, **d_idid_st}
        self.d_hpid = {**d_hpid_sc, **d_hpid_st}

        return self.d_bp, self.d_idid, self.d_hpid, self.d_color

    def identify_crossover(self):
        """ both staple and scaffold crossover
        """
        #TODO: -high- ceck co_set, add pairings, add type (single/double ,staple/scaffold)


        design_allbases = self.design.scaffold.copy()
        design_allbases.extend(
            [base for staple in self.design.staples for base in staple])

        set_co_designid = {}
        for design_base in design_allbases:
            for direct in ["up", "down"]:
                try: 
                    neighbor = design_base.up if direct == "up" else design_base.down
                    if neighbor.h != design_base.h:
                        set_co_designid[design_base.id] = {{"co": neighbor.id}, {"is_scaf": design_base.is_scaf}}
                except AttributeError:
                    pass
            
        #TODO: -mid- append single/double; loop over all and check if .p +-1  (is_scaff) -> id in list --> double else single
        self.d_co_id = set_co_designid
        return self.d_co_id, 

    def identify_nicks(self):
        """ staple nicks
       

        design_allbases = self.design.scaffold.copy()
        design_allbases.extend(
            [base for staple in self.design.staples for base in staple])

        set_co_designid = set()
        for design_base in design_allbases:
            if design_base.id in set_co_designid:
                continue

            try:  # TODO: -low- DRY
                up = design_base.up
                if up.h != design_base.h:
                    set_co_designid.update([design_base.id, up.id])
            except AttributeError:
                pass
            try:
                down = design_base.down
                if down.h != design_base.h:
                    set_co_designid.update([design_base.id, down.id])
            except AttributeError:
                pass

        self.s_co_id = set_co_designid
         """
        return self.s_co_id

    def _idx_incl(self, idx):
        """ include scaffold from idcount
        """
        if idx >= self.design.excl:
            st_idx = idx + 1
        else:
            st_idx = idx

        return st_idx


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
        staples = [strand for strand in strands if len(strand.atoms) != len(scaffold.atoms)]
        return scaffold, staples


class Design(object):

    def __init__(self, path):
        self.path = path
        self.design = self._get_design()
        self.strands = self.design.strands
        self.scaffold = self._clean_scaffold(
            self.strands)  # TODO: multiscaffold
        # needed to correc tstaple numbering
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

    d_bp, d_idid, d_hpid, d_color = linker.link()
    d_coid = linker.identify_crossover()
    pickle.dump(d_bp, open(path + "__bp-dict.p", "wb"))
    pickle.dump(d_idid, open(path + "__idid-dict.p", "wb"))
    pickle.dump(d_hpid, open(path + "__hpid-dict.p", "wb"))
    pickle.dump(d_color, open(path + "__color-dict.p", "wb"))
    pickle.dump(d_coid, open(path + "__coid-dict.p", "wb"))


if __name__ == "__main__":
    main()
