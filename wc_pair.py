#!/usr/bin/env python2
import __future__
import MDAnalysis as mda
import numpy as np
import sys
import os

from MDAnalysis.lib import mdamath
#import MDAnalysis.analysis.distances
#from math import sqrt
import pickle

import bpLinker #PYTHON3

C1P_BASEDIST = 10.7

WC_DICT = {"DC": "DG", "DG": "DC", "DT": "DA", "DA": "DT",
           "C": "G", "G": "C", "T": "A", "A": "T",
           "CYT": "GUA", "GUA": "CYT", "THY": "ADE", "ADE": "THY"}
WC_HBONDS = {"DA": ("N6", "N1"), "A": ("N6", "N1"), "ADE": ("N6", "N1"),
             "DT": ("O4", "N3"), "T": ("O4", "N3"), "THY": ("O4", "N3"),
             "DG": ("O6", "N1", "N2"), "G": ("O6", "N1", "N2"), "GUA": ("O6", "N1", "N2"),
             "DC": ("N4", "N3", "O2"), "C": ("N4", "N3", "O2"), "CYT": ("N4", "N3", "O2")}
(N6O4, N1N3_AT) = (2.95, 2.88)
(O6N4, N1N3_GC, N2O2) = (2.85, 2.85, 2.85)
WC_HBONDS_DIST = {"DA": (N6O4, N1N3_AT), "A": (N6O4, N1N3_AT), "ADE": (N6O4, N1N3_AT),
                  "DT": (N6O4, N1N3_AT), "T": (N6O4, N1N3_AT), "THY": (N6O4, N1N3_AT),
                  "DG": (O6N4, N1N3_GC, N2O2), "G": (O6N4, N1N3_GC, N2O2), "GUA": (O6N4, N1N3_GC, N2O2),
                  "DC": (O6N4, N1N3_GC, N2O2), "C": (O6N4, N1N3_GC, N2O2), "CYT": (O6N4, N1N3_GC, N2O2)}


class WcPairs(object):

    def __init__(self, universe, d_bp, d_idid, d_hpid):
        self.u = universe
        self.d_idid = d_idid
        self.d_bp = d_bp
        self.d_hpid =d_hpid
        self.wc_quality = None
        

    def eval_wc(self):
        """ return dict of wc-badness for each pair
        """
        self.wc_quality = {}
        for id1, id2 in self.d_bp.items():            
            res1 = self.u.residues[id1]
            res2 = self.u.residues[id2]

            self.wc_quality[id1], self.wc_quality[id2] = self._eval_wc_quality(res1, res2)
        return

    def _eval_wc_quality(self, res1, res2):
        quality = {}
        
        #C1p distance
        c1p1 = res1.atoms.select_atoms("name C1'")[0]
        c1p2 = res2.atoms.select_atoms("name C1'")[0]
        
        c1p_dist = mdamath.norm(c1p1.position - c1p2.position)
        bond = "C1'C1'"
        c1p_dev = abs(c1p_dist - C1P_BASEDIST)/C1P_BASEDIST
        quality[bond] = c1p_dev
        
        #H-bond distances
        atoms1 = res1.atoms.select_atoms(
            "name " + ' or name '.join(map(str, WC_HBONDS[res1.resname])))
        atoms2 = res2.atoms.select_atoms(
            "name " + ' or name '.join(map(str, WC_HBONDS[res2.resname])))
        
        
        for idx, a1 in enumerate(atoms1):
            hbond_dist = mdamath.norm(atoms2[idx].position - a1.position)
            hbond_should = WC_HBONDS_DIST[res1.resname][idx]
            hbond_dev = abs(hbond_dist - hbond_should)/hbond_should
            bond = WC_HBONDS[res1.resname][idx] + WC_HBONDS[res2.resname][idx]
            quality[bond] = hbond_dev
        
        return quality, quality


def print_usage():
    print("""
    initializes MDAnalysis universe and compoutes watson scric base pairs. they are returned as to dictionaries. this process is repeated for each Hbond-deviation criterion
    subsequently universe and dicts are stored into a pickle. each deviation criterion is stored in one pickle

    usage: projectname designname frames [dev = 0.1] [dev2] [dev3] ... 

    return: creates a pickle for each deviation: the pickle contains: (top, trj), wc_pairs, wc_index_pairs
            the tuple contains the absolute path of the files md-files (universe cannot be pickled), second and third are the two dictionaries
        """)


def proc_input():
    
    
    if len(sys.argv) < 3:
        print_usage()

    # if isinstance(sys.argv[1], str):
    project = sys.argv[1]
    # if isinstance(sys.argv[2], str):
    name = sys.argv[2]
    cwd = os.getcwd()
    base = cwd + "/" + project + "/"

    top = base + name + ".psf"
    trj = base + name + ".dcd"
    output = base + "analysis/"
    try:     
        os.mkdir(output)
    except  OSError: #FileExistsError: #PYTHON3
        pass

    n_frames = sys.argv[3]

    dev = []
    if len(sys.argv) == 4:
        dev.append(0.1)

    for av in sys.argv[4:]:
        dev.append(float(av))

    return top, trj, output, name, dev, n_frames


def main():
        top, trj, output, name, deviations, n_frames = proc_input()
        print("read ", top, trj, deviations, n_frames)
        print("output to ", output)
 


        # initialize universe and select final frame
        u = mda.Universe(top, trj)

        frames_step = int( len(u.trajectory) /  n_frames)
        frames = list(range(len(u.trajectory),0,-frames_step))

        linker = bpLinker.Linker( output + name )
        d_bp, d_idid, d_hpid =  linker.link()
            
        pickle.dump(d_bp, open( output + name + "__bp-dict.p", "wb"))
        pickle.dump(d_idid, open( output + name + "__idid-dict.p", "wb"))
        pickle.dump(d_hpid, open( output + name + "__hpid-dict.p", "wb"))

        for ts in [u.trajectory[i] for i in frames]:
            ts = ts +1
            #evaluate wc_parispair correctness

            #perform analyis

            


            #for dev in deviations:
             #   print("performing wc- analysis for dev = ", dev)
              #  wc_pairs, wc_index_pairs = get_wc_dict(u, Hbond_dev=dev)
               # pickle.dump(( deviations, (top, trj), wc_pairs, wc_index_pairs), open(
                #    str(output) + "__wc_pairs-" + str(dev) + ".p", "wb"))
    


if __name__ == "__main__":
    main()
