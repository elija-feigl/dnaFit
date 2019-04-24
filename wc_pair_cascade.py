#!/usr/bin/env python

import MDAnalysis as mda
import numpy as np
import sys
import os

from MDAnalysis.lib import mdamath
import MDAnalysis.analysis.distances
from math import sqrt
import pickle

from constants import *
from wc_pair import get_wc, get_wc_dict



def print_usage():
    print("""
    initializes MDAnalysis universe and compoutes watson scric base pairs. they are returned as to dictionaries. this process is repeated for each Hbond-deviation criterion
    subsequently universe and dicts are stored into a pickle. each deviation criterion is stored in one pickle

    usage: projectname designname [dev = 0.1] [dev2] [dev3] ... 

    return: creates a pickle for each deviation: the pickle contains: (top, trj), wc_pairs, wc_index_pairs
            the tuple contains the absolute path of the files md-files (universe cannot be pickled), second and third are the two dictionaries
        """)


def proc_input(i):
    

    
    if len(sys.argv) < 2:
        print_usage()

    # if isinstance(sys.argv[1], str):
    project = sys.argv[1]
    # if isinstance(sys.argv[2], str):
    name = sys.argv[2]
    cwd = os.getcwd()


    if i == -1:
        try:
            os.mkdir(cwd + "/" + project + "/" + str(i))
        except FileExistsError:
            pass
        j = 0
    else:
        j = i

    top = cwd + "/" + project + "/" + name + ".psf"
    trj = cwd + "/" + project + "/" + str(j) + "/" + name + ".dcd"
    output = cwd + "/" + project + "/" + str(i)+ "/" + name

    dev = []
    if len(sys.argv) == 3:
        dev.append(0.1)

    for av in sys.argv[3:]:
        dev.append(float(av))

    return top, trj, output, dev


def main():


    for i in reversed(range(-1, 8)):
        print("processing folder "+ str(i))
        # process input
        top, trj, output, deviations = proc_input(i)
        print("read ", top, trj, deviations)
        print("output to ", output)


        # initialize universe and select final frame
        u = mda.Universe(top, trj)

        if i == -1:
            u.trajectory[0]
        else:    
            u.trajectory[-1]

        for dev in deviations:
            print("performing wc- analysis for dev = ", dev)
            wc_pairs, wc_index_pairs = get_wc_dict(u, Hbond_dev=dev)
            pickle.dump(( deviations, (top, trj), wc_pairs, wc_index_pairs), open(
                str(output) + "__wc_pairs-" + str(dev) + ".p", "wb"))
    


if __name__ == "__main__":
    main()
