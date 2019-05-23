#!/usr/bin/env python

import MDAnalysis as mda
import numpy as np
import sys
import os

from MDAnalysis.lib import mdamath
import MDAnalysis.analysis.distances
import matplotlib.pyplot as plt
import pickle
from scipy.optimize import curve_fit
import bDNA 

def initialize(path, devs):

    wc_pairs_dev = []
    wc_index_pairs_dev = []
    for i in range(len(devs)):
        _, files , wc_pairs, wc_index_pairs = pickle.load( open( path + "__wc_pairs-"+ str(devs[i])+ ".p", "rb" ) )
        wc_pairs_dev.append(wc_pairs)
        wc_index_pairs_dev.append(wc_index_pairs)
    u = mda.Universe(*files)
    u.trajectory[-1]

    strands = u.segments
    scaffold = strands[0]
    staples = strands[1:]
    
    return u, strands, scaffold, staples, devs, wc_pairs_dev, wc_index_pairs_dev

def fit_gauss(data, start, stop, nbins, guess):
    
    hist, bins,_= plt.hist(data, bins=nbins, range= (start, stop), density=True)
    bin_centres = (bins[:-1] + bins[1:])/2
    coeff, _ = curve_fit(bDNA.func, bin_centres, hist, p0=guess)
    
    return coeff


def c1(file, scaffold, staples):
    C1p_scaff = scaffold.atoms.select_atoms("name C1'")
    C1p_stap = staples.atoms.select_atoms("name C1'")


    dist_C1p_scaff = mda.analysis.distances.self_distance_array(C1p_scaff.positions)
    dist_C1p_stap = mda.analysis.distances.self_distance_array(C1p_stap.positions)
    dist_C1p_diff = mda.analysis.distances.distance_array(C1p_scaff.positions, C1p_stap.positions)
    
    
    nbins = 200
    
    hist_c1p_sc = dist_C1p_scaff[np.where(dist_C1p_scaff < 12.)]
    hist_c1p_st = dist_C1p_stap[np.where(dist_C1p_stap < 12.)]
    hist_c1p_diff = dist_C1p_diff[np.where(dist_C1p_diff < 13.)]
    guess_c1p = [0.3, 4., 1, 0.2, 10., 1]
    guess_wcmg = [0.05, 9., 2., 0.2, 11., 1,  0.2, 12., 1]


    plt.subplot(412)
    coeff_sc = fit_gauss(hist_c1p_sc, 0., 11., nbins, guess_c1p) 
    coeff_st = fit_gauss(hist_c1p_st, 0., 11., nbins, guess_c1p)
    coeff_diff = fit_gauss(hist_c1p_diff, 0, 13., nbins, guess_wcmg)
    
    c1p_mu1 = (coeff_sc[1] + coeff_st[1]) * 0.5
    c1p_mu2 = (abs(coeff_sc[4]) + abs(coeff_st[4])) * 0.5
    
    c1p_s1 = (coeff_sc[2] + coeff_st[2]) * 0.5
    c1p_s2 = (abs(coeff_sc[5]) + abs(coeff_st[5])) * 0.5
        
        
    print("C1'1: ,", c1p_mu1, ", ", c1p_s1, file= file)
    print("C1'2: ,", c1p_mu2, ", ", c1p_s2, file= file)
    print("wc    ,", coeff_diff[4], ", ", abs(coeff_diff[5]), file= file)
    print("minor ,", coeff_diff[7], ", ", abs(coeff_diff[8]), file= file)
    return

def rise(file, scaffold):
    C4 = scaffold.atoms.select_atoms("name C4").positions
    C2 = scaffold.atoms.select_atoms("name C2").positions
    N1 = scaffold.atoms.select_atoms("name N1").positions

    base = bDNA.plane_unitnormal_array(C4, C2, N1)

    #calculate normal distance
    dist_C4 = bDNA.distance_array_normal(C4, base, C4)
    C4_hist = dist_C4[np.where(dist_C4 < 5.)]
    
    #redo for all three if suspect error
    #dist_C2 = bDNA.distance_array_normal(N1, base, C4)
    #C2_hist = dist_C2[np.where(dist_C2 < 5.)]
    #dist_N1 = bDNA.distance_array_normal(N1, base, C4)
    #N1_hist = dist_N1[np.where(dist_N1 < 5.)]

    nbins = 200
    guess = [0.4, 4., 1.]

    coeff = fit_gauss(C4_hist, 0., 5., nbins, guess)
    print("rise  ,", coeff[1], ", ", abs(coeff[2]), file= file)
    return

def P(file, strand):
    
    P = strand.atoms.select_atoms("name P")

    nbins = 200
    guess_P = [0.1, 6., 1.]

    dist_P =mda.analysis.distances.self_distance_array(P.positions)
    P_hist = dist_P[np.where(dist_P < 9.)]

    coeff = fit_gauss(P_hist, 0., 9., nbins, guess_P)
    
    print("P     ,", coeff[1], ", ", abs(coeff[2]), file= file)

    return

def twist(file, u, scaffold, strands, wc_index_pairs_dev, devs):
    
    bp_angles_dev = []
    for j in range(len(devs)):
        bp_angles = []
        base_C1p = scaffold.residues[0].atoms.select_atoms("name C1'")
        vector = None
        v = wc_index_pairs_dev[j][base_C1p.residues[0].resindex]
        if v[0] is not None:
            wc_C1p = u.residues[v[0]].atoms.select_atoms("name C1'")
            vector = base_C1p.positions[0] - wc_C1p.positions[0]


        for i in range(1, len(strands[0].residues)):
            next_vector = None
            next_base_C1p = scaffold.residues[i].atoms.select_atoms("name C1'")
            next_v =  wc_index_pairs_dev[j][next_base_C1p.residues[0].resindex]
            if next_v[0] is not None:
                next_wc_C1p = u.residues[next_v[0]].atoms.select_atoms("name C1'")
                next_vector = next_base_C1p.positions[0] - next_wc_C1p.positions[0]

            if vector is not None and next_vector is not None:
                bp_angles.append(bDNA.angle(vector, next_vector))

            if next_v[0] is not None:
                wc_C1p = next_wc_C1p

            base_C1p = next_base_C1p        
            vector = next_vector
        bp_angles_dev.append(bp_angles)
    
        

        twist_hist = np.asarray(bp_angles, dtype=np.float64)*360./(2.*np.pi)
        guess_twist = [0.03, 21., 10., 0.04, 36., 10., 0.02, 60., 10.]
        coeff = fit_gauss(twist_hist, 0, 90., 200, guess_twist)

        count = 0
        for _, value in wc_index_pairs_dev[j].items(): 
            if value[0] is None: 
                  count += 1 
        bp_percentage = 1. - count/(len(scaffold.residues))

        print("tA"+str(devs[j])+",", coeff[1], ", ", abs(coeff[2]), file= file)
        print("tB"+str(devs[j])+",", coeff[4], ", ", abs(coeff[5]), file= file)
        print("tC"+str(devs[j])+",", coeff[7], ", ", abs(coeff[8]), file= file)
        print("bp%"+str(devs[j])+",", bp_percentage, ", ", devs[j], file= file)


    return

def print_usage():
    print("""
    initializes MDAnalysis universe and compoutes watson scric base pairs. they are returned as to dictionaries. this process is repeated for each Hbond-deviation criterion
    subsequently universe and dicts are stored into a pickle. each devia, tion criterion is stored in one pickle

    usage: projectname designname [dev = 0.1] [dev2] [dev3] ... 

    return: creates a pickle for each deviation: the pickle contains: (top, trj), wc_pairs, wc_index_pairs
            the tuple contains the absolute path of the files md-files (universe cannot be pickled), second and third are the two dictionaries
        """)


def proc_input(folder):
    if len(sys.argv) < 2:
        print_usage()

    # if isinstance(sys.argv[1], str):
    project = sys.argv[1]
    # if isinstance(sys.argv[2], str):
    name = sys.argv[2]
    cwd = os.getcwd()

    path = cwd + "/" + project + "/" + str(folder)+ "/" + name

    dev = []
    if len(sys.argv) == 3:
        dev.append(0.1)

    for av in sys.argv[3:]:
        dev.append(float(av))

    return  path, dev


def main():

    # perform for each step
    for i in reversed(range(-1, 8)):
        print("processing folder "+ str(i))
        # process input
        path, devs = proc_input(i)
        print("output to ", path)

        u, strands, scaffold, staples, devs, _, wc_index_pairs_dev = initialize(path, devs)

        if i == -1:
            u.trajectory[0]

        with open(path + "__data.csv", mode='w') as file:
            c1(file, scaffold, staples)
            P(file, scaffold)
            rise(file, scaffold)
            twist(file, u, scaffold, strands, wc_index_pairs_dev, devs)

   
if __name__ == "__main__":
    main()
