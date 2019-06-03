#!/usr/bin/env python3
import MDAnalysis as mda
import numpy as np
import sys
import os

import ipdb
from MDAnalysis.lib import mdamath

import pickle

import bpLinker

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
DH_ATOMS = {"alpha":("O3' -","P","O5'","C5'"), "beta":("P","O5'","C5'","C4'"),
                        "gamma":("O5'","C5'","C4'","C3'"), "delta":("C5'","C4'","C3'","O3'"),
                        "epsilon":("C4'","C3'","O3'","P +"), "zeta":("C3'","O3'","P +","O5' +"),
                        "xi":{"pyr":("C4'","C1'","N1","C2"), "pur":("C4'","C1'","N9","C4")} }
BB_ATOMS = ["C1'","O3'","C3'","C4'", "O5'","C5'","P"]
PYR_ATOMS = ["N1","C2"]
PUR_ATOMS = ["N9","C4"]



class BDna(object):

    def __init__(self, universe, d_bp, d_idid, d_hpid):
        self.u = universe
        self.d_idid = d_idid
        self.d_bp = d_bp
        self.d_hpid =d_hpid

        self.d_crossoverid = None #derive from dicts? #TODO: -mid-

        self.wc_quality = None   #dict: references residue index -> bond quality 
        self.wc_geometry = None   
        self.dh_quality = None  
        self.distances = None    #only scaffold bases have complement

    def _get_next_wc(self, resindex, resindex_wc):
        """ get next residue and its complemt wc pair
            check if next exists.
            check has bp (ony scaffold residues apply here)
            TODO: -low- use resindex_wc to check 
        """
        n_resindex = resindex + 1
        
        try:
            self.u.residues[n_resindex]
            try:
                n_resindex_wc = self.d_bp[n_resindex]
            except KeyError:
                n_resindex_wc = None
        except KeyError:
            n_resindex = None
            n_resindex_wc = None

        return n_resindex, n_resindex_wc 
    
    def eval_wc(self):
        """  dict of wc-quality for each residue (key= resindex)
        """
        self.wc_quality = {}
        self.wc_geometry = {}
        for resindex, resindex_wc in self.d_bp.items():      

            res = self.u.residues[resindex]
            res_wc = self.u.residues[resindex_wc]

            self.wc_quality[resindex], self.wc_quality[resindex_wc] = self._get_wc_quality(res, res_wc)

            n_resindex, n_resindex_wc = self._get_next_wc(resindex, resindex_wc)

            if n_resindex is not None and n_resindex_wc is not None:
                n_res = self.u.residues[n_resindex]
                n_res_wc = self.u.residues[n_resindex_wc]
                self.wc_geometry[resindex], self.wc_geometry[resindex_wc] = self._get_wc_geometry(res, res_wc, n_res, n_res_wc)
 
        return

    def _get_wc_quality(self, res1, res2):
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
            quality[bond] = hbond_dev if hbond_dist > hbond_should else -hbond_dev
        
        return quality, quality

    def _get_wc_geometry(self, res, res_wc, n_res, n_res_wc):
        """  dict of twist for each wc-pair (key= resid-scaffold)
        """
        
        def _get_twist(basepair, n_basepair):
            twist = []
            for direct in ["dir-anker", "dir-center"]:
                v1 = basepair[direct]
                v2 = n_basepair[direct]
                twist.append(np.rad2deg(np.arccos(np.inner(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))))
             
            return {"anker":twist[0], "center":twist[1]}

        def _get_rise(basepair, n_basepair):
            rise = []
            n0 = basepair["n0"]
            for direct in ["center-anker", "center"]:
                P1 = basepair[direct]
                P2 = n_basepair[direct]
                rise.append(np.abs(np.inner((P2 - P1), n0)))

            return {"anker":rise[0], "center":rise[1]}

        bases = []
        for r in [res, res_wc, n_res, n_res_wc]:
            bases.append(self._get_base_plane(r))

        basepairs = []
        for i in [0,2]:
            basepairs.append({"center-anker": ((bases[i]["anker"] + bases[i+1]["anker"]) * 0.5),
                    "dir-anker": (bases[i]["anker"] - bases[i+1]["anker"]) ,
                    "center": ((bases[i]["center"] + bases[i+1]["center"]) * 0.5), 
                    "dir-center": (bases[i]["anker"] - bases[i+1]["anker"]) ,
                    "n0": ((bases[i]["n0"] + bases[i+1]["n0"]) * 0.5) })

        geometry = {"twist": _get_twist(*basepairs), "rise": _get_rise(*basepairs)}
        return geometry, geometry

    def _get_base_plane(self, res):

        atom = []
        for atom_name in ["C2", "C4", "C6", "C1'"]:
            A = res.atoms.select_atoms("name " + atom_name)[0]
            atom.append(A.position)

        nhat = np.cross((atom[1]-atom[0]), (atom[2]-atom[1]))
        n0 = nhat/np.linalg.norm(nhat)
    
        diazine_center = sum(atom[:-1]) / 3.

        plane = {"n0": n0, "anker": atom[-1], "center": diazine_center}
        
        return plane

    def eval_dh(self):
        """  dict of dihedrals for each residue (key= resid)
        """
        self.dh_quality = {}
        for res in self.u.residues: 
            self.dh_quality[res.resindex] = self._get_dihedrals(res)
        return

    def _get_dihedrals(self, res):

        atoms, logic = self._get_residue_BB(res)
        dh_valid = self._get_valid_DH(*logic)
        if res.resname in ["ADE", "GUA"]:
            pyr = False
        else:
            pyr = True

        dh = self._get_dh_for_res(atoms, pyr, dh_valid) 
        return dh

    def _get_residue_BB(self, res):
        iniSeg, terSeg, ter5 = False, False, False
        
        atoms = {}
        try:
            P = res.atoms.select_atoms("name " + "P")[0]
            atoms["P"] = P.position

        except (KeyError, IndexError):
            ter5 = True 

        for xx in BB_ATOMS[:-1]:
            atoms[xx] = res.atoms.select_atoms("name " + xx)[0].position
        if res.resname in ["ADE", "GUA"]:
            YY = PUR_ATOMS
        else:
            YY = PYR_ATOMS
        for yy in YY:
            atoms[yy] = res.atoms.select_atoms("name " + yy)[0].position

        try:
            n_res = self.u.residues[res.resindex +1 ]
            if res.segindex == n_res.segindex:
                n_P = n_res.atoms.select_atoms("name " + "P")[0]
                n_O5p = n_res.atoms.select_atoms("name " + "O5'")[0]
                atoms["P +"] = n_P.position
                atoms["O5' +"] = n_O5p.position
            else:
                terSeg = True
        except (KeyError, IndexError):
            terSeg = True

        try:
            p_res = self.u.residues[res.resindex - 1]
            if res.segindex == p_res.segindex:
                p_O3p = p_res.atoms.select_atoms("name " + "O3'")[0]
                atoms["O3' -"] = p_O3p.position
            else:
                iniSeg = True
        except (KeyError, IndexError):
            iniSeg = True
                        
        return atoms, (ter5, terSeg, iniSeg)

    def _get_valid_DH(self, ter5, terSeg, iniSeg):   
        dh_valid = ["gamma","delta", "xi"]
        if terSeg is False:
            dh_valid.extend(["epsilon","zeta"])
        if iniSeg is False:
            dh_valid.append("alpha")
        if ter5 is False:
            dh_valid.append("beta")
        return dh_valid

    def _get_dh_for_res(self, atoms, pyr, dh_valid):
        dh = {}
        for dh_name in ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "xi"]:
            if dh_name in dh_valid:
                angle = self._get_angle(atoms, pyr, dh_name)
            else:
                angle = None
            dh[dh_name] = angle
        return dh

    def _get_angle(self, atoms, pyr, dh_name):

        def _dh_angle(p1, p2, p3, p4, as_rad= False):
      
            v1 = p2 - p1
            v2 = p3 - p2
            v3 = p4 - p3
            
            n1 = np.cross(v1,v2)
            n1 /= np.linalg.norm(n1)
            n2 = np.cross(v2,v3)
            n2 /= np.linalg.norm(n2)
            m1 = np.cross(n1, v2 / np.linalg.norm(v2))
        
            x = np.dot(n1, n2)
            y = np.dot(m1, n2)

            angle = - np.arctan2(y, x)
            
            return angle if as_rad else np.rad2deg(angle)
        
        p = []
        
        if dh_name == "xi":
            if pyr:
                for i in DH_ATOMS[dh_name]["pyr"]:
                    p.append(atoms[i])
            else:
                for i in DH_ATOMS[dh_name]["pur"]:
                    p.append(atoms[i])
        else:
            for i in DH_ATOMS[dh_name]:
                p.append(atoms[i])

        angle = _dh_angle(*p) 
        
        return angle
    

    def eval_distances(self):
        """  dict of twist for each wc-pair (key= resid-scaffold)
        """
        self.distances = {}
        dist_C1p = self._get_A_distance("C1'")
        dist_P = self._get_A_distance("P")
        for resindex in range(self.u.residues.n_residues): 
            self.distances[resindex] = {"C1'": dist_C1p[resindex], "P": dist_P[resindex] }

        return 

    def _get_A_distance(self, atomname, n=5):
        d_c = {}
        for res in self.u.residues: 
            resindex = res.resindex
            A = res.atoms.select_atoms("name " + atomname)

            if len(A)==0:
                d_c[resindex] = {"strand": [], "compl": []}
                continue


            dist_A = []
            dist_A_wc = []
            #TODO: -mid- only reasonable neighbors
            for i in range(-n,n+1):
                try:
                    res_x = self.u.residues[resindex + i]
                    B = res_x.atoms.select_atoms("name " + atomname)
                    dist = np.linalg.norm(A.positions - B.positions)
                except IndexError:
                    dist = None
                dist_A.append(dist)
           
            try:
                resindex_wc = self.d_bp[resindex]
                for i in range(-n,n+1):
                    try:
                        res_x = self.u.residues[resindex_wc + i]
                        B = res_x.atoms.select_atoms("name " + atomname)
                        if len(A)==0:
                            d_c[resindex] = {"strand": dist_A, "compl": []}
                            continue
                        dist = np.linalg.norm(A.positions - B.positions)
                    except IndexError:
                        dist = None
                    dist_A_wc.append(dist)
            except KeyError:
                pass

            d_c[resindex] = {"strand": dist_A, "compl": dist_A_wc}
        return d_c

def print_usage():
    print("""
    initializes MDAnalysis universe and compoutes watson scric base pairs. they are returned as to dictionaries. this process is repeated for each Hbond-deviation criterion
    subsequently universe and dicts are stored into a pickle. each deviation criterion is stored in one pickle

    usage: projectname designname frames ####[dev = 0.1] [dev2] [dev3] ... 

    return: creates a pickle for each deviation: the pickle contains: (top, trj), wc_pairs, wc_index_pairs
            the tuple contains the absolute path of the files md-files (universe cannot be pickled), second and third are the two dictionaries
        """)


def proc_input():
    
    
    if len(sys.argv) < 3:
        print_usage()
        exit(0)

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
    except  FileExistsError:
        pass

    n_frames = int(sys.argv[3])

    dev = []
    if len(sys.argv) == 4:
        dev.append(0.1)

    for av in sys.argv[4:]:
        dev.append(float(av))

    return top, trj, base, name, output, dev, n_frames

def write_pdb(u, bDNA, PDBs):
    u.add_TopologyAttr(mda.core.topologyattrs.Tempfactors(np.zeros(len(u.atoms))))
    
    u.atoms.tempfactors = -1.
    for res in u.residues:
        try:
            res.atoms.tempfactors = bDNA.wc_quality[res.resindex]["C1'C1'"]
        except KeyError:
            pass
    PDBs["qual"].write(u.atoms)
    
    for cond in ["rise", "twist"]:
        u.atoms.tempfactors = -1.
        for res in u.residues:
            try:
                res.atoms.tempfactors = bDNA.wc_geometry[res.resindex][cond]["center"]
            except KeyError:
                pass
        PDBs[cond].write(u.atoms)
            
    u.atoms.tempfactors = -1.
    ing = 0.00
    for resindex, resindex_wc in bDNA.d_bp.items():
        u.residues[resindex].atoms.tempfactors = ing
        u.residues[resindex_wc].atoms.tempfactors = ing
        ing += 0.01
    PDBs["bp"].write(u.atoms)

    return

def main():
    top, trj, base, name ,output, deviations, n_frames = proc_input()
    print("read ", top, trj, deviations, n_frames)
    print("output to ", output)

    # initialize universe and select final frame
    u = mda.Universe(top, trj)

    frames_step = int( len(u.trajectory) /  n_frames)
    frames = list(range(len(u.trajectory)-1,0,-frames_step))

    linker = bpLinker.Linker( base + name )
    dict_bp, dict_idid, dict_hpid =  linker.link()
    
    dicts_tuple = [(dict_bp,"bp"), (dict_idid,"idid"), (dict_hpid,"hpid"), ((top, trj),"universe")]
    for dict_, dict_name in dicts_tuple:
        pickle.dump(dict_, open( output + name + "__" + dict_name + "-dict.p", "wb"))
    
    properties = []
    traj_out = output + "frames/"
    try:     
        os.mkdir(traj_out)
    except  FileExistsError:
        pass

    # open PDB files
    PDBs = {}
    for pdb_name in ["bp", "rise", "twist", "qual"]:
        PDBs[pdb_name] = mda.Writer(output + name + "__wc_" + pdb_name + ".pdb", multiframe=True)
  
    # loop over selected frames
    for i, ts in enumerate([u.trajectory[i] for i in frames]):
        print(ts)
        bDNA = BDna(u, dict_bp, dict_idid, dict_hpid)
        
        #perform analyis
        bDNA.eval_wc()
        bDNA.eval_distances()
        bDNA.eval_dh()
        properties.append(bDNA)
        props_tuple = [(bDNA.wc_geometry,"__bDNA-wc_geometry-"), (bDNA.wc_quality,"__bDNA-wc_quality-"), 
            (bDNA.dh_quality,"__bDNA-wc_quality-"), (bDNA.distances,"__bDNA-distances-")]
        for prop, prop_name in props_tuple:
            pickle.dump((ts, prop), open( traj_out + name + prop_name + str(i) + ".p", "wb"))
        write_pdb(u, bDNA, PDBs)
    
    # close PDB files
    for _,PDB in PDBs.items():
        PDB.close()
    
    # dosomething with properties
    #TODO: -low-
    return



if __name__ == "__main__":
    main()
