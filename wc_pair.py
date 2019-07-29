#!/usr/bin/env python3
import MDAnalysis as mda
import numpy as np
import sys
import os

import ipdb
from MDAnalysis.lib import mdamath

import pickle
import bpLinker

#TODO: -mid- DOC
#TODO: -mid- move CONSTANTS

C1P_BASEDIST = 10.7
TOL = 10e6
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

WC_PROPERTIES = ["rise", "slide", "shift", "twist", "tilt", "roll"]

def _norm(vector):
    return vector/ np.linalg.norm(vector)
def _proj(u, v):
    return np.inner(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
def _v_proj(u, v):
    return np.inner(u, v) / (np.linalg.norm(v) * np.linalg.norm(v)) * v

#TODO: -mid- eval nicks

class BDna(object):

    def __init__(self, universe, d_Fbp, d_DidFid, d_DhpsDid, d_Fco):
        self.u = universe
        self.d_DidFid = d_DidFid
        self.d_Fbp = d_Fbp
        self.d_Fbp_full =  {**d_Fbp, **{v: k for k, v in d_Fbp.items()}}
        self.d_DhpsDid =d_DhpsDid
        self.d_Fco = d_Fco  

        self.wc_quality = None  
        self.wc_geometry = None   
        self.dh_quality = None  
        self.distances = None    #only scaffold bases have complement
        self.co_angles = None

    def _get_next_wc(self, resindex, resindex_wc):
        """ get next residue and its complemt wc pair
            check if next exists.
            check has bp (ony scaffold residues apply here)
        """
        n_resindex = resindex + 1
        
        try:
            self.u.residues[n_resindex]
            try:
                n_resindex_wc = self.d_Fbp_full[n_resindex]
            except KeyError:
                n_resindex_wc = None
        except KeyError:
            n_resindex = None
            n_resindex_wc = None

        return n_resindex, n_resindex_wc 
    
    def eval_wc(self):
        """  dict of wc-quality for each residue (key = resindex)
        """
        self.wc_quality = {}
        self.wc_geometry = {}
        for resindex, resindex_wc in self.d_Fbp.items():      #scaffold only -> just once

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
        """  dict of geometry of basepairs for each wc-pair (key= resid-scaffold)
        """
        def _get_twist(basepair, n_basepair):
            twist = []
            for direct in ["dir-anker", "dir-center"]:
                v1 = basepair[direct]
                v2 = n_basepair[direct]
                twist.append(np.rad2deg(np.arccos(_proj(v1,v2))))
             
            return {"anker":twist[0], "center":twist[1]}

        def _get_rise(basepair, n_basepair):
            rise = []
            n0 = basepair["n0"]
            for direct in ["center-anker", "center"]:
                P1 = basepair[direct]
                P2 = n_basepair[direct]
                rise.append(np.abs(np.inner((P2 - P1), n0)))

            return {"anker":rise[0], "center":rise[1]}

        def _get_shift(basepair, n_basepair):
            shift = []
            for direct in [("center-anker","dir-anker") , ("center","dir-center") ]:
                n0 = _norm(np.cross(basepair["n0"], basepair[direct[1]]))
                P1 = basepair[direct[0]]
                P2 = n_basepair[direct[0]]
                shift.append(np.inner((P2 - P1), n0))

            return {"anker":shift[0], "center":shift[1]}

        def _get_slide(basepair, n_basepair):
            slide = []
            for direct in [("center-anker","dir-anker") , ("center","dir-center") ]:
                n0 = _norm(basepair[direct[1]])
                P1 = basepair[direct[0]]
                P2 = n_basepair[direct[0]]
                slide.append(np.inner((P2 - P1), n0))

            return {"anker":slide[0], "center":slide[1]}

        def _get_tilt(basepair, n_basepair):
            tilt = []
            n0 = basepair["n0"]
            n_n0 = n_basepair["n0"]

            for direct in ["dir-anker", "dir-center"]:
                rot_axis = np.cross(basepair[direct], n0)   
                projn0 = n0 - _v_proj(n0, rot_axis)
                projn_n0 = n_n0 - _v_proj(n_n0, rot_axis)

                dist = _proj(projn0,projn_n0)
                if 1.0 < abs(dist) < 1.0 + TOL : dist = np.sign(dist)
                if dist > 0: 
                    tilt.append(np.rad2deg(np.arccos( dist )))
                else:
                    tilt.append(np.rad2deg(- np.arccos( abs(dist) )))

            for value in tilt:
                if value > 90.: value - 180.
             
            return {"anker":tilt[0], "center":tilt[1]}
        
        def _get_roll(basepair, n_basepair):
            roll = []
            n0 = basepair["n0"]
            n_n0 = n_basepair["n0"]

            for direct in ["dir-anker", "dir-center"]:
                rot_axis = basepair[direct]
                projn0 = n0 - _v_proj(n0, rot_axis)
                projn_n0 = n_n0 - _v_proj(n_n0, rot_axis)

                dist = _proj(projn0,projn_n0)
                if  1.0 < abs(dist) < 1.0 + TOL: dist = np.sign(dist)
                if dist > 0: 
                    roll.append(np.rad2deg(np.arccos( dist )))
                else:
                    roll.append(np.rad2deg(- np.arccos( abs(dist) )))
             
            return {"anker":roll[0], "center":roll[1]}



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

        geometry = {"rise": _get_rise(*basepairs), 
                "slide": _get_slide(*basepairs),
                "shift": _get_shift(*basepairs),
                "twist": _get_twist(*basepairs),
                "tilt": _get_tilt(*basepairs),
                "roll": _get_roll(*basepairs)}
        return geometry, geometry

    def _get_base_plane(self, res):

        atom = []
        for atom_name in ["C2", "C4", "C6", "C1'"]:
            A = res.atoms.select_atoms("name " + atom_name)[0]
            atom.append(A.position)

        n0 = _norm(np.cross((atom[1]-atom[0]), (atom[2]-atom[1])))
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
        for dh_name in DH_ATOMS.keys():
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
            
            n1 = _norm(np.cross(v1,v2))
            n2 = _norm(np.cross(v2,v3))
            m1 = np.cross(n1, _norm(v2))
        
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
                resindex_wc = self.d_Fbp[resindex]
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

    def eval_co_angles(self):
        """ Definition of the angles enclosed by the four helical legs of a cross-over. 
            Vectors are computed using the coordinates of base pair midpoints at the 
            cross-over position and 2 bp away from the cross-over in each leg. The
            cross-over vector x is computed from the coordinates of the midpoints
            between the two base pairs in each of the two helices at the cross-over
            position and is normal to what we call the cross-over plane. The subscript
            “jj” indicates vectorial projections into the cross-over plane. The angle β
            is also computed as indicated for vectors C and D.(Bai 2012)
        """
     
        def get_co_baseplanes(res_index, leg_index, co_index, coleg_index):
            bases = []
            for r in [res_index, leg_index, co_index, coleg_index]:
                try:
                    wc_r = self.d_Fbp_full[r]
                    bases.append(self._get_base_plane(self.u.residues[r]))
                    bases.append(self._get_base_plane(self.u.residues[wc_r]))
                except KeyError:
                    bases.append(self._get_base_plane(self.u.residues[r]))
                    bases.append(self._get_base_plane(self.u.residues[r])) 
            bp_planes = []
            for i in [0, 2, 4, 6]:
                bp_planes.append({"center-anker": ((bases[i]["anker"] + bases[i+1]["anker"]) * 0.5),
                        "dir-anker": (bases[i]["anker"] - bases[i+1]["anker"]) ,
                        "center": ((bases[i]["center"] + bases[i+1]["center"]) * 0.5), 
                        "dir-center": (bases[i]["anker"] - bases[i+1]["anker"]) ,
                        "n0": ((bases[i]["n0"] + bases[i+1]["n0"]) * 0.5) })
            return bp_planes

        def get_co_angles_end(bpplanes):  #TODO: -low- cleanup
            
            def project_to_plane(vect, n0): 
                pro = []
                for x in vect: 
                    d = np.inner(x, n0)
                    p = [d * n0[i] for i in range(len(n0))]
                    pro.append([x[i] - p[i] for i in range(len(x))])
                return pro


            a1 = bpplanes[0]["center"]
            a2 = bpplanes[1]["center"]
            c1 = bpplanes[2]["center"]
            c2 = bpplanes[3]["center"]

            center = (a1 + c1) * 0.5
            a = a2 - a1
            c = c2 - c1
            
            n0 = _norm((a1 - c1))
            proj_ac = project_to_plane([a,c], n0)

            #dist = _proj(proj_ac[0], proj_ac[0])
            #if 1.0 < abs(dist) < 1.0 + TOL : dist = np.sign(dist)
            #gamma1 = np.rad2deg(np.arccos(dist))
            ang_temp = np.rad2deg(np.arccos(_proj(a, n0))) #unprojected
            beta = 90. - ang_temp

            return {"angles": {"co_beta": beta}, "center": center, "plane": n0}


        def get_co_angles_full(bpplanes, double_bpplanes): #TODO: -low- cleanup -high- check!

            def project_to_plane(vect, n0): 
                pro = []
                for x in vect: 
                    d = np.inner(x, n0)
                    p = [d * n0[i] for i in range(len(n0))]
                    pro.append([x[i] - p[i] for i in range(len(x))])
                return pro

            points = []
            
            for x in [bpplanes[0:2], double_bpplanes[0:2],bpplanes[2:], double_bpplanes[2:]]: #abcd
                    ins = x[0]["center"]
                    out = x[1]["center"]
                    points.append((ins, out))
        
            # acbd
            center = sum( p[0] for p in points) * 0.25 #((a0 +b0)/2 +  (c0 +d0)/2)/2
            n0 = _norm(points[0][0] + points[2][0] - points[1][0] - points[3][0]  ) # n0((a0 +b0)/2 -  (c0 +d0)/2)

            acbd = []
            for p_in, p_out in points:
                acbd.append(p_out - p_in)            
        
            proj_acbd = project_to_plane(acbd, n0)
            
            d1 = _proj(proj_acbd[0], proj_acbd[1])
            if 1.0 < abs(d1) < 1.0 + TOL : d1 = np.sign(d1)
            gamma1 = np.rad2deg(np.arccos(d1))
            d2 = _proj(proj_acbd[2], proj_acbd[3])
            if 1.0 < abs(d2) < 1.0 + TOL : d2 = np.sign(d2)
            gamma2 = np.rad2deg(np.arccos(d2))
    
            a1 = _proj(proj_acbd[0], [ -x for x in proj_acbd[2]] )
            if 1.0 < abs(a1) < 1.0 + TOL : a1 = np.sign(a1)
            alpha1 = np.rad2deg(np.arccos(a1))
            a2 = _proj(proj_acbd[1], [ -x for x in proj_acbd[3]] )
            if 1.0 < abs(a2) < 1.0 + TOL : a2 = np.sign(a2)
            alpha2 = np.rad2deg(np.arccos(a2))
           
            ang_temp1 = np.rad2deg(np.arccos(_proj(acbd[0], n0))) #unprojected
            ang_temp2 = np.rad2deg(np.arccos(_proj(acbd[2], n0)))
            
        
            beta = 180. - ang_temp1 - ang_temp2

            return { "angles": {"co_beta": beta, "co_gamma1":gamma1, "co_gamma2":gamma2,
                    "co_alpha1":alpha1, "co_alpha2":alpha2}, "center": center, "plane": n0}
            

        self.co_angles = {}
        
        co_done = set()
        for res_index, co in self.d_Fco.items():     
            if res_index not in co_done:
                leg_index = co["leg"]
                co_index = co["co"]
                coleg_index = self.d_Fco[co_index]["leg"]

                co_done.update([res_index, co_index])
                # res -> a, co -> c
                bpplanes = get_co_baseplanes(res_index, leg_index, co_index, coleg_index)

                co_type = co["type"][0]
                if co_type == "double":
                    double_res_index = co["type"][1]
                    double_leg_index = self.d_Fco[double_res_index]["leg"]
                    double_co_index = self.d_Fco[double_res_index]["co"]
                    double_coleg_index = self.d_Fco[double_co_index]["leg"]
                    co_done.update([double_res_index, double_co_index])
                    # double -> b, d_co -> d
                    double_bpplanes = get_co_baseplanes( double_res_index, double_leg_index, double_co_index, double_coleg_index)
                    co_data = get_co_angles_full(bpplanes, double_bpplanes) # acbd
                    crossover_ids = (res_index,  double_res_index, co_index, double_co_index)
                elif co_type == "single":
                    single_res_index = co["type"][1]
                    single_leg_index = co["type"][2]
                    single_co_index = self.d_Fco[co_index]["type"][1]
                    single_coleg_index = self.d_Fco[co_index]["type"][2]
                    co_done.update([single_res_index, single_co_index])
                    # single -> b, s_co -> d
                    single_bpplanes = get_co_baseplanes( single_res_index, single_leg_index, single_co_index, single_coleg_index)
                    co_data = get_co_angles_full(bpplanes, single_bpplanes) # acbd
                    crossover_ids = (res_index,  single_res_index, co_index, single_co_index)
                else:
                    co_data = get_co_angles_end(bpplanes)
                    crossover_ids = (res_index, co_index)
                
                self.co_angles[co["co_index"]] = {"ids(abcd)": crossover_ids, "type": co_type, "is_scaffold": co["is_scaffold"], "angles": co_data["angles"], "center": co_data["center"]} 
        return 

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

    project = sys.argv[1]
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
    
    for cond in WC_PROPERTIES:
        u.atoms.tempfactors = -1.
        for res in u.residues:
            try:
                res.atoms.tempfactors = bDNA.wc_geometry[res.resindex][cond]["center"]
            except KeyError:
                pass
        PDBs[cond].write(u.atoms)

    for dh in DH_ATOMS.keys():
        u.atoms.tempfactors = -1.
        for res in u.residues:
            try:
                res.atoms.tempfactors = bDNA.dh_quality[res.resindex][dh]
            except KeyError:
                pass
        PDBs[dh].write(u.atoms)
            
    u.atoms.tempfactors = -1.
    ing = 0.00
    for resindex, resindex_wc in bDNA.d_Fbp.items():
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
    universe = (top, trj)
    u = mda.Universe(top, trj)

    frames_step = int( len(u.trajectory) /  n_frames)
    frames = list(range(len(u.trajectory)-1,0,-frames_step))

    linker = bpLinker.Linker( base + name)
    dict_bp, dict_idid, dict_hpid, _ = linker.link()
    dict_co = linker.identify_crossover()
    for dict_name in [ "dict_bp", "dict_idid", "dict_hpid", "dict_co", "universe"]:
        pickle.dump(eval(dict_name), open( output + name + "__" + dict_name + ".p", "wb"))
    
    properties = []
    traj_out = output + "frames/"
    try:     
        os.mkdir(traj_out)
    except  FileExistsError:
        pass

    # open PDB files
    PDBs = {}
    for pdb_name in [*WC_PROPERTIES, "bp", "qual" ]:
        PDBs[pdb_name] = mda.Writer(output + name + "__wc_" + pdb_name + ".pdb", multiframe=True)
    for pdb_name in DH_ATOMS.keys():
        PDBs[pdb_name] = mda.Writer(output + name + "__dh_" + pdb_name + ".pdb", multiframe=True)
  
    # loop over selected frames
    for i, ts in enumerate([u.trajectory[i] for i in frames]):
        print(ts)
        bDNA = BDna(u, dict_bp, dict_idid, dict_hpid, dict_co)
        
        #perform analyis
        print("eval_wc", name)
        bDNA.eval_wc()
        print("eval_distances", name)
        bDNA.eval_distances()
        print("eval_dh", name)
        bDNA.eval_dh()
        print("eval_co_angles", name)
        bDNA.eval_co_angles()
        #ipdb.set_trace()
        
        print("write pickle", name)
        properties.append(bDNA)
        props_tuple = [(bDNA.wc_geometry,"wc_geometry"), (bDNA.wc_quality,"wc_quality"), 
            (bDNA.dh_quality,"dh_quality"), (bDNA.distances,"distances"), (bDNA.co_angles,"co_angles")]
        for prop, prop_name in props_tuple:
            pickle.dump((ts, prop), open( traj_out + name + "__bDNA-" + prop_name + "-" + str(i)+ ".p", "wb"))
        print("write pdbs", name)
        write_pdb(u, bDNA, PDBs)
    
    # close PDB files
    for _,PDB in PDBs.items():
        PDB.close()
    
    return


if __name__ == "__main__":
    main()
