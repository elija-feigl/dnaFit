#!/usr/bin/env python
import mrcfile
import numpy as np
import MDAnalysis as mda
import sys
import pickle
import os
import ipdb


def mrc_segment(atoms_selection, path_in, path_out, context=2, clipping=0.):

    if not len(atoms_selection):
        print("EXIT - no atoms in this selection")
        return

    u = atoms_selection.universe
    u.trajectory[-1]

    print(mrcfile)
    with mrcfile.open(path_in + ".mrc", mode='r+') as mrc:
        m_o = np.array(mrc.header["origin"])
        m_origin = np.array([m_o["x"], m_o["y"], m_o["z"]])
        m_c = np.array(mrc.header["cella"])
        m_cell = np.array([m_c["x"], m_c["y"], m_c["z"]])
        m_grid = np.array(
            [mrc.header["nz"], mrc.header["ny"], mrc.header["nx"]])
        m_spacing = m_cell/m_grid
        m_data = mrc.data

    data_mask = np.zeros(m_grid, dtype=np.float32)

    grid_positions = np.rint(
        ((atoms_selection.positions-m_origin) / m_spacing)).astype(int)
    for pos in grid_positions:
        x, y, z = pos[2], pos[1], pos[0]  # fast to slow axis
        data_mask[x-context:x+context, y-context:y+context, z-context:z +
                  context] = 1.

    data = m_data * data_mask

    # get rid of zero-padding
    idx_data = np.nonzero(data)

    x_min, x_max = np.min(idx_data[0]),  np.max(idx_data[0])
    y_min, y_max = np.min(idx_data[1]),  np.max(idx_data[1])
    z_min, z_max = np.min(idx_data[2]),  np.max(idx_data[2])

    data_small = np.zeros(
        (x_max-x_min,  y_max-y_min, z_max-z_min), dtype=np.float32)
    data_small = data[x_min: x_max, y_min: y_max, z_min: z_max]

    grid = np.shape(data_small)
    origin = m_origin + (z_min, y_min, x_min) * m_spacing
    cell = grid * m_spacing

    # clipp below threshold
    data_small[data_small < clipping] = 0.

    with mrcfile.new(path_out + ".mrc", overwrite=True) as mrc_out:
        mrc_out.set_data(data_small)
        mrc_out._set_voxel_size(*(cell/grid))
        mrc_out.header["origin"] = tuple(origin)

    return

def proc_input():
    if len(sys.argv) < 2:
        print_usage()
    name = sys.argv[1]
    cwd = os.getcwd()
    path = cwd + "/"
    print(len(sys.argv) )

    if len(sys.argv) > 2:
        rang = int(sys.argv[2])
    else:
        rang = 3

    if len(sys.argv) > 3:
        context = int(sys.argv[3])
    else:
        context = 2

    return path, name, rang, context


def _categorise_lists(topo, plus=3):
    # TODO. names set list etc
    inv_dict_bp = {v: k for k, v in topo["dict_bp"].items()}
    id_ds = set()
    id_coplus = set()

    for wc_id1, wc_id2 in topo["dict_bp"].items():
        id_ds.add(wc_id1)
        id_ds.add(wc_id2)
    id_ss = set(topo["dict_idid"].values()) - id_ds

    id_co = set()
    id_co_init = {
        id_design for id_design in topo["dict_co"].keys() if id_design not in id_ss}
    for base in id_co_init:
        co = topo["dict_co"][base]["co"]
        try:
            co_bp = topo["dict_bp"][co]
        except KeyError:
            co_bp = inv_dict_bp[co]
        try:
            bp = topo["dict_bp"][base]
        except KeyError:
            bp = inv_dict_bp[base]

        if topo["dict_co"][base]["type"][0] == "double":
            dou = topo["dict_co"][base]["type"][1]
            dou_co = topo["dict_co"][dou]["co"]
            try:
                dou_co_bp = topo["dict_bp"][dou_co]
            except KeyError:
                dou_co_bp = inv_dict_bp[dou_co]
            try:
                dou_bp = topo["dict_bp"][dou]
            except KeyError:
                dou_bp = inv_dict_bp[dou]
            tup = (base, bp, co, co_bp, dou, dou_bp, dou_co, dou_co_bp)
        else:
            tup = (base, bp, co, co_bp)

        tup_plus = []
        for x in tup:
            # TODO: fix quick and dirty: wrong if staple ends or another co
            for i in range(-plus, plus):
                tup_plus.append(x+i)

        id_co.add(tup)
        id_coplus.add(tuple(tup_plus))

    #ipdb.set_trace()
    return id_co, id_coplus


def _topology(name, path):
    DICTS = ["dict_bp", "dict_idid", "dict_hpid", "dict_co", "universe"]
    # read general info
    topo = {}
    for pickle_name in DICTS:
        topo[pickle_name] = pickle.load(
            open(path + name + "__" + pickle_name + ".p", "rb"))

    return topo


def print_usage():
    print("""
    initializes MDAnalysis universe and compoutes watson scric base pairs. they are returned as to dictionaries. this process is repeated for each Hbond-deviation criterion
    subsequently universe and dicts are stored into a pickle. each devia, tion criterion is stored in one pickle

    usage: designname [range = 3] [context = 2]  ... 

    return: creates a pickle for each deviation: the pickle contains: (top, trj), wc_pairs, wc_index_pairs
            the tuple contains the absolute path of the files md-files (universe cannot be pickled), second and third are the two dictionaries
        """)

def main():

    path_in, name, rang, context = proc_input()
    print("output to ", path_in)

    path_analysis = path_in + "/analysis/"

    topo = _topology(name, path_analysis)
    _, id_coplus_lists = _categorise_lists(topo, plus=rang)

    # initialize universe and select final frame
    u = mda.Universe(*topo["universe"])
    u.trajectory[-1]

    # full mapp
    mrc_segment(u.atoms, path_in + name, path_analysis +
                name + "-all", context=context)

    # crossovers
    path_out = path_analysis + "/seg-co/"
    try:
        os.mkdir(path_out)
    except FileExistsError:
        pass

    intig = 0
    for co_select in id_coplus_lists:
 
        atoms_select = mda.AtomGroup([], u)
        for base_id in co_select:
            atoms_select += u.residues[base_id].atoms
        mrc_segment(atoms_select, path_in + name, path_out +
                    name + "-co" + str(intig), context=context)
        intig += 1


if __name__ == "__main__":
    main()
