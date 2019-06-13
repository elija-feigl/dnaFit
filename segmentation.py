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


def print_usage():
    print("""
    initializes MDAnalysis universe and compoutes watson scric base pairs. they are returned as to dictionaries. this process is repeated for each Hbond-deviation criterion
    subsequently universe and dicts are stored into a pickle. each devia, tion criterion is stored in one pickle

    usage: designname [context = 2]  ... 

    return: creates a pickle for each deviation: the pickle contains: (top, trj), wc_pairs, wc_index_pairs
            the tuple contains the absolute path of the files md-files (universe cannot be pickled), second and third are the two dictionaries
        """)


def proc_input():
    if len(sys.argv) < 1:
        print_usage()
    # if isinstance(sys.argv[2], str):
    name = sys.argv[1]
    cwd = os.getcwd()
    path = cwd + "/"

    if len(sys.argv) == 2:
        context = 2
    else:
        context = int(sys.argv[2])

    return path, name, context


def main():

    path_in, name, context = proc_input()
    print("output to ", path_in)

    path_analysis = path_in + "/analysis/"

    input_pickle = ["dict_bp", "dict_idid",
                    "dict_hpid", "dict_coid", "universe"]
    input_linkage = {}
    for pickle_name in input_pickle:
        input_linkage[pickle_name] = pickle.load(
            open(path_analysis + name + "__" + pickle_name + ".p", "rb"))

    # initialize universe and select final frame
    u = mda.Universe(*input_linkage["universe"])
    u.trajectory[-1]

    #full mapp
    mrc_segment(u.atoms, path_in + name, path_analysis + name + "-all", context=context)

    #crossovers
    path_out = path_analysis + "/seg-co/"
    try:     
        os.mkdir(path_out)
    except  FileExistsError:
        pass
    intig = 0
    for id_design, co in input_linkage["dict_coid"].items():
        co_set = set()
        co_context = 6
        id_co_design = co["co"]

        co_set.add(input_linkage["dict_idid"][id_design])
        co_set.add(input_linkage["dict_idid"][id_co_design])
        co_set.add(input_linkage["dict_bp"][   input_linkage["dict_idid"][id_design]])
        co_set.add(input_linkage["dict_bp"][   input_linkage["dict_idid"][id_co_design]])

        co_set_fit_add = set()
        for co in co_set: #TODO: -high- better neighbor search
            for i in range(-co_context,co_context):
                co_set_fit_add.add(co+i)

        co_set_fit = co_set | co_set_fit_add


        #ipdb.set_trace()
        atoms_select = mda.AtomGroup([], u)
        for base_id in co_set_fit:
            atoms_select += u.residues[base_id].atoms
        mrc_segment(atoms_select, path_in + name, path_out + name + "-co" +str(intig), context=context)
        intig += 1
    

if __name__ == "__main__":
    main()
