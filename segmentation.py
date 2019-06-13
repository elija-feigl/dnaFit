#!/usr/bin/env python
import mrcfile
import numpy as np
import MDAnalysis as mda
import sys
import os

def mrc_ccc(atoms_selection, path, mrc, context=2, clipping=0.):

    if not len(atoms_selection):
        print("EXIT - no atoms in this selection")
        return

    u = atoms_selection.universe
    u.trajectory[-1]

    print(mrcfile)
    with mrcfile.open(mrc + ".mrc", mode='r+') as mrc:
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

    with mrcfile.new(path + "-masked.mrc", overwrite=True) as mrc_out:
        mrc_out.set_data(data_small)
        mrc_out._set_voxel_size(*(cell/grid))
        mrc_out.header["origin"] = tuple(origin)

    return


def print_usage():
    print("""
    initializes MDAnalysis universe and compoutes watson scric base pairs. they are returned as to dictionaries. this process is repeated for each Hbond-deviation criterion
    subsequently universe and dicts are stored into a pickle. each devia, tion criterion is stored in one pickle

    usage: projectname designname folder mrcname [context = 2]  ... 

    return: creates a pickle for each deviation: the pickle contains: (top, trj), wc_pairs, wc_index_pairs
            the tuple contains the absolute path of the files md-files (universe cannot be pickled), second and third are the two dictionaries
        """)


def proc_input():
    if len(sys.argv) < 2:
        print_usage()

    # if isinstance(sys.argv[1], str):
    project = sys.argv[1]
    # if isinstance(sys.argv[2], str):
    name = sys.argv[2]
    cwd = os.getcwd()
    mrc = sys.argv[4]
    folder = sys.argv[3]

    path = cwd + "/" + project + "/" + name
    top = cwd + "/" + project + "/" + name + ".psf"
    trj = cwd + "/" + project + "/" + folder + "/" + name + ".dcd"
    mrcname = cwd + "/" + project + "/" + name + mrc

    if len(sys.argv) == 5:
        context = 2
    else:
        context = int(sys.argv[5])

    return top, trj, path, context, mrcname


def main():

    top, trj, path, context, mrcfile = proc_input()
    print("read ", top, trj, context)
    print("output to ", path)

    # initialize universe and select final frame
    u = mda.Universe(top, trj)
    u.trajectory[-1]

    

    mrc_ccc(u.atoms, path, mrcfile, context=context)


if __name__ == "__main__":
    main()
