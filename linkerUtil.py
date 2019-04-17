import mrcfile
import numpy as np
import ipywidgets as widgets
import MDAnalysis as mda
import nglview as nv
from IPython.display import display


def mrc(atoms_selection, path, context=2):
    
    if not len(atoms_selection):
        print("EXIT - no atoms in this selection")
        return
        
    
    
    with mrcfile.open( path +".mrc", mode='r+')  as mrc:
        m_o = np.array(mrc.header["origin"])
        m_origin = np.array([m_o["x"],m_o["y"],m_o["z"]])
        m_c = np.array(mrc.header["cella"])
        m_cell = np.array([m_c["x"],m_c["y"],m_c["z"]])
        m_grid = np.array([mrc.header["nz"], mrc.header["ny"],mrc.header["nx"]])
        m_spacing = m_cell/m_grid
        m_data = mrc.data


    data_mask = np.zeros(m_grid, dtype=np.float32)

    grid_positions = np.rint(((atoms_selection.positions-m_origin)/ m_spacing)).astype(int)
    for pos in grid_positions:
        x, y, z = pos[2], pos[1], pos[0] # fast to slow axis
        data_mask[x-context:x+context,y-context:y+context,z-context:z+context] = 1.
        
    data = m_data * data_mask

    
    
    #get rid of zero-padding
    idx_data = np.nonzero(data)

    x_min, x_max = np.min(idx_data[0]),  np.max(idx_data[0])
    y_min, y_max = np.min(idx_data[1]),  np.max(idx_data[1])
    z_min, z_max = np.min(idx_data[2]),  np.max(idx_data[2])


    data_small = np.zeros((x_max-x_min,  y_max-y_min, z_max-z_min), dtype=np.float32)
    data_small = data[x_min: x_max,y_min: y_max,z_min: z_max] 


    grid = np.shape(data_small)
    origin = m_origin + (z_min, y_min, x_min) * m_spacing
    cell = grid * m_spacing


    with mrcfile.new(path +"-temp.mrc", overwrite=True) as mrc_out:
        mrc_out.set_data(data_small)
        mrc_out._set_voxel_size(*(cell/grid))
        mrc_out.header["origin"] = tuple(origin)

    return


def fit_widget(atoms_selection, atoms_scaffold, atoms_staple, path, opacity_map=0.9, isolevel=5., backbone=True, fast=False):
    
    #TODO get cadnano staple colors not genreic ones
    colors= ["red","green","black", "brown", "gold", "navy", "teal", "orange", "pink"]
    
    if not len(atoms_selection):
        print("NOTHING TO DISPLAY - no atoms in this selection")
        return
    
    if backbone:
        scaffold_view = atoms_scaffold.select_atoms("name C1' C2' C3' C4' C5' O1P O2P O3' O4' O5' P")
        staple_view = atoms_staple.select_atoms("name C1' C2' C3' C4' C5' O1P O2P O3' O4' O5' P")

    else:
        scaffold_view = atoms_scaffold
        staple_view = atoms_staple
    
    view = nv.show_mdanalysis(scaffold_view)
    
    view.clear(component = 0 )
    view.add_representation("spacefill",component = 0,  radius = 0.6 , color="blue")#add_representation("spacefill", component = 0,  radius = 0.6 , color="blue")


    view.add_component(path +"-temp.mrc")
    view.clear(component = 1 )
    view.add_representation("surface", component = 1, color='grey', wireframe=True, opacity=opacity_map, isolevel=isolevel)
    if fast: #no staple coloring
            view.add_trajectory(staple_view)
            view.clear(component = 2 )
            view.add_representation("spacefill", component =2, radius = 0.6 , color="red")
    else:    
        segnames = [seg.segid  for seg in staple_view.segments]
        for i, segid in enumerate(segnames):
            s_view = staple_view.select_atoms("segid "+str(segid))
            view.add_trajectory(s_view)
            view.clear(component = (i+ 2) )
            view.add_representation("spacefill", component =(i+ 2), radius = 0.6 , color=colors[i % len(colors)])

    display(view)
    
    return