import __future__ 

import ipywidgets as widgets
import MDAnalysis as mda
import numpy as np

import mrcfile
import nanodesign as nd 
from nanodesign.converters import Converter

class Linker(object):
    def __init__(self, path):
        self.fit = Fit(path)
        self.design = Design(path)
        self.selection_scaffold = None
        self.selection_staples = None
        
    
    def _parse_selection(self, base_selection, helix_selection):
        # get helices from design
        helices = []
        for idx in helix_selection:
            helix = self.design.design.structure_helices_map[idx]
            helices.append(helix)
        # get bases from helices
        self.selection_scaffold = []
        self.selection_staples = []
        for helix in helices:
            
            selection_scaffold_add = [ base for base in helix.scaffold_bases if base.p in base_selection and base.num_deletions == 0 and base.num_insertions == 0]
            selection_staples_add = [ base for base in helix.staple_bases if base.p in base_selection and base.num_deletions == 0 and base.num_insertions == 0]

            self.selection_scaffold += selection_scaffold_add
            self.selection_staples += selection_staples_add

        return 

    def _link_scaffold(self):
        """idea: collect position in scaffold (0-x) by comparing to index in list of scaffold_design positions"""
        
        idx_bases = []
        #get all possible positions for a scaffold base
        design_idx = [s.id for s in self.design.scaffold if s.num_deletions == 0 and s.num_insertions == 0]
        for base in self.selection_scaffold:
            idx = design_idx.index(base.id)
            idx_bases.append(idx)

        atoms = self.fit.scaffold.residues[idx_bases].atoms

        return atoms


    def _link_staples(self):
        """ loop over all staples, then perform the same procesdure as for scaffold
        """
        
        atoms = mda.AtomGroup([],self.fit.u)
        
        for i, staple in enumerate(self.design.staples):
            idx_segment = self.design.s_dict[i]
            idx_staple = self._idx_incl(i) 

            
            # get all bases indices of that staple that correspond to selection
            selection =  [b for b in self.selection_staples if  b.strand == idx_staple and b.num_deletions == 0 and b.num_insertions == 0]
            staple_selected = (len(selection) != 0)
            
            if staple_selected:
                idx_bases = []
                #get all possible positions for a  base in this specific staple
                design_idx = [s.id for s in staple if s.num_deletions == 0 and s.num_insertions == 0]

                for base in selection:
                    idx = design_idx.index(base.id)
                    idx_bases.append(idx)
                
                
                #TODO: cleanup after final testing
                for j in idx_bases:
                    try:    
                        atoms += self.fit.staples[idx_segment].residues[j].atoms
                    except:
                        print("something went wrong - excluded base ", j, "on helix ", base.h)
                        pass

        return atoms



    def link(self, helices, bases):

        self._parse_selection(bases, helices)

        atoms_scaffold = self._link_scaffold()
        atoms_staple = self._link_staples()

        atoms_selection = atoms_scaffold + atoms_staple  

        return atoms_selection, atoms_scaffold, atoms_staple  
    


    def _idx_incl(self, idx):
        """ include scaffold from idcount
        """
        if idx >= self.design.excl:
            st_idx = idx + 1 
        else:
            st_idx = idx

        return st_idx
    

    def make_selection(self):
        
        layout_w = widgets.Layout(width= '500px', height = '70px')
        layout_h = widgets.Layout(width= '50px')
        layout_label = widgets.Layout(width= '500px')
        layout_box = widgets.Layout(width= '600px', border= '1px solid black')
        layout_button = widgets.Layout(width='30px', height='30px')

    

        base_p = [base.p for strand in self.design.strands for base in strand.tour]
        maxb = max(base_p)
        minb = min(base_p)
        row = [coor[0] for coor, _ in self.design.design.structure_helices_coord_map.items()]
        minr = min(row)
        maxr = max(row)
        col = [coor[1] for coor, _ in self.design.design.structure_helices_coord_map.items()]
        minc = min(col)
        maxc = max(col) 







        helix_buttons = []
        for r in range(minr, maxr+1):
            row = []
            for c in range(minc, maxc+1):
                try:
                    id = str(self.design.design.structure_helices_coord_map[(r,c)].id)
                except KeyError:
                    id=""

                but = widgets.ToggleButton(description=id, layout=layout_button)
                row.append(but)
            helix_buttons.append(row)



        HBoxes = []
        for row in helix_buttons:
            HBoxes.append(widgets.HBox(row))

        buttons_h = widgets.VBox(HBoxes)
        
        
        slider_b = widgets.IntRangeSlider(
            value=[minb, maxb],
            min=minb,
            max=maxb,
            step=1,
            description='base:',
            disabled=False,
            continuous_update=True,
            orientation='horizontal',
            readout=True,
            readout_format='d',
            layout = layout_w
        )
        
        slider_c = widgets.IntSlider(
            value=2,
            min=1,
            max=10,
            step=1,
            description='context:',
            disabled=False,
            continuous_update=True,
            orientation='vertical',
            readout=True,
            readout_format='d',
            layout = layout_h
        )

        sliders_w = widgets.VBox([buttons_h, slider_b])
        sliders = widgets.Box([sliders_w, slider_c], layout= layout_box)
        display(sliders)



        return (helix_buttons, slider_b, slider_c)
    
    
    def eval_selection(self, helix_buttons, slider_b, slider_c):
    
        
        selection_bases  = range(slider_b.lower, slider_b.upper+1)
        
        
        selection_helices = []
        flat_list = [item for sublist in helix_buttons for item in sublist]
        for button in flat_list:
            if button.value:
                try:
                    selection_helices.append(int(button.description))
                except ValueError:
                    pass
        return (selection_helices, selection_bases),slider_c.value

    def writepdb(self, atoms_selection, path):
        with mda.Writer(path +"-temp.pdb", multiframe=True,  n_atoms=atoms_selection.n_atoms) as W:
            for ts in self.fit.u.trajectory:
                W.write(atoms_selection)
        return
    



class Fit(object):
    
    def __init__(self, path):
        self.path = path
        self.u = self._get_universe()
        self.strands = self.u.segments
        self.scaffold = self.strands[0]
        self.staples = self.strands[1:]

    def _get_universe(self):
        top = self.path + ".psf"
        trj = self.path  + "-out.dcd"
        
        u = mda.Universe(top, trj)
        u.trajectory[-1]
        return u

class Design(object):
        
    def __init__(self, path):
        self.path = path
        self.design = self._get_design()
        self.strands = self.design.strands
        self.scaffold = [ strand.tour for strand in self.strands if strand.is_scaffold ][0]
        self.excl = self.scaffold[0].strand #needed to correc tstaple numbering
        self.staples = [ strand.tour for strand in self.strands if not strand.is_scaffold ]
        self.h_dict = self._create_helix_order()
        self.s_dict = self._create_staple_order()
    
    
    
    def _get_design(self):
        file_name =  self.path + ".json"
        seq_file = self.path + ".seq"
        seq_name = None
        
        converter = Converter()
        converter.read_cadnano_file(file_name, seq_file, seq_name)
        return converter.dna_structure
    
    
    def _create_helix_order(self):
        """ helices are not prcessed in the same order as they are listed by idx
        """
        
        helices_dict = self.design.structure_helices_map
        h_dict = {i:h.load_order for (i, h) in helices_dict.items()}

        return h_dict
    
   

    def _create_staple_order(self):
        """ exchange staple id with load_id 
            map design-staple-order to universe-staple-order
            staple order is not the same for mda and nd
            enrg/MDanalysis: sort left to right top to bottom; only counting 3' ends
            nanodesigns/auto: sort left to right top to bottom; counting every staple piece
        """        
        list_hp = [(self.h_dict[s[0].h], s[0].p) for s in self.staples]
        list_hp_s = sorted(list_hp, key=lambda x: (x[0], x[1]))

        idx_list = [list_hp.index(list_hp_s[i]) for i in range(len(list_hp))]

        s_dict = dict(zip( idx_list, range(len(idx_list))))

        return s_dict
    