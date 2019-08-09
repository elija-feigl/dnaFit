#!/usr/bin/env python3#
import sys
import os

import numpy as np
import pandas as pd
import pickle
from statistics import mean

import ipywidgets as widgets

# TODO: -low get dynamically from dicts
DICTS = ["dict_bp", "dict_idid", "dict_hpid", "dict_co", "dict_nicks",
         "list_skips", "dict_idseq"]
CATEGORIES = ["co", "co_plus", "ss", "ds", "all", "clean", "nick", "A", "T",
              "G", "C"]
CO_CATEGORIES = ["single", "double", "end", "all"]
PROP_TYPE = ["wc_geometry", "wc_quality",
             "dh_quality", "distances"]
WCGEOMETRY = ["twist", "rise", "tilt", "roll", "shift", "slide"]
DIST = ["C1'", "P"]
DIHEDRALS = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "xi"]
COANGLES = ["co_alpha1", "co_alpha2", "co_gamma1", "co_gamma2", "co_beta"]


class DataPrep(object):
    def __init__(self, path, name, plus=3):
        self.name = name
        self.path = path
        self.topo = self._topology()
        self.categories = self._categorise(plus)
        self.columns = ["categories", "position"]

        self.df = None
        self.df_co = None

    def _topology(self):
        topo = {}
        for pickle_name in DICTS:
            topo[pickle_name] = pickle.load(
                open(self.path + self.name + "__" + pickle_name + ".p", "rb"))

        self.inv_dict_idid = {v: k for k, v in topo["dict_idid"].items()}
        self.inv_dict_hpid = {v: k for k, v in topo["dict_hpid"].items()}
        self.inv_dict_bp = {v: k for k, v in topo["dict_bp"].items()}

        return topo

    def _traj_frame(self, frame):
        traj_path = self.path + "frames/"
        data = {}
        for pickle_name in PROP_TYPE + ["co_angles"]:
            nnn = (traj_path + self.name + "__bDNA-" +
                   pickle_name + "-" + str(frame) + ".p")
            ts, prop = pickle.load(open(nnn, "rb"))
            data[pickle_name] = prop
        return data, ts

    def _categorise(self, plus):
        # TODO: -high sequence
        id_A = set(resindex for resindex, base in
                   self.topo["dict_idseq"].items() if base == "A")
        id_T = set(resindex for resindex, base in
                   self.topo["dict_idseq"].items() if base == "T")
        id_G = set(resindex for resindex, base in
                   self.topo["dict_idseq"].items() if base == "G")
        id_C = set(resindex for resindex, base in
                   self.topo["dict_idseq"].items() if base == "C")

        id_ds = set()
        id_co_plus = set()

        for wc_id1, wc_id2 in self.topo["dict_bp"].items():
            id_ds.add(wc_id1)
            id_ds.add(wc_id2)
        id_ss = set(self.topo["dict_idid"].values()) - id_ds

        id_co_bp = set()
        id_co = {id_fit for id_fit in self.topo["dict_co"].keys()
                 if id_fit not in id_ss}
        for resindex in id_co:
            try:
                id_co_bp.add(self.topo["dict_bp"][resindex])
            except KeyError:
                id_co_bp.add(self.inv_dict_bp[resindex])

        id_co = id_co | id_co_bp
        id_ds = id_ds - id_co

        for resindex in id_co:  # TODO: fix quick and dirty
            id_co_plus.add(resindex)
            for i in range(-plus, plus):
                id_co_plus.add(resindex+i)

        id_all = id_ss | id_ds | id_co_plus
        id_clean = id_all - (id_co | id_co_plus | id_ss)
        id_nick = (list(self.topo["dict_nicks"].values()) +
                   list(self.topo["dict_nicks"].keys()))

        return {"co": id_co, "co_plus": id_co_plus, "ss": id_ss, "ds": id_ds,
                "all": id_all, "clean": id_clean, "nick": id_nick,
                "A": id_A, "T": id_T, "G": id_G, "C": id_C}

    def create_df(self, frame=0):
        data, ts = self._traj_frame(frame)
        max_res = max(self.categories["all"])
        id_prop_dict = {}
        for resindex in range(0, max_res+1):
            position = self.inv_dict_hpid[self.inv_dict_idid[resindex]]
            categories = [
                cat for cat in CATEGORIES if resindex in self.categories[cat]]
            strand_type = "scaffold" if position[2] else "staple"
            categories.append(strand_type)
            id_prop_dict[resindex] = [categories, position]

            for prop in PROP_TYPE:
                if prop == "wc_quality":
                    try:  # not all are basepaired
                        bp_qual = mean(data[prop][resindex].values())
                    except KeyError:
                        bp_qual = None
                    id_prop_dict[resindex].append(bp_qual)
                    if prop not in self.columns:
                        self.columns.append(prop)
                elif prop == "dh_quality":
                    for dh in DIHEDRALS:
                        angle = data[prop][resindex][dh]
                        id_prop_dict[resindex].append(angle)
                        if ("dh-" + dh) not in self.columns:
                            self.columns.append("dh-" + dh)
                elif prop == "wc_geometry":
                    for geom in WCGEOMETRY:
                        try:
                            cent_ank = data[prop][resindex][geom]
                        except KeyError:
                            cent_ank = {"C1p": None,
                                        "diazine": None, "C6C8": None}
                        for loc, value in cent_ank.items():
                            id_prop_dict[resindex].append(value)
                            if (geom + "-" + loc) not in self.columns:
                                self.columns.append(geom + "-" + loc)
                elif prop == "distances":
                    for atom in DIST:
                        dist_dicts = data[prop][resindex][atom]
                        for loc, dist_list in dist_dicts.items():
                            idx_self = int(len(dist_list)*0.5)
                            if loc == "strand":
                                if (dist_list[idx_self+1] is not None and
                                        dist_list[idx_self-1] is not None):
                                    dist = 0.5 * (dist_list[idx_self-1] +
                                                  dist_list[idx_self+1])
                                elif dist_list[idx_self+1] is not None:
                                    dist = dist_list[idx_self+1]
                                else:
                                    dist = dist_list[idx_self-1]
                            else:
                                dist = dist_list[idx_self]

                            id_prop_dict[resindex].append(dist)
                            if (atom + "-" + loc) not in self.columns:
                                self.columns.append(atom + "-" + loc)
        self.df = pd.DataFrame.from_dict(
            id_prop_dict, orient='index', columns=self.columns)

        max_co = max(data["co_angles"].keys())
        id_co_dict = {}
        for co_id in range(1, max_co+1):
            try:
                co_type = data["co_angles"][co_id]["type"]
            except KeyError:
                continue  # running index not continuos
            strand_type = ("scaffold" if
                           data["co_angles"][co_id]["is_scaffold"]
                           else "staple")
            id_co_dict[co_id] = [[co_type, strand_type]]

            for name in COANGLES:
                try:
                    angle = data["co_angles"][co_id]["angles"][name]
                except KeyError:
                    angle = None
                id_co_dict[co_id].append(angle)

        self.df_co = pd.DataFrame.from_dict(
            id_co_dict, orient='index', columns=(["type"] + COANGLES))
        return self.df, self.df_co, ts


class FileBrowser(object):
    """
    https://gist.github.com/DrDub/6efba6e522302e43d055
    """

    def __init__(self):
        self.path = os.getcwd()
        self._update_files()

    def _update_files(self):
        self.files = list()
        self.dirs = list()
        if(os.path.isdir(self.path)):
            for f in os.listdir(self.path):
                ff = self.path + "/" + f
                if os.path.isdir(ff):
                    self.dirs.append(f)
                else:
                    self.files.append(f)

    def widget(self):
        box = widgets.VBox()
        self._update(box)
        return box

    def _update(self, box):

        layout_button = widgets.Layout(
            width='800px', height='30px', text_align="start",
            border="1px solid black")

        def on_click(b):
            if b.description == '..':
                self.path = os.path.split(self.path)[0]
            else:
                self.path = self.path + "/" + b.description
            self._update_files()
            self._update(box)

        buttons = []
        if self.files:
            button = widgets.Button(description='..',
                                    background_color='#d0d0ff',
                                    layout=layout_button)
            button.on_click(on_click)
            buttons.append(button)
        for f in self.dirs:
            button = widgets.Button(description=f,
                                    background_color='#d0d0ff',
                                    layout=layout_button)
            button.on_click(on_click)
            buttons.append(button)
        # for f in self.files:
        #    button = widgets.Button(description=f)
        #    button.on_click(on_click)
        #    buttons.append(button)
        box.children = tuple(
            [widgets.HTML("<h3>%s</h3>" % (self.path,))] + buttons)


def proc_input():
    if len(sys.argv) < 2:
        sys.exit(0)

    project = sys.argv[1]
    name = sys.argv[2]
    cwd = os.getcwd()

    output = cwd + "/" + project + "/analysis/"

    return output, name


def main():

    # process input
    path, name = proc_input()
    print("output to ", path)
    prep = DataPrep(path, name, plus=3)
    df, df_co, ts = prep.create_df()
    # import ipdb; ipdb.set_trace()

if __name__ == "__main__":
    main()
