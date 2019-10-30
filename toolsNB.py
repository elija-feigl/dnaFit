#!/usr/bin/env python3#
import sys
import os

import pandas as pd
import numpy as np
import pickle
from pathlib import Path

from statistics import mean

import ipywidgets as widgets

from linkage import Linkage
from project import Project

# TODO: -low get dynamically from dicts
CATEGORIES = ["co", "co_plus", "ss", "ds", "clean", "nick"]
CO_CATEGORIES = ["single", "double", "end", "all"]
PROP_TYPE = ["bp_geometry_local", "bp_geometry_global", "bp_quality",
             "dh_quality", "distances", "localres"]
WCGEOMETRY = ["twist", "rise", "tilt", "roll", "shift", "slide"]
DIHEDRALS = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "xi"]
COANGLES = ["co_alpha1", "co_alpha2", "co_gamma1", "co_gamma2", "co_beta"]


class DataPrep(object):
    def __init__(self, path, name, plus=3):
        self.name = name
        self.path = path
        self.link = self._load_linkage()
        self.categories = self._categorise(plus)
        self.columns = ["categories", "position"]

        self.df = None
        self.df_co = None

    def _load_linkage(self):
        self.project = Project(input=Path(self.path),
                               output=Path(self.path) / "analysis",
                               name=self.name,
                               )
        link = Linkage()
        link.load_linkage(project=self.project)
        link.reverse()
        return link

    def _traj_frame(self, frame):
        traj_path = self.project.output / "frames/"
        data = {}
        for pickle_name in PROP_TYPE + ["co_angles"]:
            if pickle_name == "localres":
                nnn = "{}/{}__{}.p".format(self.project.output,
                                           self.name,
                                           pickle_name,
                                           )
                try:
                    prop = pickle.load(open(nnn, "rb"))
                except FileNotFoundError:
                    print("flie {} not found".format(nnn))
                    prop = None
            else:
                nnn = "{}/{}__bDNA-{}-{}.p".format(traj_path,
                                                   self.name,
                                                   pickle_name,
                                                   frame
                                                   )
                ts, prop = pickle.load(open(nnn, "rb"))
            data[pickle_name] = prop
        return data, ts

    def _categorise(self, plus):
        id_ds = set()

        for wc_id1, wc_id2 in self.link.Fbp.items():
            id_ds.add(wc_id1)
            id_ds.add(wc_id2)
        id_ss = set(self.link.DidFid.values()) - id_ds

        id_co = set()
        for co in self.link.Fco.values():
            for P in co.Ps:
                if P is None:
                    continue
                if P.sc is not None:
                    id_co.add(P.sc.resindex)
                if P.st is not None:
                    id_co.add(P.st.resindex)

        id_ds = id_ds - id_co

        id_co_plus = set()
        for resindex in id_co:  # TODO: fix quick and dirty
            id_co_plus.add(resindex)
            for i in range(-plus, plus):
                id_co_plus.add(resindex + i)

        id_all = id_ss | id_ds
        id_clean = id_all - (id_co | id_co_plus | id_ss)
        id_nick = [*self.link.Fnicks.values(),
                   *self.link.Fnicks.keys()]

        return {"co": id_co, "co_plus": id_co_plus, "ss": id_ss, "ds": id_ds,
                "all": id_all, "clean": id_clean, "nick": id_nick}

    def create_df(self, frame=0):
        data, ts = self._traj_frame(frame)
        max_res = max(self.categories["all"])
        id_prop_dict = {}
        for resindex in range(0, max_res + 1):
            position = self.link.DidDhps[self.link.FidDid[resindex]]
            categories = [
                cat for cat in CATEGORIES if resindex in self.categories[cat]]
            strand_type = "scaffold" if position[2] else "staple"
            categories.append(strand_type)
            id_prop_dict[resindex] = [tuple(categories), position]

            sequence = self.link.FidSeq_local[resindex]
            for i, _ in enumerate(sequence):
                if "seq" + str(i) + "-loc" not in self.columns:
                    self.columns.append("seq" + str(i) + "-loc")
                id_prop_dict[resindex].append(sequence[:(i + 1)])
            if "sequence-loc" not in self.columns:
                self.columns.append("sequence-loc")
            id_prop_dict[resindex].append(sequence)

            sequence = self.link.FidSeq_global[resindex]
            for i, _ in enumerate(sequence):
                if "seq" + str(i) + "-glo" not in self.columns:
                    self.columns.append("seq" + str(i) + "-glo")
                id_prop_dict[resindex].append(sequence[:(i + 1)])
            if "sequence-glo" not in self.columns:
                self.columns.append("sequence-glo")
            id_prop_dict[resindex].append(sequence)

            num_HN = self.link.FidHN[resindex]
            for i, _ in enumerate(num_HN):
                if "num_HN" + str(i) not in self.columns:
                    self.columns.append("num_HN" + str(i))
                id_prop_dict[resindex].append(str(num_HN[i]))

            if "strand" not in self.columns:
                self.columns.append("strand")
            id_prop_dict[resindex].append(strand_type)

            if "motif" not in self.columns:
                self.columns.append("motif")
            if resindex in self.categories["co"]:
                motif = "co"
            elif resindex in self.categories["nick"]:
                motif = "nick"
            elif resindex in self.categories["ds"]:
                motif = "ds"
            else:
                motif = "ss"
            id_prop_dict[resindex].append(motif)

            for prop in PROP_TYPE:
                if prop == "bp_quality":
                    try:  # not all are basepaired
                        bp_qual = mean(data[prop][resindex].values())
                    except KeyError:
                        bp_qual = np.nan
                    id_prop_dict[resindex].append(bp_qual)
                    if prop not in self.columns:
                        self.columns.append(prop)
                elif prop == "dh_quality":
                    for dh in DIHEDRALS:
                        angle = data[prop][resindex][dh]
                        id_prop_dict[resindex].append(angle)
                        if ("dh-" + dh) not in self.columns:
                            self.columns.append("dh-" + dh)
                elif prop in ["bp_geometry_local", "bp_geometry_global"]:
                    for geom in WCGEOMETRY:
                        try:
                            cent_ank = data[prop][resindex][geom]
                        except KeyError:
                            cent_ank = {"C1'": np.nan,
                                        "diazine": np.nan, "C6C8": np.nan}
                        for loc, value in cent_ank.items():
                            id_prop_dict[resindex].append(value)
                            col_name = "{}-{}-{}".format(geom,
                                                         loc,
                                                         prop[12:15]
                                                         )
                            if (col_name) not in self.columns:
                                self.columns.append(col_name)
                elif prop == "distances":
                    for atom in ["C1'", "P"]:
                        dist_dicts = data[prop][resindex][atom]
                        for loc, dist in dist_dicts.items():
                            id_prop_dict[resindex].append(dist)
                            if (atom + "-" + loc) not in self.columns:
                                self.columns.append(atom + "-" + loc)
                elif prop == "localres" and data[prop] is not None:
                    try:  # not all are basepaired
                        localres = data[prop][resindex]
                    except KeyError:
                        localres = np.nan
                    id_prop_dict[resindex].append(localres)
                    if prop not in self.columns:
                        self.columns.append(prop)
        self.df = pd.DataFrame.from_dict(
            id_prop_dict, orient='index', columns=self.columns)

        id_co_dict = {}
        for co_id in self.link.Fco.keys():
            try:
                co_type = data["co_angles"][co_id]["type"]
            except KeyError:
                continue  # running index not continuos
            strand_type = ("scaffold" if
                           data["co_angles"][co_id]["is_scaffold"]
                           else "staple")
            id_co_dict[co_id] = [[co_type, strand_type]]

            co_resids = data["co_angles"][co_id]["resindices"]
            if data["localres"] is not None:
                locres = 0.
                for resid in co_resids:
                    locres += (data["localres"][resid] / len(co_resids))

                id_co_dict[co_id].append(locres)

            for name in COANGLES:
                try:
                    angle = data["co_angles"][co_id]["angles"][name]
                except KeyError:
                    angle = None
                id_co_dict[co_id].append(angle)

        if data["localres"] is not None:
            self.df_co = pd.DataFrame.from_dict(
                id_co_dict, orient='index',
                columns=(["type", "localres"] + COANGLES)
            )
        else:
            self.df_co = pd.DataFrame.from_dict(
                id_co_dict, orient='index',
                columns=(["type"] + COANGLES)
            )
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
    _, _, _ = prep.create_df()


if __name__ == "__main__":
    main()
