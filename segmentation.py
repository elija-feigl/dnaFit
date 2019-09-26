#!/usr/bin/env python
import mrcfile as mrc
import numpy as np
import MDAnalysis as mda
import pickle
import argparse
import attr
import os

from pathlib import Path
from typing import List, Set, Dict, Tuple


@attr.s(slots=True)
class Project(object):
    input: Path = attr.ib()
    output: Path = attr.ib()
    name: str = attr.ib()
    context: int = attr.ib()
    range: int = attr.ib()
    halfmap: bool = attr.ib()
    localres: bool = attr.ib()
    star: bool = attr.ib()


def mrc_segment(atoms: "mda.atomgroup",
                path_in: str,
                path_out: str,
                context: int=3,
                star: bool= False,
                ) -> None:

    def remove_padding(data: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        idx_data = np.nonzero(data)

        pos_min = np.min(idx_data, axis=1)
        pos_max = np.max(idx_data, axis=1)
        s_cell = np.max(pos_max-pos_min)
        pad = 0.5 * (s_cell + pos_min - pos_max)
        pos_low = [int(p) if (p % 1.) == 0. else int(p) + 1 for p in pad]
        pos_high = -pad.astype(int)

        data_small = np.zeros(np.full(3, s_cell), dtype=np.float32)
        data_small[pos_low[0]: pos_high[0] or None,
                   pos_low[1]: pos_high[1] or None,
                   pos_low[2]: pos_high[2] or None
                   ] = data[pos_min[0]: pos_max[0],
                            pos_min[1]: pos_max[1],
                            pos_min[2]: pos_max[2]
                            ]
        v_origin_small = pos_min - pos_low
        return data_small, v_origin_small

    if not len(atoms):
        print("EXIT - no atoms in this selection")
        return
    u = atoms.universe
    u.trajectory[-1]

    with mrc.open(path_in, mode='r') as mrc_in:
        o = np.array(mrc_in.header["origin"])
        origin = np.array([o["x"], o["y"], o["z"]])
        c = np.array(mrc_in.header["cella"])
        cellA = np.array([c["x"], c["y"], c["z"]])
        shape = np.array([mrc_in.header["nx"],
                          mrc_in.header["ny"],
                          mrc_in.header["nz"]
                          ])
        voxel_size = cellA / shape
        v_context = np.full(3, context / voxel_size).astype(int) + 1
        data_all = np.swapaxes(mrc_in.data, 0, 2)  # TODO: -low- faster wtht swap

    data_mask = np.zeros(shape, dtype=np.float32)
    atoms_voxel = np.rint((atoms.positions - origin) / voxel_size)
    for voxel in atoms_voxel.astype(int):
        low = voxel - v_context
        high = voxel + v_context
        data_mask[low[0]:high[0], low[1]:high[1], low[2]:high[2]] = 1.

    data = data_all * data_mask
    data_small, v_origin_small = remove_padding(data=data)

    shape_small = np.shape(data_small)
    origin_small = origin + (v_origin_small * voxel_size)
    center_small = np.divide(shape_small, 2).astype(int)

    with mrc.new(path_out, overwrite=True) as mrc_out:
        mrc_out.set_data(np.swapaxes(data_small, 0, 2))  # TODO
        mrc_out._set_voxel_size(*(voxel_size))
        mrc_out.header["origin"] = tuple(origin_small)

    if star:
        path_out_split = path_out.split("/")
        path_star = "Tomograms/seg-co/" + path_out_split[-1]
        star_header = """data_\n
                        loop_
                        _rlnMicrographName #1
                        _rlnCoordinateX #2
                        _rlnCoordinateY #3
                        _rlnCoordinateZ #4
                        """.replace(" ", "")
        with open(path_out + ".star", mode="w") as star_out:
            star_out.write(star_header)
            star_out.write("{}.star {} {} {}".format(path_star, *center_small))
    return


def _mrc_localres(atoms: "mda.atomgroup", path_in: str) -> dict:
    def get_localres(atoms: "mda.atomgroup",
                     data: np.ndarray,
                     origin: np.ndarray,
                     voxel_size: np.ndarray,
                     ) -> None:
        locres = 0.
        for atom in atoms:
            atom_voxel = np.rint((atom.position - origin) / voxel_size)
            locres += data[tuple(atom_voxel.astype(int))]
        locres /= len(atoms)
        return locres

    if not len(atoms):
        print("EXIT - no atoms in this selection")
        return

    u = atoms.universe
    u.trajectory[-1]

    with mrc.open(path_in, mode='r') as mrc_in:
        o = np.array(mrc_in.header["origin"])
        origin = np.array([o["x"], o["y"], o["z"]])
        c = np.array(mrc_in.header["cella"])
        cellA = np.array([c["x"], c["y"], c["z"]])
        shape = np.array([mrc_in.header["nx"],
                          mrc_in.header["ny"],
                          mrc_in.header["nz"]
                          ])
        voxel_size = cellA / shape
        data = np.swapaxes(mrc_in.data, 0, 2)

    dict_localres = {}
    for res in atoms.residues:
        localres = get_localres(res.atoms, data, origin, voxel_size)
        dict_localres[res.resindex] = localres

    return dict_localres


def _categorise_lists(topo, plus=3):
    # TODO. check: names set list etc
    dict_FidDid = {v: k for k, v in iter(topo["dict_idid"].items())}
    dict_DidDhps = {v: k for k, v in iter(topo["dict_hpid"].items())}
    dict_bpFULL = {**topo["dict_bp"], **{v: k for k, v in iter(topo["dict_bp"].items())}}

    id_ds = set()
    id_coplus = set()

    for wc_id1, wc_id2 in iter(topo["dict_bp"].items()):
        id_ds.add(wc_id1)
        id_ds.add(wc_id2)
    id_ss = set(topo["dict_idid"].values()) - id_ds

    id_co = set()
    id_co_init = {id_design for id_design in topo["dict_coid"].keys()
                  if id_design not in id_ss}
    allready_done = set()
    for base in id_co_init:
        typ = topo["dict_coid"][base]["type"][0]
        co_index = topo["dict_coid"][base]["co_index"]

        if base not in allready_done:
            allready_done.add(base)
            co = topo["dict_coid"][base]["co"]
            allready_done.add(co)

            co_bp = dict_bpFULL[co]
            bp = dict_bpFULL[base]

            if topo["dict_coid"][base]["type"][0] == "double":
                dou = topo["dict_coid"][base]["type"][1]
                allready_done.add(dou)
                dou_co = topo["dict_coid"][dou]["co"]
                allready_done.add(dou_co)

                dou_co_bp = dict_bpFULL[dou_co]
                dou_bp = dict_bpFULL[dou]

                tup = (base, bp, co, co_bp, dou,
                       dou_bp, dou_co, dou_co_bp, co_index, typ)
            else:
                tup = (base, bp, co, co_bp, co_index, typ)

            tup_plus = []
            for x in tup[:-2]:
                h, p, is_scaf = dict_DidDhps[dict_FidDid[x]]
                for i in range(-plus, plus):
                    try:
                        tup_plus.append(
                            topo["dict_idid"][topo["dict_hpid"][(h,
                                                                 p+i,
                                                                 is_scaf)]])
                    except KeyError:
                        pass  # helix end

            tup_plus.append(co_index)
            tup_plus.append(typ)
            id_co.add(tup)
            id_coplus.add(tuple(tup_plus))

    nick_allready_done = set()
    id_nick = set()
    for id1, id2 in iter(topo["dict_nicks"].items()):
        if id1 not in nick_allready_done:
            nick_allready_done.add(id1)
            nick_allready_done.add(id2)
            tup = (id1, id2, dict_bpFULL[id1], dict_bpFULL[id2])
            id_nick.add(tup)

    id_nick_plus = []
    for nick in id_nick:
        tup_plus = []
        for x in nick:
            h, p, is_scaf = dict_DidDhps[dict_FidDid[x]]
            for i in range(-plus, plus):
                try:
                    tup_plus.append(
                        topo["dict_idid"][topo["dict_hpid"][(h,
                                                             p+i,
                                                             is_scaf)]])
                except KeyError:
                    pass  # helix end
        id_nick_plus.append(tup_plus)
    return id_co, id_coplus, id_nick, id_nick_plus


def _topology(project: Project) -> dict:
    DICTS = ["dict_bp", "dict_idid", "dict_hpid", "dict_color",
             "dict_coid", "dict_nicks", "list_skips", "universe"]
    topo = {}
    for dic in DICTS:
        pickle_file = project.output / "{}__{}.p".format(project.name, dic)
        topo[dic] = pickle.load(open(pickle_file, "rb"))

    return topo


def get_description() -> str:
    return """cut subset from map according to atoms belongign to sepcific
              motif. Also produces minimal box map. can also segment halfmaps
              and evaluate local-resolution per residue -> dict and pdb
              """


def proc_input() -> Project:
    parser = argparse.ArgumentParser(
        description=get_description(),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument("--folder",
                        help="input folder",
                        type=str,
                        default="./",
                        )
    parser.add_argument("--name",
                        help="name of design and files",
                        type=str,
                        required=True,
                        default=argparse.SUPPRESS,
                        )
    parser.add_argument("--context",
                        help="context aroud atom in Angstrom",
                        type=int,
                        default=5,
                        )
    parser.add_argument("--range",
                        help="number of addit. basepairs in helix of motif",
                        type=int,
                        default=10,
                        )
    parser.add_argument("--halfmap",
                        help="also segment halfmaps",
                        action="store_true"
                        )
    parser.add_argument("--localres",
                        help="compute localres per molecule",
                        action="store_true"
                        )
    parser.add_argument("--star",
                        help="create starfile",
                        action="store_true"
                        )
    args = parser.parse_args()
    project = Project(input=Path(args.folder),
                      output=Path(args.folder) / "analysis",
                      name=args.name,
                      context=args.context,
                      range=args.range,
                      halfmap=args.halfmap,
                      localres=args.localres,
                      star=args.star,
                      )
    return project


def mask_minimal_box(u, project):
    path_in = project.input / "{}.mrc".format(project.name)
    path_out = project.output / "{}-masked.mrc".format(project.name)
    mrc_segment(atoms=u.atoms,
                path_in=path_in,
                path_out=path_out,
                context=project.context,
                )
    return


def local_res(u, path_color, project):
    dict_localres = _mrc_localres(atoms=u.atoms,
                                  path_in=path_color,
                                  )
    path_colorpickle = project.output / "{}__localres.p".format(project.name)
    pickle.dump(dict_localres, open(path_colorpickle, "wb"))
    path_colorpdb = project.output / "{}_localres.pdb".format(project.name)
    pdb = mda.Writer(path_colorpdb, multiframe=True)
    u.add_TopologyAttr(mda.core.topologyattrs.Tempfactors(np.zeros(len(u.atoms))))
    u.atoms.tempfactors = -1.
    for res in u.residues:
            res.atoms.tempfactors = dict_localres[res.resindex]
    pdb.write(u.atoms)
    return


def main():
    project = proc_input()
    print("input from ", project.input)
    topo = _topology(project=project)
    _, id_coplus_lists, _, id_nickplus_list = _categorise_lists(
        topo,
        plus=project.range,
    )
    # initialize universe and select final frame
    u = mda.Universe(*topo["universe"])
    u.trajectory[-1]

    print("mask minimal box")
    mask_minimal_box(u, project)

    path_color = project.input / "{}_localres.mrc".format(project.name)
    if os.path.isfile(path_color):
        print("compute per residue resolution")
        local_res(u, path_color, project)

    motif_cat = {"co": id_coplus_lists, "nick": id_nickplus_list}
    for motif in ["nick", "co"]:
        path_motif = project.output / motif
        print("output to ", path_motif)
        try:
            os.mkdir(path_motif)
        except FileExistsError:
            pass

        for index, co_select_typ in enumerate(motif_cat[motif]):
            if motif == "co":
                co_select = co_select_typ[:-2]
                typ = co_select_typ[-1]
                index = co_select_typ[-2]
                atoms_select = mda.AtomGroup([], u)
                for base_id in co_select:
                    atoms_select += u.residues[base_id].atoms

            elif motif == "nick":
                typ = ""
                atoms_select = mda.AtomGroup([], u)
                for base_id in co_select_typ:
                    atoms_select += u.residues[base_id].atoms

            if project.halfmap:
                print("segmenting halfmaps")
                specs = {"": "", "_unfil_half1": "h1-", "_unfil_half2": "h2-"}
            else:
                specs = {"": ""}
            for inp, out in specs.items():
                path_in = project.input / "{}{}.mrc".format(project.name,
                                                            inp,
                                                            )
                path_out = path_motif / "{}__{}{}{}{}.mrc".format(project.name,
                                                                  out,
                                                                  typ,
                                                                  motif,
                                                                  index,
                                                                  )
                mrc_segment(atoms=atoms_select,
                            path_in=path_in,
                            path_out=path_out,
                            context=project.context,
                            star=project.star,
                            )
    return


if __name__ == "__main__":
    main()
