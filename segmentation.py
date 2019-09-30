#!/usr/bin/env python
import mrcfile as mrc
import numpy as np
import MDAnalysis as mda
import pickle
import argparse
import attr
import os
import contextlib

from pathlib import Path
from typing import List, Set, Dict, Tuple, Any, FrozenSet
from itertools import chain


class UnexpectedCaseError(Exception):
    """Raised when a case occurs that makes no sense in the programs context"""
    pass


@contextlib.contextmanager
def ignored(*exceptions):
    try:
        yield
    except exceptions:
        pass


@attr.s(slots=True)
class Project(object):
    input: Path = attr.ib()
    output: Path = attr.ib()
    name: str = attr.ib()
    # specific
    context: int = attr.ib()
    range: int = attr.ib()
    halfmap: bool = attr.ib()
    localres: bool = attr.ib()
    star: bool = attr.ib()


@attr.s(auto_attribs=True)
class Linkage(object):
    Fbp: Dict[int, int] = {}
    DidFid: Dict[int, int] = {}
    DhpsDid: Dict[Tuple[int, int, bool], int] = {}
    Dcolor: Dict[int, int] = {}
    Dskips: Set[Tuple[int, int, bool]] = set()
    Fnicks: Dict[int, int] = {}
    FidSeq: Dict[int, str] = {}
    Fco: Dict[int, Any] = {}
    universe: Tuple[str, str] = ("", "")

    def dump_linkage(self, project: Project) -> None:
        for name, link in vars(self).items():
            output = project.output / "{}__{}.p".format(project.name, name)
            pickle.dump(link, open(output, "wb"))
        return

    def load_linkage(self, project: Project) -> None:
        for name in vars(self).keys():
            input = project.output / "{}__{}.p".format(project.name, name)
            value = pickle.load(open(input, "rb"))
            setattr(self, name, value)
        return

    def reverse(self) -> None:
        def reverse_d(dict: dict) -> dict:
            return {v: k for k, v in iter(dict.items())}

        self.FidDid = reverse_d(self.DidFid)
        self.DidDhps = reverse_d(self.DhpsDid)
        self.Fbp_rev = reverse_d(self.Fbp)
        self.Fbp_full = {**self.Fbp, **self.Fbp_rev}
        return


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
        raise UnexpectedCaseError("no atoms in this selection")

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
        data_all = np.swapaxes(mrc_in.data, 0, 2)  # TODO: -low- faster wo swap

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


def _mrc_localres(atoms: "mda.atomgroup", path_in: str) -> Dict[int, float]:
    def get_localres(atoms: "mda.atomgroup",
                     data: np.ndarray,
                     origin: np.ndarray,
                     voxel_size: np.ndarray,
                     ) -> float:
        locres = 0.
        for atom in atoms:
            atom_voxel = np.rint((atom.position - origin) / voxel_size)
            locres += data[tuple(atom_voxel.astype(int))]
        locres /= len(atoms)
        return locres

    if not len(atoms):
        raise UnexpectedCaseError("no atoms in this selection")

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


def categorise(link: Linkage,
               plus: int=3
               ) -> Tuple[Set[Tuple[FrozenSet[int], int, str]],
                          Set[FrozenSet[int]]
                          ]:

    def _get_co(res: int, link: Linkage) -> List[int]:
        ser = link.Fco[res]["co"]
        res_bp = link.Fbp_full[res]
        ser_bp = link.Fbp_full[ser]
        return [res, ser, res_bp, ser_bp]

    def _expand_co(co: Tuple[List[int], int, str],
                   link: Linkage,
                   ) -> Tuple[FrozenSet[int], int, str]:
        bases = co[0][:len(co)]
        co_plus = _expand_selection(selection=bases, link=link)
        return (frozenset(co_plus), co[1], co[2])

    def _expand_nick(n: Tuple[int, int, int, int],
                     link: Linkage
                     ) -> FrozenSet[int]:
        n_plus = _expand_selection(selection=list(n), link=link)
        return frozenset(n_plus)

    def _expand_selection(selection: List[int], link: Linkage) -> Set[int]:
        expand = set()
        for resindex in selection:
            h, p, is_scaf = link.DidDhps[link.FidDid[resindex]]
            for i in range(-plus, plus):
                position = (h, p+i, is_scaf)
                with ignored(KeyError):
                    expand.add(link.DidFid[link.DhpsDid[position]])
        return expand

    link.reverse()
    ds = set(chain.from_iterable((a, b) for a, b in iter(link.Fbp.items())))
    ss = set(link.DidFid.values()) - ds

    co = set()
    co_plus = set()
    co_done: Set[int] = set()
    co_init = set(resid for resid in link.Fco.keys() if resid not in ss)
    for res in co_init:
        if res in co_done:
            continue

        typ = link.Fco[res]["type"][0]
        co_index = link.Fco[res]["co_index"]

        half = _get_co(res=res, link=link)

        if typ == "double":
            f_res = link.Fco[res]["type"][1]
            flah = _get_co(res=f_res, link=link)
            full = half[:2] + flah[:2] + half[2:] + flah[2:]
            crossover = (frozenset(full), co_index, typ)
            co_done.update(full[:4])
        else:
            crossover = (frozenset(half), co_index, typ)
            co_done.update(half[:2])
        co.add(crossover)
        co_plus.add(_expand_co(co=(half, co_index, typ), link=link))

    nick = set()
    nick_plus = set()
    nick_done: Set[int] = set()
    for res, ser in iter(link.Fnicks.items()):
        if res in nick_done:
            continue
        nick_done.update([res, ser])
        res_bp = link.Fbp_full[res]
        ser_bp = link.Fbp_full[ser]
        n = (res, ser, res_bp, ser_bp)
        nick.add(frozenset(n))
        nick_plus.add(_expand_nick(n=n, link=link))

    return co_plus, nick_plus


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
    empty_TopologyAttr = np.zeros(len(u.atoms))
    u.add_TopologyAttr(mda.core.topologyattrs.Tempfactors(empty_TopologyAttr))
    u.atoms.tempfactors = -1.
    for res in u.residues:
            res.atoms.tempfactors = dict_localres[res.resindex]
    pdb.write(u.atoms)
    return


def main():
    project = proc_input()

    print("input from ", project.input)
    linkage = Linkage()
    linkage.load_linkage(project=project)
    co, nick = categorise(link=linkage, plus=project.range)
    u = mda.Universe(*linkage.universe)
    u.trajectory[-1]

    print("mask minimal box")
    mask_minimal_box(u, project)

    if project.localres:
        print("compute per residue resolution")
        path_color = project.input / "{}_localres.mrc".format(project.name)
        local_res(u, path_color, project)

    motifs = {"co": co, "nick": nick}
    if project.halfmap:
        print("segmenting halfmaps")
    for motif_name, motif in motifs.items():
        path_motif = project.output / motif_name
        print("output to ", path_motif)
        with ignored(FileExistsError):
            os.mkdir(path_motif)

        for index, subset in enumerate(motif):
            if motif_name == "co":
                base_selection, index, typ = subset # TOdO: co-index
                atoms_select = mda.AtomGroup([], u)
                for resindex in base_selection:
                    atoms_select += u.residues[resindex].atoms
            elif motif_name == "nick":
                typ = ""
                atoms_select = mda.AtomGroup([], u)
                for base_id in subset:
                    atoms_select += u.residues[base_id].atoms

            if project.halfmap:
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
                                                                  motif_name,
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
