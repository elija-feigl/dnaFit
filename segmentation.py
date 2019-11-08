#!/usr/bin/env python
import mrcfile as mrc
import numpy as np
import MDAnalysis as mda
# import attr

from typing import List, Set, Tuple, FrozenSet
from itertools import chain

from utils import UnexpectedCaseError, ignored
from linkage import Linkage


def mrc_segment(atoms: "mda.atomgroup",
                path_in: str,
                path_out: str,
                context: int = 3,
                star: bool = False,
                ) -> None:

    def remove_padding(data: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        idx_data = np.nonzero(data)

        pos_min = np.min(idx_data, axis=1)
        pos_max = np.max(idx_data, axis=1)
        s_cell = np.max(pos_max - pos_min)
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


def categorise(link: Linkage,
               plus: int = 3
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
                position = (h, p + i, is_scaf)
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
