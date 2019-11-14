#!/usr/bin/env python
import mrcfile as mrc
import numpy as np
import MDAnalysis as mda
# import attr

from typing import Dict, Set, Tuple, FrozenSet

from utils import UnexpectedCaseError
from linkage import Linkage
from project import Project
from design import Design


STAR_HEADER = """data_\nloop_\n_rlnMicrographName #1\n_rlnCoordinateX #2
_rlnCoordinateY #3\n_rlnCoordinateZ #4"""


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
        path_star = "Tomograms/seg-co/" + path_out.stem + ".star"
        path_out_star = path_out.parent / (path_out.stem + ".star")
        with open(path_out_star, mode="w") as star_out:
            star_out.write(STAR_HEADER)
            star_out.write("{} {} {} {}".format(path_star, *center_small))
    return


def categorise(link: Linkage,
               project: Project,
               ) -> Dict[str, Set[Tuple[FrozenSet[int], str, str]]]:

    def _expand_selection(selection: Set[int],
                          link: Linkage,
                          plus: int,
                          ) -> FrozenSet[int]:
        expand = set()
        for resindex in selection:
            h, p, is_scaf = link.DidDhps[link.FidDid[resindex]]
            for i in range(-plus, plus):
                position = (h, p + i, is_scaf)
                expand.add(link.DidFid[link.DhpsDid[position]])
        return frozenset(expand)

    categories = dict()

    plus = project.range
    co_segment = set()
    for key, co in link.Fco.items():
        co_res = set()
        for bp in co.Ps:
            if bp is None:
                continue
            if bp.sc is not None:
                co_res.add(bp.sc.resindex)
            if bp.st is not None:
                co_res.add(bp.st.resindex)

        co_res_plus = _expand_selection(selection=co_res, link=link, plus=plus)
        typ = "co-{}".format(co.typ)
        idf = key.strip("[]").replace(" ", "").replace(",", "-")
        identifier = idf.replace(")-(", "_").strip("()")
        co_segment.add(tuple([co_res_plus, identifier, typ]))
    categories["co"] = co_segment

    nick_segment = set()
    nick_done: Set[int] = set()
    for res, ser in iter(link.Fnicks.items()):
        if res in nick_done:
            continue
        nick_done.update([res, ser])

        res_bp = link.Fbp_full[res]
        ser_bp = link.Fbp_full[ser]
        nick = set([res, ser, res_bp, ser_bp])

        nick_plus = _expand_selection(selection=nick, link=link, plus=plus)
        h, p, _ = link.DidDhps[link.FidDid[res]]
        idenifier = "{}-{}".format(h, p)
        typ = "nick"
        nick_segment.add(tuple([nick_plus, idenifier, typ]))
    categories["nick"] = nick_segment

    design = Design(project)
    ds_domain = set()
    for domain in design.design.domain_list:
        bases = domain.base_list
        across_id = domain.connected_domain
        is_long_ds_staple = (len(bases) >= 14
                             and across_id != -1
                             and not bases[0].is_scaf
                             )
        if not is_long_ds_staple:
            continue
        else:
            across = design.design.domain_list[across_id]
            bases += across.base_list
            domain_resindices = frozenset({link.DidFid[base.id]
                                           for base in bases
                                           if base.id in link.DidFid  # !skip
                                           }
                                          )
            idenifier = str(domain.id)
            typ = "ds_domain"
            ds_domain.add(tuple([domain_resindices, idenifier, typ]))
    categories["ds"] = ds_domain

    return categories
