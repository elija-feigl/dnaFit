from typing import Dict, Set, Tuple, FrozenSet, Any
from pathlib import Path

from ..link.linkage_pickle import Linkage
from ..data.design import Design

""" collection of scripts to allow creating subsets of a cryo-EM map.

    COMMENTS:
    code not well maintained
"""


STAR_HEADER = """data_\n\nloop_\n_rlnMicrographName #1\n_rlnCoordinateX #2
_rlnCoordinateY #3\n_rlnCoordinateZ #4\n"""


def categorise(link: Linkage,
               json: Path,
               seq: Path,
               plus: int,
               ) -> Dict[str, Set[Tuple[Any, ...]]]:

    def _expand_selection(selection: Set[int],
                          link: Linkage,
                          plus: int,
                          ) -> FrozenSet[int]:
        expand = set()
        for resindex in selection:
            h, p, is_scaf = link.DidDhps[link.FidDid[resindex]]
            # NOTE: insertions have negative position values
            p = abs(p)  # hotfix for single insertions
            for i in range(-plus, plus):
                position = (h, p + i, is_scaf)
                if position in link.DhpsDid:  # skips included as None
                    expand.add(link.DidFid[link.DhpsDid[position]])
        return frozenset(expand)

    categories = dict()

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
        typ = f"co-{co.typ}"
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
        # NOTE: insertions have negative position values
        p = abs(p)  # hotfix for single insertions
        identifier = f"{h}-{p}"
        typ = "nick"
        nick_segment.add(tuple([nick_plus, identifier, typ]))
    categories["nick"] = nick_segment

    design = Design(json=json, seq=seq,)
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
            identifier = str(domain.id)
            typ = "ds_domain"
            ds_domain.add(tuple([domain_resindices, identifier, typ]))
    categories["ds"] = ds_domain

    return categories
