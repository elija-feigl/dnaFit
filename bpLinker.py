#!/usr/bin/env python3#
import os

import MDAnalysis as mda  # type:ignore
import numpy as np
import pickle
import contextlib
import argparse
import attr
import nanodesign as nd

from pathlib import Path
from itertools import chain
from operator import attrgetter
from typing import List, Set, Dict, Tuple, Any, TextIO
from nanodesign.converters import Converter


@contextlib.contextmanager
def ignored(*exceptions):
    try:
        yield
    except exceptions:
        pass


class UnexpectedCaseError(Exception):
    """Raised when a case occurs that makes no sense in the programs context"""
    pass


@attr.s(slots=True)
class Project(object):
    input: Path = attr.ib()
    output: Path = attr.ib()
    name: str = attr.ib()
    # specific
    ENmodify: bool = attr.ib()
    EN: str = attr.ib()


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


@attr.s
class Linker(object):
    """ Linker class
    """
    # TODO: move categorize to linker?
    project: Project = attr.ib()
    Fbp: Dict[int, int] = {}
    DidFid: Dict[int, int] = {}
    DhpsDid: Dict[Tuple[int, int, bool], int] = {}
    Fnicks: Dict[int, int] = {}
    FidSeq: Dict[int, str] = {}
    Dskips: Set[Tuple[int, int, bool]] = set()
    Fco: Dict[int, Any] = {}

    def __attrs_post_init__(self):
        self.fit: Fit = Fit(self.project)
        self.design: Design = Design(self.project)
        self.Dskips = self._eval_skips()
        self.FidSeq = self._eval_sequence()

    def _eval_sequence(self) -> Dict[int, str]:
        FidSeq = {res.resindex: res.resname[0]
                  for res in self.fit.u.residues
                  }
        return FidSeq

    def _eval_skips(self) -> Set[Tuple[int, int, bool]]:
        """ Returns
            -------
                (base.h, base.p, base.is_scaf) of all skips
        """
        design_allbases = [base
                           for strand in self.design.strands
                           for base in strand.tour
                           ]
        Dskips = set((base.h, base.p, base.is_scaf)
                     for base in design_allbases
                     if self._is_del(base)
                     )
        return Dskips

    def _is_del(self, base: "nd.base") -> bool:
        return base.num_deletions != 0

    def create_linkage(self) -> Linkage:
        """ invoke _link_scaffold, _link_staples, _link_bp to compute mapping
            of every base design-id to fit-id as well as the basepair mapping.
            basepairs are mapped from scaffold to staple, unique (invertable).
            updates linker attributes corresponding to the respective mapping
            and returns them.
        -------
         Returns
            -------
                Linkage
        """

        self._link()
        self._identify_bp()
        self._identify_crossover()
        self._identify_nicks()
        universe = self._get_universe_tuple()

        return Linkage(Fbp=self.Fbp,
                       DidFid=self.DidFid,
                       DhpsDid=self.DhpsDid,
                       Dcolor=self.Dcolor,
                       Dskips=self.Dskips,
                       Fco=self.Fco,
                       Fnicks=self.Fnicks,
                       FidSeq=self.FidSeq,
                       universe=universe,
                       )

    def _link(self) -> Tuple[Dict[int, int],
                             Dict[Tuple[int, int, bool], int],
                             ]:
        def link_scaffold() -> Tuple[Dict[int, int],
                                     Dict[Tuple[int, int, bool], int],
                                     ]:
            """ collect position in scaffold (0-x) by comparing index in list
                of scaffold_design positions
            -------
                Returns
                -------
                DidFid
                    design-id -> fit-id
                DhpsDid
                    helix-number, base-position, is_scaffold -> design-id
            """
            Dscaffold = self.design.scaffold
            Did = [base.id for base in Dscaffold]
            Dhp = [(base.h, base.p, True) for base in Dscaffold]
            Fid_local = [Did.index(base.id) for base in Dscaffold]
            Fid_global = self.fit.scaffold.residues[Fid_local].resindices

            DidFid = dict(zip(Did, Fid_global))
            DhpsDid = dict(zip(Dhp, Did))
            return (DidFid, DhpsDid)

        def link_staples() -> Tuple[Dict[int, int],
                                    Dict[Tuple[int, int, bool], int],
                                    Dict[int, int],
                                    ]:
            """same procedure as scaffold for each
            -------
            Returns
                -------
                DidFid
                    design-id -> fit-id
                DhpsDid
                    helix-number, base-position, is_scaffold -> design-id
                dict color
                    fit-segment-id -> color
            """
            def get_resid(segindex: int, resindex_local: int) -> int:
                segment = self.fit.staples[segindex]
                return segment.residues[resindex_local].resindex

            DidFid: Dict[int, int] = {}
            DhpsDid: Dict[Tuple[int, int, bool], int] = {}
            color: Dict[int, int] = {}

            for i, staple in enumerate(self.design.staples):
                seg_id = self.design.stapleorder[i]

                Did = [base.id for base in staple]
                Dhp = [(base.h, base.p, False) for base in staple]

                Fid_local = [Did.index(base.id)for base in staple]
                Fid_global = [get_resid(seg_id, resid) for resid in Fid_local]

                icolor = self.design.design.strands[staple[0].strand].icolor
                segidxforcolor = self.fit.staples[seg_id].segindex
                color[segidxforcolor] = icolor

                DidFid_add = dict(zip(Did, Fid_global))
                DhpsDid_add = dict(zip(Dhp, Did))
                DidFid = {**DidFid, **DidFid_add}
                DhpsDid = {**DhpsDid, **DhpsDid_add}

            return (DidFid, DhpsDid, color)

        DidFid_sc, DhpsDid_sc = link_scaffold()
        DidFid_st, DhpsDid_st, self.Dcolor = link_staples()

        self.DidFid = {**DidFid_sc, **DidFid_st}
        self.DhpsDid = {**DhpsDid_sc, **DhpsDid_st}
        return (self.DidFid, self.DhpsDid)

    def _identify_bp(self) -> Dict[int, int]:
        """ link basepairs by mapping indices according to json (cadnano).
            basepairs are mapped from scaffold to staple, unique (invertable).
        -------
         Returns
            -------
            self.Fbp
                fit-id -> fit-id
        """
        def Fid(Did: int) -> int:
            return self.DidFid[Did]

        self.Fbp = {Fid(base.id): Fid(base.across.id)
                    for base in self.design.scaffold
                    if base.across is not None
                    }
        return self.Fbp

    def _get_nextInHelix(self, h: int, p: int, is_scaf: bool, i: int
                         ) -> Tuple[int, int, bool]:
        Dhps = (h, p + i, is_scaf)
        while Dhps in self.Dskips:
            i += np.sign(i)
            Dhps = (h, p + i, is_scaf)
        return Dhps

    def _identify_crossover(self) -> Dict[int, Any]:
        """ for every base id that is involved in a crossover
            updates linker attribute of crossovers and returns it
        -------
         Returns
            -------
            self.Fco
                fit-id ->
                "co_index": co_index, "co": co_id,
                "leg": leg_id,
                "is_scaffold": is_scaf,
                "position": (h, p, is_scaf)
                if end, single, dpouble
                    "type": ("end", None)
                    "type": ("single", single, single_leg)
                    "type": ("double", double)
        """
        # TODO: double, single -> full, half
        def add_co_type() -> None:

            def get_co_leg_p(n_Dhps: Tuple[int, int, bool], i: int) -> int:
                """determine leg base by Dhps and +1/-1 (def: 2 bases away)"""
                h, p, is_scaf = n_Dhps
                Dhps = self._get_nextInHelix(h, p, is_scaf, i)
                leg_Did = self.DhpsDid[Dhps]
                return self.DidFid[leg_Did]

            def is_nextInStrand(b_Dhps: Tuple[int, int, bool],
                                a_Dhps: Tuple[int, int, bool]
                                ) -> bool:
                a_Fid = self.DidFid[self.DhpsDid[a_Dhps]]
                b_Fid = self.DidFid[self.DhpsDid[b_Dhps]]
                return True if abs(a_Fid - b_Fid) == 1 else False

            for co in iter(self.Fco.values()):
                h, p, is_scaf = co["position"]
                for i in [-1, 1]:
                    n_Dhps = self._get_nextInHelix(h, p, is_scaf, i)
                    n_Did = self.DhpsDid.get(n_Dhps, None)
                    n_Fid = self.DidFid.get(n_Did, None)

                    if n_Fid is None:
                        co["type"] = ("end", None, None)
                    elif is_nextInStrand(co["position"], n_Dhps):
                        continue
                    else:
                        is_double = (n_Fid in iter(self.Fco.keys()))
                        if is_double:
                            co["type"] = ("double", n_Fid, None)
                        else:
                            leg = get_co_leg_p(n_Dhps, i)
                            co["type"] = ("single", n_Fid, leg)
                co.pop("base")
            return

        def get_co_leg_id(base: "nd.base", direct: str) -> int:
            """determine leg base base and up/down (def: 2 bases away)"""
            m = base.down.p if direct == "up" else base.up.p
            i = (m - base.p) * 2.
            Dhps = self._get_nextInHelix(base.h, base.p, base.is_scaf, i)
            leg_Did = self.DhpsDid[Dhps]
            return self.DidFid[leg_Did]

        def get_next(base: "nd.base", direct: str) -> "nd.base":
            n = (base.up if direct == "up" else base.down)
            if n is None:
                return None
            else:
                n_position = (n.h, n.p, n.is_scaf)
                while n_position in self.Dskips:
                    n = (n.up if direct == "up" else n.down)
                    if n is None:
                        return None
                    n_position = (n.h, n.p, n.is_scaf)
                return n

        def is_co(base: "nd.base", neighbor: "nd.base", direct: str) -> bool:
            if neighbor is None:
                return False
            else:
                while self._is_del(neighbor):
                    neighbor = get_next(neighbor, direct)
                return neighbor.h != base.h

        def get_co(base: "nd.base",
                   neighbor: "nd.base",
                   direct: str,
                   run_id: int
                   ) -> Tuple[dict, int]:
            co_id = self.DidFid[neighbor.id]
            leg_id = get_co_leg_id(base, direct)
            Dhps = (base.h, base.p, base.is_scaf)

            try:  # TODO -high cleanup co_index
                index = self.Fco[co_id]["co_index"]
            except KeyError:
                run_id += 1
                index = run_id
            data = {"co_index": index,
                    "co": co_id,
                    "leg": leg_id,
                    "is_scaffold": base.is_scaf,
                    "position": Dhps,
                    "base": base,
                    }
            return data, run_id

        allbases = self.design.scaffold.copy()
        staples = [base for staple in self.design.staples for base in staple]
        allbases.extend(staples)

        run_id = 0
        for base in allbases:
            base_Fid = self.DidFid[base.id]
            for direct in ["up", "down"]:
                neighbor = get_next(base, direct)
                if is_co(base, neighbor, direct):
                    co_data, run_id = get_co(base, neighbor, direct, run_id)
                    self.Fco[base_Fid] = co_data
                    break
        add_co_type()
        return self.Fco

    def _identify_nicks(self) -> Dict[int, int]:
        """ for every nick, provide id of base accross nick, bidirectional
        -------
         Returns
            -------
            self.Fnicks
                fit-id -> fit_id
        """
        def is_nick(candidate: "nd.base", base: "nd.base") -> bool:
            is_onhelix = (candidate.h == base.h)
            is_neighbor = (abs(base.p - candidate.p) <= 2)  # skip = 2
            is_base = (candidate is base)
            b_Fid = self.DidFid[base.id]
            c_Fid = self.DidFid[candidate.id]
            is_ds_staple = (b_Fid in self.Fbp.values() and
                            c_Fid in self.Fbp.values()
                            )
            return (is_onhelix and
                    is_neighbor and
                    not is_base and
                    is_ds_staple
                    )

        def Fid(Did: int) -> int:
            return self.DidFid[Did]

        staple_end_bases = list(chain.from_iterable((s[0], s[-1])
                                for s in self.design.staples)
                                )
        self.Fnicks = {Fid(base.id): Fid(candidate.id)
                       for base in staple_end_bases
                       for candidate in staple_end_bases
                       if is_nick(candidate, base)
                       }
        return self.Fnicks

    def _get_universe_tuple(self) -> Tuple[str, str]:
        infile = self.project.input / self.project.name
        top = infile.with_suffix(".psf")
        trj = infile.with_suffix(".dcd")
        return (str(top), str(trj))


@attr.s(slots=True, cmp=False, auto_attribs=True)
class ENBond(object):
    """ Elatic Networt Bond class
        harmonic bond of the from k*(r(ab)-r0)**2
        http://www.ks.uiuc.edu/Research/namd/2.7/ug/node26.html
        -------
         Input
            -------
                a1: atomnumber 1
                a2: atomnumber 2
                k: force konstant
                r0: equilibrium distance [A]
                btype: keywords that indicate topology
    """
    a1: int = 0
    a2: int = 0
    k: float = 0.
    r0: float = 0.
    btype: Set[str] = attr.Factory(set)

    def __str__(self) -> str:
        return "bond {b.a1} {b.a2} {b.k} {b.r0}".format(b=self)


@attr.s(slots=True, auto_attribs=True)
class Logic(object):
    """ Logic base class for elastic networks
    """
    long: bool = False
    strand: bool = False
    Hbond: bool = False
    crossstack: bool = False
    nick: bool = False
    co: bool = False
    ssDNA: bool = False
    dihedral: bool = False


@attr.s
class ElaticNetwortModifier(object):
    """ Elatic Networt Modifier class
    """
    linker: Linker = attr.ib()

    def __attrs_post_init__(self):
        self.u: "mda.universe" = self.linker.fit.u
        self.Fbp_full: Dict[int, int] = {**self.linker.Fbp,
                                         **{v: k for k, v
                                            in self.linker.Fbp.items()
                                            }
                                         }
        self.network: set = self._get_network()

    def _get_network(self) -> set:
        infile = self.linker.project.input / self.linker.project.name
        exb_filepath = infile.with_suffix(".exb")
        if exb_filepath.exists():
            with open(exb_filepath) as exb_file:
                network = self._process_exb(exb_file)
        else:
            raise FileNotFoundError
        return network

    def _categorize_bond(self, a1: int, a2: int, r0: float) -> Set[str]:
        def categorize_logic(a1: int, a2: int, r0: float) -> Logic:
            bond_logic = Logic()
            if r0 > 10.:
                bond_logic.long = True
            else:
                base = [self.u.atoms[b].resindex for b in [a1, a2]]
                pair = [self.Fbp_full.get(b, None) for b in base]
                res = base + pair

                is_neighbor = (abs(base[0] - base[1]) == 1)
                bond_logic.crossstack = True if is_neighbor else False

                is_ddDNA = all(bp is not None for bp in pair)
                if is_ddDNA:
                    is_crossstack = any(abs(b - p) == 1
                                        for (b, p) in zip(base, reversed(pair))
                                        )
                    bond_logic.crossstack = True if is_crossstack else False
                    is_Hbond = all(b == p
                                   for (b, p) in zip(base, reversed(pair))
                                   )
                    bond_logic.Hbond = True if is_Hbond else False
                    is_nick = any(resA == self.linker.Fnicks.get(resB, None)
                                  for resA in res
                                  for resB in res
                                  )
                    bond_logic.nick = True if is_nick else False
                    is_co = any(b in self.linker.Fco for b in base)
                    bond_logic.co = True if is_co else False
                else:
                    bond_logic.ssDNA = True
            return bond_logic

        bond_logic = categorize_logic(a1=a1, a2=a2, r0=r0)
        bond_type = set()
        if bond_logic.long:
            bond_type.add("long")
        else:
            bond_type.add("short")
            if bond_logic.strand:
                bond_type.add("strand")
            elif bond_logic.Hbond:
                bond_type.add("Hbond")
            elif bond_logic.crossstack:
                bond_type.add("crossstack")
            elif bond_logic.ssDNA:
                bond_type.add("ssDNA")

            if bond_logic.nick:
                bond_type.add("nick")
            if bond_logic.co:
                bond_type.add("co")

        return bond_type

    def _process_exb(self, exb_file: TextIO) -> set:
        """ create elastic network from file
        -------
         Returns
            -------
            elastic_network
        """
        network = set()
        for line in exb_file:
            split_line = line.split()
            if split_line[0] == "bond":
                a1 = int(split_line[1])
                a2 = int(split_line[2])
                k = float(split_line[3])
                r0 = float(split_line[4])
                bond_type = self._categorize_bond(a1=a1, a2=a2, r0=r0)
                network.add(ENBond(a1=a1, a2=a2, k=k, r0=r0, btype=bond_type))
            elif split_line[0] == "dihedral":
                raise NotImplementedError
            else:
                raise UnexpectedCaseError
        return network

    def _change_modify_logic(self) -> None:
        logic_string = self.linker.project.EN
        self.modify_logic = Logic(long=bool(logic_string[0]),
                                  strand=bool(logic_string[1]),
                                  Hbond=bool(logic_string[2]),
                                  crossstack=bool(logic_string[3]),
                                  nick=bool(logic_string[4]),
                                  co=bool(logic_string[5]),
                                  ssDNA=bool(logic_string[6]),
                                  dihedral=bool(logic_string[7]),
                                  )
        return

    def _modify_en(self) -> set:
        """ create reduced elastic network according to boolean flags
        -------
         Returns
            -------
            reduced_elastic_network
        """
        self._change_modify_logic()
        if self.modify_logic.dihedral:
            _ = self._compute_dihedral()
        logic = attr.asdict(self.modify_logic)
        exclude_type = {name for name, is_active in logic.items() if is_active}

        mod_network = {bond for bond in self.network
                       if any(ex in exclude_type for ex in bond.btype)
                       }
        return mod_network

    def write_en(self) -> None:
        """ write the  modified (by logic) network to file
        """
        mod_network = self._modify_en()
        outfile = self.linker.project.input / self.linker.project.name
        exb_filepath = "{}_{}.exb".format(outfile, self.linker.project.EN)

        with open(exb_filepath, mode="w+") as mod_exb_file:
            for bond in mod_network:
                mod_exb_file.write("{}\n".format(bond))
        return

    def _compute_dihedral(self):
        """ compute restraints corresponding to backbone dihedral
        -------
         Returns
            -------
            dihedral_elastic_network
        """
        return NotImplementedError


@attr.s
class Fit(object):
    project: Project = attr.ib()

    def __attrs_post_init__(self):
        self.infile = self.project.input / self.project.name
        self.u = self._get_universe()
        self.scaffold, self.staples = self._split_strands()

    def _get_universe(self) -> "nd.universe":
        top = self.infile.with_suffix(".psf")
        trj = self.infile.with_suffix(".dcd")
        # TODO: -mid- if pdb, try invoke vmd animate write dcd
        if top.exists() and trj.exists():
            u = mda.Universe(str(top), str(trj))
        else:
            raise FileNotFoundError
        return u

    def _split_strands(self) -> Tuple["nd.segment", List["nd.segment"]]:
        # TODO: -low- multiscaffold
        strands = self.u.segments
        scaffold = max(strands, key=attrgetter("residues.n_residues"))
        staples = [strand for strand in strands if len(
            strand.atoms) != len(scaffold.atoms)]
        return scaffold, staples


@attr.s
class Design(object):
    project: Project = attr.ib()

    def __attrs_post_init__(self):
        self.infile = self.project.input / self.project.name
        self.design = self._get_design()
        self.strands = self.design.strands
        self.scaffold = self._clean_scaffold(self.strands)
        self.excl = self.scaffold[0].strand
        self.staples = self._clean_staple(self.strands)
        self.helixorder = self._create_helix_order()
        self.stapleorder = self._create_staple_order()

    def _is_del(self, base: "nd.base") -> bool:
        return base.num_deletions != 0

    def _clean_scaffold(self, strands: List["nd.base"]) -> List["nd.base"]:
        # TODO: -low- multiscaffold
        # TODO: -low- insertions
        scaffold = [s.tour for s in strands if s.is_scaffold][0]
        scaffold_clean = [b for b in scaffold if not self._is_del(b)]
        return scaffold_clean

    def _clean_staple(self, strands: List["nd.base"]) -> List[List["nd.base"]]:
        # TODO: -low- insertions
        staples = [s.tour for s in strands if not s.is_scaffold]
        staples_clean = [[b for b in s if not self._is_del(b)]
                         for s in staples
                         ]
        return staples_clean

    def _get_design(self) -> Any:
        fil = self.infile.with_suffix(".json")
        seq = self.infile.with_suffix(".seq")
        converter = Converter()
        if fil.exists() and seq.exists():
            converter.read_cadnano_file(file_name=str(fil),
                                        seq_file_name=str(seq),
                                        seq_name=None,
                                        )
        else:
            raise FileNotFoundError
        return converter.dna_structure

    def _create_helix_order(self) -> Dict[int, int]:
        """ helices are not listed in the json in same order as they are listed
            by helix-index. NOTE: enrgMD is helix reorder sensitive.
        -------
            Returns
            -------
            helixorder
        """
        helices_dict = self.design.structure_helices_map
        helixorder = {i: h.load_order for (i, h) in iter(helices_dict.items())}
        return helixorder

    def _create_staple_order(self) -> Dict[int, int]:
        """ enrgMD and nanodesign number staples differently.
            enrgMD: first occurence of staple sorted by h, p
            nanodesign: 5' end of staple sorted by h, p
        -------
            Returns
            -------
            stapleorder
                nanodesign -> enrgMD
        """
        Dhps = [(self.helixorder[s[0].h], s[0].p)
                for s in self.staples
                ]
        Dhps_sorted = sorted(Dhps, key=lambda x: (x[0], x[1]))
        order_ND = [Dhps.index(Dhps_sorted[i])
                    for i, _ in enumerate(Dhps)
                    ]
        stapleorder = {nd: idx for (idx, nd) in enumerate(order_ND)}
        return stapleorder


def get_description() -> str:
    return """links structural information of the cadnano designfile
              [design.json] to fitted atomic model [design.psf, design.dcd].
              stores dictionaries as pickles containg mapping for motifs,
              residue-id, lattice position and base-pairing."""


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
    parser.add_argument("--ENmodify",
                        help="reset nomenclature enrgMD",
                        action="store_true"
                        )
    parser.add_argument("--EN",
                        help="""bool-string:
                                long,strand,Hbond,xstack,nick,co,ss,dhdrl""",
                        type=str,
                        default="11111110",
                        )
    args = parser.parse_args()
    project = Project(input=Path(args.folder),
                      output=Path(args.folder) / "analysis",
                      name=args.name,
                      ENmodify=args.ENmodify,
                      EN=args.EN,
                      )
    with ignored(FileExistsError):
        os.mkdir(project.output)
    return project


def main():
    project = proc_input()

    linker = Linker(project)
    linkage = linker.create_linkage()

    if not project.ENmodify:
        print("linkage output to {}".format(project.output))
        linkage.dump_linkage(project=project)

    else:
        print("modifying extrabonds")
        en = ElaticNetwortModifier(linker)
        en.write_en()


if __name__ == "__main__":
    main()
