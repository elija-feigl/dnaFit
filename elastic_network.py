#!/usr/bin/env python
# -*- coding: utf-8 -*-3#
import MDAnalysis as mda  # type:ignore
import attr

from typing import Set, Dict, TextIO

from linker import Linker
from utils import UnexpectedCaseError


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
        self.modify_logic = Logic(long=bool(int(logic_string[0])),
                                  strand=bool(int(logic_string[1])),
                                  Hbond=bool(int(logic_string[2])),
                                  crossstack=bool(int(logic_string[3])),
                                  nick=bool(int(logic_string[4])),
                                  co=bool(int(logic_string[5])),
                                  ssDNA=bool(int(logic_string[6])),
                                  dihedral=bool(int(logic_string[7])),
                                  )

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

    def _compute_dihedral(self):
        """ compute restraints corresponding to backbone dihedral
        -------
         Returns
            -------
            dihedral_elastic_network
        """
        return NotImplementedError
