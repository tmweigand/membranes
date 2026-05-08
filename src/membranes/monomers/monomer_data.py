"""
monomer_data.py
"""

import json
from dataclasses import dataclass, asdict
from typing import List


@dataclass
class Atom:
    id: int
    mol: int
    type: int
    charge: float
    eps: float
    sigma: float
    x: float
    y: float
    z: float
    name: str
    mass: float = 0.0


@dataclass
class Bond:
    id: int
    type: tuple[str, str]
    atom_ids: tuple[int, int]
    atom_names: tuple[str, str]
    k: float
    r_eq: float


@dataclass
class Angle:
    id: int
    type: int
    atom_ids: tuple[int, int, int]
    atom_names: tuple[str, str, str]
    k: float
    theta_eq: float


@dataclass
class Dihedral:
    id: int
    type: int
    atom_ids: tuple[int, int, int, int]
    atom_names: tuple[str, str, str, str]
    k: float
    n: int
    d: int
    weight: float


@dataclass
class Improper:
    id: int
    type: int
    atom_ids: tuple[int, int, int, int]
    atom_names: tuple[str, str, str, str]
    k: float
    n: int
    d: int


@dataclass
class MoleculeFile:
    name: str
    atoms: List[Atom]
    bonds: List[Bond]
    angles: List[Angle]
    dihedrals: List[Dihedral]
    impropers: List[Improper]

    @classmethod
    def from_molecule(cls, mol) -> "MoleculeFile":
        return cls(
            name=mol.name,
            atoms=mol.atoms,
            bonds=mol.bonds,
            angles=mol.angles,
            dihedrals=mol.dihedrals,
            impropers=mol.impropers,
        )

    def save(self, path: str) -> None:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(asdict(self), f, indent=2)
