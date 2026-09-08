"""
molecule.py
"""

from typing import List, Dict, Optional
from collections import defaultdict

from .monomer_data import Atom, Bond, Angle, Dihedral, Improper


class Molecule:

    def __init__(self, name, atoms, bonds, angles, dihedrals, impropers):

        self.name = name
        self.atoms = atoms
        self.bonds = bonds
        self.angles = angles
        self.dihedrals = dihedrals
        self.impropers = impropers

        self.n_atoms = len(atoms)
        self.n_bonds = len(bonds)
        self.n_angles = len(angles)
        self.n_dihedrals = len(dihedrals)
        self.n_impropers = len(impropers)

        self.n_atom_types = len({atom.type for atom in atoms})
        self.n_bond_types = len({bond.type for bond in bonds})
        self.n_angle_types = len({angle.type for angle in angles})
        self.n_dihedral_types = len({dihedral.type for dihedral in dihedrals})
        self.n_improper_types = len({improper.type for improper in impropers})

        self.masses = None
        self.pair_coeffs = None
        self.bond_coeffs = None
        self.angle_coeffs = None
        self.dihedral_coeffs = None
        self.improper_coeffs = None

        # atom lookup tables
        self.atom_id_to_atom = {atom.id: atom.name for atom in atoms}
        self.atom_name_to_ids = defaultdict(list)
        for atom in atoms:
            self.atom_name_to_ids[atom.name].append(atom.id)

    def get_atom_ids(self, atom_name: str) -> list[int]:
        """Return all atom IDs that match a given atom name."""
        return list(self.atom_name_to_ids.get(atom_name, []))

    def get_atom(self, atom_id: int) -> Atom:
        """Return Atom by ID."""
        return self.atom_id_to_atom[atom_id]

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def write_mol(self, path: str) -> None:
        """Write a LAMMPS molecule file with type labels."""
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(self._mol_lines()))
            fh.write("\n")

    def _mol_lines(self) -> list[str]:
        def _lbl(t: tuple) -> str:
            return "-".join(t)

        lines: list[str] = [
            f"LAMMPS molecule file: {self.name}",
            "",
            f"{self.n_atoms} atoms",
            f"{self.n_bonds} bonds",
            f"{self.n_angles} angles",
            f"{self.n_dihedrals} dihedrals",
            f"{self.n_impropers} impropers",
            "",
        ]

        # Types
        lines += ["Types", ""]
        for atom in self.atoms:
            lines.append(f"  {atom.id:>6}  {atom.type}")
        lines.append("")

        # Charges
        lines += ["Charges", ""]
        for atom in self.atoms:
            lines.append(f"  {atom.id:>6}  {atom.charge:>14.7f}")
        lines.append("")

        # Coords
        lines += ["Coords", ""]
        for atom in self.atoms:
            lines.append(
                f"  {atom.id:>6}" f"  {atom.x:>13.7f}  {atom.y:>13.7f}  {atom.z:>13.7f}"
            )
        lines.append("")

        # Bonds
        if self.n_bonds > 0:
            lines += ["Bonds", ""]
            for bond in self.bonds:
                lines.append(
                    f"  {bond.id:>6}  {_lbl(bond.type)}"
                    f"  {bond.atom_ids[0]:>6}  {bond.atom_ids[1]:>6}"
                )
            lines.append("")

        # Angles
        if self.n_angles > 0:
            lines += ["Angles", ""]
            for angle in self.angles:
                lines.append(
                    f"  {angle.id:>6}  {_lbl(angle.type)}"
                    f"  {angle.atom_ids[0]:>6}  {angle.atom_ids[1]:>6}  {angle.atom_ids[2]:>6}"
                )
            lines.append("")

        # Dihedrals
        if self.n_dihedrals > 0:
            lines += ["Dihedrals", ""]
            for dihedral in self.dihedrals:
                lines.append(
                    f"  {dihedral.id:>6}  {_lbl(dihedral.type)}"
                    f"  {dihedral.atom_ids[0]:>6}  {dihedral.atom_ids[1]:>6}"
                    f"  {dihedral.atom_ids[2]:>6}  {dihedral.atom_ids[3]:>6}"
                )
            lines.append("")

        # Impropers
        if self.n_impropers > 0:
            lines += ["Impropers", ""]
            for improper in self.impropers:
                lines.append(
                    f"  {improper.id:>6}  {_lbl(improper.type)}"
                    f"  {improper.atom_ids[0]:>6}  {improper.atom_ids[1]:>6}"
                    f"  {improper.atom_ids[2]:>6}  {improper.atom_ids[3]:>6}"
                )
            lines.append("")

        return lines


def collect_atom_data(
    pmd_structure,
    name,
    mol_id: int = 1,
    custom_charges: Optional[Dict[int, float]] = None,
) -> list:
    """
    From self.pmd_structure collect the Atom information.

    Assume LJ pair style

    """
    cc = custom_charges or {}
    atoms = []
    for a in pmd_structure.atoms:
        _2_1_6 = 2.0 ** (1.0 / 6.0)
        atom = Atom(
            id=a.idx + 1,
            mol=mol_id,
            type=a.atom_type.name,
            charge=cc.get(a.idx + 1, a.charge),
            eps=getattr(a.atom_type, "epsilon"),
            sigma=2.0 * getattr(a.atom_type, "rmin") / _2_1_6,
            x=a.xx,
            y=a.xy,
            z=a.xz,
            name=name + "-" + a.atom_type.name,
            mass=a.mass,
        )
        atoms.append(atom)

    return atoms


def collect_bond_data(pmd_structure) -> list:
    """
    From self.pmd_structure collect the Bond information

    Assume harmonic style given as:

        E = k*(r - r_eq)^2

    """
    bonds = []
    for i, b in enumerate(pmd_structure.bonds, 1):
        key = tuple(sorted((b.atom1.atom_type.name, b.atom2.atom_type.name)))
        if b.type:
            bond = Bond(
                id=i,
                type=key,
                atom_ids=(b.atom1.idx + 1, b.atom2.idx + 1),
                atom_names=(b.atom1.name, b.atom2.name),
                k=getattr(b.type, "k"),
                r_eq=getattr(b.type, "req"),
            )
            bonds.append(bond)

    return bonds


def collect_angle_data(pmd_structure):
    """From self.pmd_structure collect the Angle information.

    Order matters! - atom2 is the center atom

    Assumes harmonic angle style given as::

        E = k*(theta - theta_eq)^2

    """
    angles = []
    for i, a in enumerate(pmd_structure.angles, 1):
        if not a.type:
            continue
        atype = a.type
        angle_type_key = (
            a.atom1.atom_type.name,
            a.atom2.atom_type.name,
            a.atom3.atom_type.name,
        )
        angle = Angle(
            id=i,
            type=angle_type_key,
            atom_ids=(a.atom1.idx + 1, a.atom2.idx + 1, a.atom3.idx + 1),
            atom_names=(a.atom1.name, a.atom2.name, a.atom3.name),
            k=getattr(a.type, "k"),
            theta_eq=getattr(a.type, "theteq"),
        )
        angles.append(angle)
    return angles


def collect_dihedral_data(pmd_structure):
    """From self.pmd_structure collect the Dihedral information

    Assumes charmm style given as:

        E = K*(1 + cos(n*theta - d))

    """
    dihedrals = []
    for i, d in enumerate(pmd_structure.dihedrals, 1):
        dihedral_type_key = (
            d.atom1.atom_type.name,
            d.atom2.atom_type.name,
            d.atom3.atom_type.name,
            d.atom4.atom_type.name,
        )
        dihedral = Dihedral(
            id=i,
            type=dihedral_type_key,
            atom_ids=(
                d.atom1.idx + 1,
                d.atom2.idx + 1,
                d.atom3.idx + 1,
                d.atom4.idx + 1,
            ),
            atom_names=(
                d.atom1.name,
                d.atom2.name,
                d.atom3.name,
                d.atom4.name,
            ),
            k=getattr(d.type, "phi_k"),
            n=getattr(d.type, "per"),
            d=int(getattr(d.type, "phase")),
            weight=getattr(d.type, "scee"),
        )
        dihedrals.append(dihedral)
    return dihedrals


def collect_improper_data(pmd_structure):
    """From self.pmd_structure collect the Improper information

    GAFF impropers use a periodic (cvff) style:

        E = K*(1 + cos(n*chi - d))
    """
    impropers = []
    for i, im in enumerate(pmd_structure.impropers, 1):
        improper_type_key = (
            im.atom1.atom_type.name,
            im.atom2.atom_type.name,
            im.atom3.atom_type.name,
            im.atom4.atom_type.name,
        )
        improper = Improper(
            id=i,
            type=improper_type_key,
            atom_ids=(
                im.atom1.idx + 1,
                im.atom2.idx + 1,
                im.atom3.idx + 1,
                im.atom4.idx + 1,
            ),
            atom_names=(
                im.atom1.name,
                im.atom2.name,
                im.atom3.name,
                im.atom4.name,
            ),
            k=getattr(im.type, "phi_k"),
            n=getattr(im.type, "per"),
            d=int(getattr(im.type, "phase")),
        )
        impropers.append(improper)
    return impropers


def convert_to_molecule(
    name,
    pmd_structure,
    mol_id: int = 1,
    custom_charges: Optional[Dict[int, float]] = None,
) -> Molecule:

    cc = custom_charges or {}
    atoms = collect_atom_data(pmd_structure, name, mol_id, cc)
    bonds = collect_bond_data(pmd_structure)
    angles = collect_angle_data(pmd_structure)
    dihedrals = collect_dihedral_data(pmd_structure)
    impropers = collect_improper_data(pmd_structure)
    return Molecule(name, atoms, bonds, angles, dihedrals, impropers)
