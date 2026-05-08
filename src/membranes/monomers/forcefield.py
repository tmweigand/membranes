"""
forcefield.py
"""

from __future__ import annotations


def _type_label(t: tuple) -> str:
    """Join a type-key tuple into a human-readable string, e.g. ('ca','c3') -> 'ca-c3'."""
    return "-".join(t)


class ForceField:
    """
    Merges force-field parameters from one or more Molecule objects and
    writes a LAMMPS data file with type label sections, coefficient sections,
    and full topology (Atoms, Bonds, Angles, Dihedrals, Impropers).

    Parameters
    ----------
    molecules : Molecule | list[Molecule]
        One or more fully-populated Molecule objects (atoms must have mass).
    """

    def __init__(self, molecules):
        from .molecule import Molecule  # local import avoids circular dependency

        if isinstance(molecules, Molecule):
            molecules = [molecules]
        self._build(list(molecules))

    # ------------------------------------------------------------------
    # Build
    # ------------------------------------------------------------------

    def _build(self, molecules) -> None:
        self.molecules = molecules
        self.source_names: list[str] = [mol.name for mol in molecules]

        # ordered dicts: label -> params  (first occurrence wins)
        atom_params: dict[str, tuple] = {}  # label -> (mass, eps, sigma)
        bond_params: dict[str, tuple] = {}  # label -> (k, r_eq)
        angle_params: dict[str, tuple] = {}  # label -> (k, theta_eq)
        dihedral_params: dict[str, tuple] = {}  # label -> (k, n, d, weight)
        improper_params: dict[str, tuple] = {}  # label -> (k, n, d)

        for mol in molecules:
            for atom in mol.atoms:
                if atom.type not in atom_params:
                    atom_params[atom.type] = (atom.mass, atom.eps, atom.sigma)

            for bond in mol.bonds:
                label = _type_label(bond.type)
                if label not in bond_params:
                    bond_params[label] = (bond.k, bond.r_eq)

            for angle in mol.angles:
                label = _type_label(angle.type)
                if label not in angle_params:
                    angle_params[label] = (angle.k, angle.theta_eq)

            for dihedral in mol.dihedrals:
                label = _type_label(dihedral.type)
                if label not in dihedral_params:
                    dihedral_params[label] = (
                        dihedral.k,
                        dihedral.n,
                        dihedral.d,
                        dihedral.weight,
                    )

            for improper in mol.impropers:
                label = _type_label(improper.type)
                if label not in improper_params:
                    improper_params[label] = (improper.k, improper.n, improper.d)

        # 1-based integer ID -> label
        self.atom_type_labels: dict[int, str] = {
            i + 1: k for i, k in enumerate(atom_params)
        }
        self.bond_type_labels: dict[int, str] = {
            i + 1: k for i, k in enumerate(bond_params)
        }
        self.angle_type_labels: dict[int, str] = {
            i + 1: k for i, k in enumerate(angle_params)
        }
        self.dihedral_type_labels: dict[int, str] = {
            i + 1: k for i, k in enumerate(dihedral_params)
        }
        self.improper_type_labels: dict[int, str] = {
            i + 1: k for i, k in enumerate(improper_params)
        }

        # label -> 1-based integer ID  (reverse lookups)
        self._atom_label_to_id: dict[str, int] = {
            lbl: tid for tid, lbl in self.atom_type_labels.items()
        }
        self._bond_label_to_id: dict[str, int] = {
            lbl: tid for tid, lbl in self.bond_type_labels.items()
        }
        self._angle_label_to_id: dict[str, int] = {
            lbl: tid for tid, lbl in self.angle_type_labels.items()
        }
        self._dihedral_label_to_id: dict[str, int] = {
            lbl: tid for tid, lbl in self.dihedral_type_labels.items()
        }
        self._improper_label_to_id: dict[str, int] = {
            lbl: tid for tid, lbl in self.improper_type_labels.items()
        }

        # 1-based integer ID -> params
        self.atom_params: dict[int, tuple] = {
            i + 1: v for i, v in enumerate(atom_params.values())
        }
        self.bond_params: dict[int, tuple] = {
            i + 1: v for i, v in enumerate(bond_params.values())
        }
        self.angle_params: dict[int, tuple] = {
            i + 1: v for i, v in enumerate(angle_params.values())
        }
        self.dihedral_params: dict[int, tuple] = {
            i + 1: v for i, v in enumerate(dihedral_params.values())
        }
        self.improper_params: dict[int, tuple] = {
            i + 1: v for i, v in enumerate(improper_params.values())
        }

    # ------------------------------------------------------------------
    # Public lookups
    # ------------------------------------------------------------------

    def atom_type_id(self, label: str) -> int:
        return self._atom_label_to_id[label]

    def bond_type_id(self, label: str) -> int:
        return self._bond_label_to_id[label]

    def angle_type_id(self, label: str) -> int:
        return self._angle_label_to_id[label]

    def dihedral_type_id(self, label: str) -> int:
        return self._dihedral_label_to_id[label]

    def improper_type_id(self, label: str) -> int:
        return self._improper_label_to_id[label]

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def write(self, path: str) -> None:
        """Write a LAMMPS data file with type labels, coefficients, and topology."""
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(self._build_lines()))
            fh.write("\n")

    def _build_lines(self) -> list[str]:
        total_atoms = sum(mol.n_atoms for mol in self.molecules)
        total_bonds = sum(mol.n_bonds for mol in self.molecules)
        total_angles = sum(mol.n_angles for mol in self.molecules)
        total_dihedrals = sum(mol.n_dihedrals for mol in self.molecules)
        total_impropers = sum(mol.n_impropers for mol in self.molecules)

        mols = ", ".join(self.source_names)
        lines: list[str] = [
            "# LAMMPS master force field",
            f"# Source molecules: {mols}",
            "",
            f"{total_atoms} atoms",
            f"{total_bonds} bonds",
            f"{total_angles} angles",
            f"{total_dihedrals} dihedrals",
            f"{total_impropers} impropers",
            "",
            f"{len(self.atom_type_labels)} atom types",
            f"{len(self.bond_type_labels)} bond types",
            f"{len(self.angle_type_labels)} angle types",
            f"{len(self.dihedral_type_labels)} dihedral types",
            f"{len(self.improper_type_labels)} improper types",
            "",
        ]

        def _label_section(title: str, labels: dict[int, str]) -> None:
            if not labels:
                return
            lines.append(title)
            lines.append("")
            for tid, lbl in labels.items():
                lines.append(f"  {tid}  {lbl}")
            lines.append("")

        _label_section("Atom Type Labels", self.atom_type_labels)
        _label_section("Bond Type Labels", self.bond_type_labels)
        _label_section("Angle Type Labels", self.angle_type_labels)
        _label_section("Dihedral Type Labels", self.dihedral_type_labels)
        _label_section("Improper Type Labels", self.improper_type_labels)

        # Masses  -- label as identifier
        if self.atom_params:
            lines += ["Masses", ""]
            for tid, (mass, eps, sigma) in self.atom_params.items():
                lbl = self.atom_type_labels[tid]
                lines.append(f"  {lbl}  {mass:.4f}")
            lines.append("")

        # Pair Coeffs  lj/cut:  epsilon  sigma  -- label as identifier
        if self.atom_params:
            lines += ["Pair Coeffs  # lj/cut", ""]
            for tid, (mass, eps, sigma) in self.atom_params.items():
                lbl = self.atom_type_labels[tid]
                lines.append(f"  {lbl}  {eps:.6f}  {sigma:.6f}")
            lines.append("")

        # Bond Coeffs  harmonic:  k  r0  -- label as identifier
        if self.bond_params:
            lines += ["Bond Coeffs  # harmonic", ""]
            for tid, (k, r_eq) in self.bond_params.items():
                lbl = self.bond_type_labels[tid]
                lines.append(f"  {lbl}  {k:.4f}  {r_eq:.4f}")
            lines.append("")

        # Angle Coeffs  harmonic:  k  theta0  -- label as identifier
        if self.angle_params:
            lines += ["Angle Coeffs  # harmonic", ""]
            for tid, (k, theta_eq) in self.angle_params.items():
                lbl = self.angle_type_labels[tid]
                lines.append(f"  {lbl}  {k:.4f}  {theta_eq:.4f}")
            lines.append("")

        # Dihedral Coeffs  charmm:  k  n  d  weight  -- label as identifier
        if self.dihedral_params:
            lines += ["Dihedral Coeffs  # charmm", ""]
            for tid, (k, n, d, weight) in self.dihedral_params.items():
                lbl = self.dihedral_type_labels[tid]
                lines.append(f"  {lbl}  {k:.4f}  {int(n)}  {int(d)}  {weight:.4f}")
            lines.append("")

        # Improper Coeffs  cvff:  k  d  n  -- label as identifier
        # AMBER phase 0° -> d = +1,  180° -> d = -1
        if self.improper_params:
            lines += ["Improper Coeffs  # cvff", ""]
            for tid, (k, n, phase) in self.improper_params.items():
                lbl = self.improper_type_labels[tid]
                d = -1 if round(phase) == 180 else 1
                lines.append(f"  {lbl}  {k:.4f}  {d}  {int(n)}")
            lines.append("")

        # ------------------------------------------------------------------
        # Topology sections
        # ------------------------------------------------------------------

        # Atoms  format: id  mol  type_label  charge  x  y  z
        if total_atoms > 0:
            lines += ["Atoms  # full", ""]
            atom_id = 1
            for mol_idx, mol in enumerate(self.molecules, 1):
                for atom in mol.atoms:
                    lines.append(
                        f"  {atom_id:>6}  {mol_idx:>6}  {atom.type}"
                        f"  {atom.charge:>12.7f}"
                        f"  {atom.x:>13.7f}  {atom.y:>13.7f}  {atom.z:>13.7f}"
                    )
                    atom_id += 1
            lines.append("")

        # Bonds  format: id  type_label  atom1  atom2
        if total_bonds > 0:
            lines += ["Bonds", ""]
            bond_id = 1
            atom_offset = 0
            for mol in self.molecules:
                for bond in mol.bonds:
                    lbl = _type_label(bond.type)
                    a1 = atom_offset + bond.atom_ids[0]
                    a2 = atom_offset + bond.atom_ids[1]
                    lines.append(f"  {bond_id:>6}  {lbl}  {a1:>6}  {a2:>6}")
                    bond_id += 1
                atom_offset += mol.n_atoms
            lines.append("")

        # Angles  format: id  type_label  atom1  atom2  atom3
        if total_angles > 0:
            lines += ["Angles", ""]
            angle_id = 1
            atom_offset = 0
            for mol in self.molecules:
                for angle in mol.angles:
                    lbl = _type_label(angle.type)
                    a1 = atom_offset + angle.atom_ids[0]
                    a2 = atom_offset + angle.atom_ids[1]
                    a3 = atom_offset + angle.atom_ids[2]
                    lines.append(f"  {angle_id:>6}  {lbl}  {a1:>6}  {a2:>6}  {a3:>6}")
                    angle_id += 1
                atom_offset += mol.n_atoms
            lines.append("")

        # Dihedrals  format: id  type_label  a1  a2  a3  a4
        if total_dihedrals > 0:
            lines += ["Dihedrals", ""]
            dihedral_id = 1
            atom_offset = 0
            for mol in self.molecules:
                for dihedral in mol.dihedrals:
                    lbl = _type_label(dihedral.type)
                    a1 = atom_offset + dihedral.atom_ids[0]
                    a2 = atom_offset + dihedral.atom_ids[1]
                    a3 = atom_offset + dihedral.atom_ids[2]
                    a4 = atom_offset + dihedral.atom_ids[3]
                    lines.append(
                        f"  {dihedral_id:>6}  {lbl}"
                        f"  {a1:>6}  {a2:>6}  {a3:>6}  {a4:>6}"
                    )
                    dihedral_id += 1
                atom_offset += mol.n_atoms
            lines.append("")

        # Impropers  format: id  type_label  a1  a2  a3  a4
        if total_impropers > 0:
            lines += ["Impropers", ""]
            improper_id = 1
            atom_offset = 0
            for mol in self.molecules:
                for improper in mol.impropers:
                    lbl = _type_label(improper.type)
                    a1 = atom_offset + improper.atom_ids[0]
                    a2 = atom_offset + improper.atom_ids[1]
                    a3 = atom_offset + improper.atom_ids[2]
                    a4 = atom_offset + improper.atom_ids[3]
                    lines.append(
                        f"  {improper_id:>6}  {lbl}"
                        f"  {a1:>6}  {a2:>6}  {a3:>6}  {a4:>6}"
                    )
                    improper_id += 1
                atom_offset += mol.n_atoms
            lines.append("")

        return lines
