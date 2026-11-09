import math

# System settings

n_water = 1000
n_side = 10
box_length = 31.1  # Angstrom

water_coords = [
    (0.000,  0.000,  0.000),  # O
    (1.100,  0.023, -0.026),  # H
    (-0.399, -0.686, 0.762),  # H
]

charges = [-0.834, 0.417, 0.417]

if n_side ** 3 != n_water:
    raise ValueError("n_water must equal n_side**3 for this grid.")

spacing = box_length / n_side

atoms = []
bonds = []
angles = []

atom_id = 0
bond_id = 0
angle_id = 0
mol_id = 0

# Place water molecules on a cubic lattice.

for ix in range(n_side):
    for iy in range(n_side):
        for iz in range(n_side):

            mol_id += 1

            ox = (ix + 0.5) * spacing
            oy = (iy + 0.5) * spacing
            oz = (iz + 0.5) * spacing

            molecule_atom_ids = []

            for atom_type, coord, charge in zip(
                [1, 2, 2], water_coords, charges
            ):
                atom_id += 1
                molecule_atom_ids.append(atom_id)

                x = ox + coord[0]
                y = oy + coord[1]
                z = oz + coord[2]

                atoms.append(
                    (atom_id, mol_id, atom_type, charge, x, y, z)
                )

            oxygen_id, h1_id, h2_id = molecule_atom_ids

            bond_id += 1
            bonds.append((bond_id, 1, oxygen_id, h1_id))

            bond_id += 1
            bonds.append((bond_id, 1, oxygen_id, h2_id))

            angle_id += 1
            angles.append(
                (angle_id, 1, h1_id, oxygen_id, h2_id)
            )

# Write LAMMPS data file

with open("water.data", "w") as f:
    f.write("Water-only system generated from tip3p.mol(lammps_parameters) \n\n")

    f.write(f"{len(atoms)} atoms\n")
    f.write(f"{len(bonds)} bonds\n")
    f.write(f"{len(angles)} angles\n")
    f.write("0 dihedrals\n")
    f.write("0 impropers\n\n")

    f.write("2 atom types\n")
    f.write("1 bond types\n")
    f.write("1 angle types\n")
    f.write("0 dihedral types\n")
    f.write("0 improper types\n\n")

    f.write(f"0.0 {box_length:.6f} xlo xhi\n")
    f.write(f"0.0 {box_length:.6f} ylo yhi\n")
    f.write(f"0.0 {box_length:.6f} zlo zhi\n\n")

    f.write("Masses\n\n")
    f.write("1 15.9994  # O\n")
    f.write("2 1.00800   # H\n\n")

    f.write("Atoms # full\n\n")
    for row in atoms:
        aid, mid, atype, q, x, y, z = row
        f.write(
            f"{aid} {mid} {atype} {q:.8f} "
            f"{x:.6f} {y:.6f} {z:.6f}\n"
        )

    f.write("\nBonds\n\n")
    for row in bonds:
        f.write("{} {} {} {}\n".format(*row))

    f.write("\nAngles\n\n")
    for row in angles:
        f.write("{} {} {} {} {}\n".format(*row))

print("Created water.data")
print("Water molecules:", n_water)
print("Atoms:", len(atoms))
print("Bonds:", len(bonds))
print("Angles:", len(angles))
print("Box length:", box_length, "Angstrom")
