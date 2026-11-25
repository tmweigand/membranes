#!/usr/bin/env python3

import argparse
import os


def get_forcefield_parameters(forcefield):
    if forcefield == "riley":
        return {
            "charges": (-0.84844, 0.42422, 0.42422),
            "masses": (16.000000, 1.008000),
        }

    if forcefield == "tip3p_orig":
        return {
            "charges": (-0.834, 0.417, 0.417),
            "masses": (15.9994, 1.00800),
        }

    raise ValueError("Unsupported force field: {}".format(forcefield))


def write_water_data(output_file, forcefield):
    n_side = 10  ### n_atoms = n_side^3 = 1000, can be modified to have a bigger system but box_length should also be changed.
    box_length = (
        31.1  ### Box length changes should be based on the density of the water.
    )
    spacing = box_length / n_side

    oxygen = (0.0, 0.0, 0.0)
    hydrogen_1 = (1.100, 0.023, -0.026)
    hydrogen_2 = (-0.399, -0.686, 0.762)

    params = get_forcefield_parameters(forcefield)
    q_o, q_h1, q_h2 = params["charges"]

    atoms = []
    bonds = []
    angles = []

    atom_id = 1
    bond_id = 1
    angle_id = 1
    molecule_id = 1

    for ix in range(n_side):
        for iy in range(n_side):
            for iz in range(n_side):
                origin = (
                    (ix + 0.5) * spacing,
                    (iy + 0.5) * spacing,
                    (iz + 0.5) * spacing,
                )

                coordinates = []
                for offset in (oxygen, hydrogen_1, hydrogen_2):
                    coordinates.append(tuple(origin[d] + offset[d] for d in range(3)))

                atoms.append((atom_id, molecule_id, 1, q_o, *coordinates[0]))
                atoms.append((atom_id + 1, molecule_id, 2, q_h1, *coordinates[1]))
                atoms.append((atom_id + 2, molecule_id, 2, q_h2, *coordinates[2]))

                bonds.append((bond_id, 1, atom_id, atom_id + 1))
                bonds.append((bond_id + 1, 1, atom_id, atom_id + 2))

                angles.append((angle_id, 1, atom_id + 1, atom_id, atom_id + 2))

                atom_id += 3
                bond_id += 2
                angle_id += 1
                molecule_id += 1

    output_file = os.path.abspath(output_file)
    output_dir = os.path.dirname(output_file)
    os.makedirs(output_dir, exist_ok=True)

    with open(output_file, "w") as f:
        f.write("Water system generated for {}\n\n".format(forcefield))
        f.write("{} atoms\n".format(len(atoms)))
        f.write("{} bonds\n".format(len(bonds)))
        f.write("{} angles\n\n".format(len(angles)))

        f.write("2 atom types\n")
        f.write("1 bond types\n")
        f.write("1 angle types\n\n")

        f.write("0.0 {:.6f} xlo xhi\n".format(box_length))
        f.write("0.0 {:.6f} ylo yhi\n".format(box_length))
        f.write("0.0 {:.6f} zlo zhi\n\n".format(box_length))

        f.write("Atoms # full\n\n")
        for atom in atoms:
            f.write(
                "{} {} {} {:.6f} {:.6f} {:.6f} {:.6f}\n".format(
                    atom[0],
                    atom[1],
                    atom[2],
                    atom[3],
                    atom[4],
                    atom[5],
                    atom[6],
                )
            )

        f.write("\nBonds\n\n")
        for bond in bonds:
            f.write("{} {} {} {}\n".format(*bond))

        f.write("\nAngles\n\n")
        for angle in angles:
            f.write("{} {} {} {} {}\n".format(*angle))


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--forcefield",
        required=True,
        choices=("riley", "tip3p_orig"),
    )
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    write_water_data(args.output, args.forcefield)

    print("Wrote {} data file: {}".format(args.forcefield, args.output))


if __name__ == "__main__":
    main()
