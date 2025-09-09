import glob
import gzip
import numpy as np
import matplotlib.pyplot as plt
import pmmoto


import rdf_helpers


def py_read_lammps_atoms_velocity(input_file, include_mass=False):
    """
    Read position of atoms from LAMMPS file
    atom_map must sync with LAMMPS ID
    """

    pmmoto.io.io_utils.check_file(input_file)

    if input_file.endswith(".gz"):
        domain_file = gzip.open(input_file, "rt")
    else:
        domain_file = open(input_file, "r", encoding="utf-8")

    charges = {}

    lines = domain_file.readlines()
    domain_data = np.zeros([3, 2], dtype=np.double)
    count_atom = 0
    for n_line, line in enumerate(lines):
        if n_line == 1:
            time_step = float(line)
        elif n_line == 3:
            num_objects = int(line)
            atom_position = np.zeros([num_objects, 3], dtype=np.double)
            atom_velocity = np.zeros([num_objects, 3], dtype=np.double)
            atom_type = np.zeros(num_objects, dtype=int)
            if include_mass:
                masses = np.zeros(num_objects, dtype=float)
        elif 5 <= n_line <= 7:
            domain_data[n_line - 5, 0] = float(line.split(" ")[0])
            domain_data[n_line - 5, 1] = float(line.split(" ")[1])
        elif n_line >= 9:
            split = line.split(" ")

            type = int(split[2])
            atom_type[count_atom] = type
            charge = float(split[4])
            if type in charges:
                if charge not in charges[type]:
                    charges[type].append(charge)
            else:
                charges[type] = [charge]

            if include_mass:
                masses[count_atom] = float(split[3])

            for count, n in enumerate([5, 6, 7]):
                atom_position[count_atom, count] = float(split[n])  # x,y,z,atom_id

            for count, n in enumerate([11, 12, 13]):
                atom_velocity[count_atom, count] = float(split[n])  # x,y,z,atom_id

            count_atom += 1

    domain_file.close()

    if include_mass:
        return atom_position, atom_type, atom_velocity, masses, domain_data
    else:
        return atom_position, atom_type, atom_velocity, domain_data


def water_velocity(water_file):
    """
    Count the water atoms
    """
    water_positions, water_atom_type, atom_velocity, _ = py_read_lammps_atoms_velocity(
        water_file
    )

    print(atom_velocity[0, :])
    vel_magnitudes = np.linalg.norm(atom_velocity, axis=1)
    print(np.max(vel_magnitudes), np.average(vel_magnitudes))


if __name__ == "__main__":

    bridges = False

    # Grab Files
    if bridges:
        voxels_in = (3520, 3520, 4000)
        membrane_files, _ = rdf_helpers.get_bridges_files()
    else:
        voxels_in = (1200, 1200, 1200)
        water_files = glob.glob("data/water_data/*")

    for water_file in water_files:
        water_velocity(water_file)
