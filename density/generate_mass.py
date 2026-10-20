"""test_generate_mass.py"""

import glob
from mpi4py import MPI
import numpy as np
import pmmoto
import os
import time


# This maps from from the lammps input file id and charge
# to a unique id. Atom_map.txt relates ids to atom types
# type:mass
atom_id_mass_map = {
    1: 12.01,
    3: 12.01,
    4: 16,
    5: 1.008,
    7: 14.01,
    8: 1.008,
    12: 16,
    14: 1.008,
}

water_id_mass_map = {15: 16, 16: 1.008}


def initialize_domain():
    """
    Initialize the membrane domain
    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Full domain with reservoirs
    box = [
        [0.0, 176],
        [0.0, 176],
        [-287, 237],  # Ignore water reservoirs
    ]

    sd = pmmoto.initialize(
        voxels=(100, 100, 100),
        box=box,
        rank=rank,
        subdomains=(2, 2, 2),
        boundary_types=((2, 2), (2, 2), (2, 2)),
    )

    return sd


def check_files(files_1, files_2):
    """
    Ensure the membrane and water files match in length and timestamps.
    """
    if len(files_1) != len(files_2):
        raise ValueError("File lengths do not match.")

    files_1.sort()
    files_2.sort()

    for f1, f2 in zip(files_1, files_2):
        value = -1
        if f1.split(".")[-1] == "gz":
            value = -2
        time_1 = int(f1.split(".")[value])
        time_2 = int(f2.split(".")[value])

        if time_1 != time_2:
            raise ValueError(f"These times dont match! {time_1} {time_2}")


def generate_masses(bridges):
    """
    Generates masses.
    """

    sd = initialize_domain()

    membrane_atom_label_file = "rdf/atom_map.txt"
    atom_labels_to_name = pmmoto.io.data_read.read_atom_map(membrane_atom_label_file)

    # Collect unique elements for van der Waals
    elements = list(
        {atom_data["element"] for atom_data in atom_labels_to_name.values()}
    )

    if bridges:
        membrane_files = glob.glob(
            "/ocean/projects/cts200024p/rvickers/RV-P3/64x/03_90/perm_v2/membrane/*"
        )
        water_files = glob.glob(
            "/ocean/projects/cts200024p/rvickers/RV-P3/64x/03_90/perm_v2/pressure/*"
        )

        # Remove specific unwanted files
        unwanted_files = [
            "/ocean/projects/cts200024p/rvickers/RV-P3/64x/03_90/perm_v2/pressure/pressuredata.100000000.gz",
            "/ocean/projects/cts200024p/rvickers/RV-P3/64x/03_90/perm_v2/pressure/pressuredata.110005000.gz",
            "/ocean/projects/cts200024p/rvickers/RV-P3/64x/03_90/perm_v2/pressure/pressuredata.99995000.gz",
        ]
        for file in unwanted_files:
            if file in water_files:
                water_files.remove(file)

    else:
        membrane_files = glob.glob("data/membrane_data/*")
        water_files = glob.glob("data/water_data/*")

    check_files(membrane_files, water_files)
    num_bins = 500
    start = sd.domain.box[2][0]
    end = sd.domain.box[2][1]

    # Initialize collection of bins
    membrane_bins = pmmoto.analysis.bins.Bin(start, end, num_bins, name="membrane")
    water_bins = pmmoto.analysis.bins.Bin(start, end, num_bins, name="water")

    membrane_radii = {}
    for label in atom_id_mass_map.keys():
        membrane_radii[label] = 1e-6

    water_radii = {}
    for label in water_id_mass_map.keys():
        water_radii[label] = 1e-6

    for n_file, (membrane_file, water_file) in enumerate(
        zip(membrane_files, water_files)
    ):

        if sd.rank == 0:
            iter_time = time.time()

        membrane_positions, membrane_atom_type, _, _ = (
            pmmoto.io.data_read.read_lammps_atoms(membrane_file)
        )

        water_positions, water_atom_type, _, _ = pmmoto.io.data_read.read_lammps_atoms(
            water_file
        )

        membrane = pmmoto.particles.initialize_atoms(
            sd,
            membrane_positions,
            membrane_radii,
            membrane_atom_type,
            atom_id_mass_map,
            by_type=False,
            trim_within=True,
        )

        coords = membrane.return_coordinates()
        masses = membrane.return_masses()

        pmmoto.analysis.bins.sum_masses(
            coordinates=coords,
            dimension=2,
            bin=membrane_bins,
            masses=masses,
            subdomain=sd,
        )

        water = pmmoto.particles.initialize_atoms(
            sd,
            water_positions,
            water_radii,
            water_atom_type,
            water_id_mass_map,
            by_type=False,
            trim_within=True,
        )

        coords = water.return_coordinates()
        masses = water.return_masses()

        pmmoto.analysis.bins.sum_masses(
            coordinates=coords,
            dimension=2,
            bin=water_bins,
            masses=masses,
            subdomain=sd,
        )

        if sd.rank == 0:
            print(
                f"Processed file {n_file} in {time.time() - iter_time} seconds. File {membrane_file}",
                flush=True,
            )

            if n_file > 0 and n_file % 250 == 0:
                if sd.rank == 0:
                    print(
                        f"Saving results after {n_file} with filename {membrane_file}",
                        flush=True,
                    )
                membrane_bins.save_bin(sd, "data_out/membrane_bins/")
                water_bins.save_bin(sd, "data_out/water_bins/")

    # Final save
    membrane_bins.save_bin(sd, "data_out/membrane_bins/")
    water_bins.save_bin(sd, "data_out/water_bins/")


if __name__ == "__main__":
    bridges = True
    generate_masses(bridges)
