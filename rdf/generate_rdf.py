"""test_generate_rdf.py"""

import glob
from mpi4py import MPI
import numpy as np
import pmmoto
import os
import time


# This maps from from the lammps input file id and charge
# to a unique id. Atom_map.txt relates ids to atom types
atom_id_charge_map = {
    (1, 0.6797): 1,
    (1, 0.743425): 2,
    (3, -0.23): 3,
    (3, -0.1956): 4,
    (3, -0.1565): 5,
    (3, 0.014): 6,
    (3, 0.1716): 7,
    (4, -0.587509): 8,
    (5, 0.10745): 9,
    (5, 0.131): 10,
    (5, 0.1816): 11,
    (7, -0.4621): 12,
    (7, -0.398375): 13,
    (8, 0.23105): 14,
    (12, -0.5351): 15,
    (14, 0.4315): 16,
}


def initialize_domain():
    """
    Initialize the membrane domain
    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Full domain with reservoirs
    box = [
        [0.0, 175.768],
        [0.0, 175.768],
        [-150, 175],  # Ignore water reservoirs
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


def save_rdf(sd, bins, rdf_bins, atom_labels_to_name, time):
    """
    Save the RDF data to files.

    Parameters
    ----------
    sd : SubDomain
        The subdomain object containing rank information
    bins : Bins
        The bins object containing bin centers
    rdf_bins : dict
        Dictionary containing RDF values for each atom type
    atom_labels_to_name : dict
        Dictionary mapping atom labels to atom names
    """
    if sd.rank == 0:
        # Create output directory if it doesn't exist
        os.makedirs("data_out/rdf_data", exist_ok=True)

        # Save data for each atom type
        for label, rdf in rdf_bins.rdf_bins.items():
            atom_name = atom_labels_to_name[label]["label"]
            centers = bins.bin_centers[label]

            # Stack the data into columns
            data = np.column_stack((centers, rdf))

            # Save with header
            np.savetxt(
                f"data_out/rdf_data/rdf_{atom_name}.txt",
                data,
                header=f"Distance(Å) RDF at {time}",
                delimiter="\t",
                comments="#",
            )


def generate_rdf(bridges):
    """
    Test for generating a radial distribution function from LAMMPS data
    """

    sd = initialize_domain()

    membrane_atom_label_file = "rdf/atom_map.txt"
    atom_labels_to_name = pmmoto.io.data_read.read_atom_map(membrane_atom_label_file)

    # Collect unique elements for van der Waals
    elements = list(
        {atom_data["element"] for atom_data in atom_labels_to_name.values()}
    )

    uff_radii = pmmoto.particles.uff_radius(atom_names=elements)

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

    num_labels = len(atom_labels_to_name)
    num_bins = [280] * num_labels
    water_radius = 1.4

    membrane_radii = np.zeros(num_labels)
    m_radii = {}
    labels = np.zeros(num_labels, dtype=int)
    names = []
    for n, atom_data in enumerate(atom_labels_to_name.values()):
        element = atom_data["element"]
        element_number = pmmoto.particles.convert_atoms_elements_to_ids([element])
        membrane_radii[n] = uff_radii[element_number[0]] + water_radius
        m_radii[n + 1] = uff_radii[element_number[0]] + water_radius
        labels[n] = n + 1
        names.append(atom_data["label"])

    bins = pmmoto.analysis.bins.Bins(
        starts=np.zeros(num_labels),
        ends=membrane_radii,
        num_bins=num_bins,
        labels=labels,
        names=names,
    )

    for n_file, (membrane_file, water_file) in enumerate(
        zip(membrane_files, water_files)
    ):

        if sd.rank == 0:
            iter_time = time.time()

        membrane_positions, membrane_atom_type, _, _ = (
            pmmoto.io.data_read.read_lammps_atoms(membrane_file, atom_id_charge_map)
        )

        water_positions, water_atom_type, _, _ = pmmoto.io.data_read.read_lammps_atoms(
            water_file
        )

        water_radii = {15: water_radius, 16: water_radius}

        membrane = pmmoto.particles.initialize_atoms(
            sd,
            membrane_positions,
            m_radii,
            membrane_atom_type,
            by_type=True,
            add_periodic=True,
            set_own=True,
        )

        water = pmmoto.particles.initialize_atoms(
            sd,
            water_positions,
            water_radii,
            water_atom_type,
            by_type=True,
            trim_within=True,
        )

        water = water.return_list(15)

        pmmoto.domain_generation.rdf.bin_distances(
            subdomain=sd, probe_atom_list=water, atoms=membrane, bins=bins
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
            bins.save_bins(sd, "data_out/bins/")

    # Final save
    bins.save_bins(sd, "data_out/bins/")


if __name__ == "__main__":
    bridges = True
    generate_rdf(bridges)
