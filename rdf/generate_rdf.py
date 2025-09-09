"""test_generate_rdf.py"""

import glob
from mpi4py import MPI
import numpy as np
import pmmoto
import time
import rdf_helpers

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
        [0.0, 176],
        [0.0, 176],
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


def generate_rdf(bridges, extended):
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
        membrane_files, water_files = rdf_helpers.get_bridges_files()
    else:
        membrane_files = glob.glob("data/membrane_data/*")
        water_files = glob.glob("data/water_data/*")

    rdf_helpers.check_files(membrane_files, water_files)

    num_labels = len(atom_labels_to_name)
    num_bins = [280] * num_labels

    water_radius = 1.4
    if extended:
        nnn = 5
        water_radius = water_radius * nnn
        num_bins = [280 * nnn] * num_labels

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

        membrane_atom_counts = membrane.get_own_count(sd, True)

        _sum = 0
        for counts in membrane_atom_counts.values():
            _sum += counts

        water = pmmoto.particles.initialize_atoms(
            sd,
            water_positions,
            water_radii,
            water_atom_type,
            by_type=True,
            trim_within=True,
        )

        # Remove the hydrogens. The oxygens are type 15
        # labeling is so confusing here
        water = water.return_list(15)

        water_oxygen_counts = water.get_own_count(sd, True)

        volume = sd.domain.get_volume()

        print(membrane_atom_counts, water_oxygen_counts, volume)  # / volume)

        # density = {}
        # for atom, count in membrane_atom_counts.items():
        #     density[atom] = count / volume

        pmmoto.domain_generation.rdf.bin_distances(
            subdomain=sd, probe_atom_list=water, atoms=membrane, bins=bins
        )

        if sd.rank == 0:
            _time = time.time() - iter_time
            print(
                f"Processed file {n_file} in {_time} seconds. File {membrane_file}",
                flush=True,
            )

        if n_file > 0 and n_file % 250 == 0:
            if sd.rank == 0:
                print(
                    f"Saving results after {n_file} with filename {membrane_file}",
                    flush=True,
                )
            # if extended:
            #     bins.save_bins(sd, "data_out/bins_extended/")
            # else:
            #     bins.save_bins(sd, "data_out/bins/")

    # Final save
    # if extended:
    #     bins.save_bins(sd, "data_out/bins_extended/")
    # else:
    #     bins.save_bins(sd, "data_out/bins/")


if __name__ == "__main__":
    bridges = False
    extended = True
    generate_rdf(bridges, extended)
