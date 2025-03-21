"""test_generate_rdf.py"""

import glob
from mpi4py import MPI
import numpy as np
import pmmoto
import matplotlib.pyplot as plt
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


def gen_radii(atom_ids, value):
    radii = {}
    for _id in atom_ids:
        radii[_id] = value
    return radii


def generate_binned_distance_plots(sd, bins, binned_distances, atom_labels_to_name):

    if sd.rank == 0:
        pmmoto.io.io_utils.check_file_path("data_out/binned_distance/")
        # Set style for publication-quality figures
        # plt.style.use("seaborn")

        # Create figure with specific size
        plt.figure(figsize=(10, 6))

        # Plot RDFs for each atom type
        for label, rdf in binned_distances.rdf_bins.items():
            plt.plot(
                bins.bin_centers[label],
                rdf,
                label=f"Atom type {atom_labels_to_name[label]}",
                linewidth=2,
            )

        # Customize plot
        plt.xlabel("Distance (Å)", fontsize=12)
        plt.ylabel("Observations", fontsize=12)
        plt.title("Distance Observation by Atom Type", fontsize=14)
        plt.legend(frameon=True)
        plt.grid(True, alpha=0.3)

        # Set axis limits and ticks
        plt.xlim(left=0)
        plt.ylim(bottom=0)

        # Adjust layout to prevent label clipping
        plt.tight_layout()

        # Save figure with high DPI
        plt.savefig(
            f"data_out/binned_distance/distance_all_types.pdf",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()

        # Also save individual plots
        for label, rdf in binned_distances.rdf_bins.items():
            plt.figure(figsize=(8, 5))
            plt.plot(bins.bin_centers[label], rdf, linewidth=2, color="navy")

            plt.xlabel("Distance (Å)", fontsize=12)
            plt.ylabel("Observations", fontsize=12)
            plt.title(
                f"Distance Observation for Atom Type {atom_labels_to_name[label]}",
                fontsize=14,
            )
            plt.grid(True, alpha=0.3)
            plt.xlim(left=0)
            plt.ylim(bottom=0)

            plt.tight_layout()
            plt.savefig(
                f"data_out/binned_distance/binned_distance_{atom_labels_to_name[label]}.pdf",
                dpi=300,
                bbox_inches="tight",
            )
            plt.close()


def generate_rdf_plots(sd, bins, rdf_bins, atom_labels_to_name):

    if sd.rank == 0:

        pmmoto.io.io_utils.check_file_path("data_out/rdf/")
        # Set style for publication-quality figures
        # plt.style.use("seaborn")

        # Create figure with specific size
        plt.figure(figsize=(10, 6))

        # Plot RDFs for each atom type
        for label, rdf in rdf_bins.items():
            plt.plot(
                bins.bin_centers[label],
                rdf,
                label=f"Atom type {atom_labels_to_name[label]}",
                linewidth=2,
            )

        # Customize plot
        plt.xlabel("Distance (Å)", fontsize=12)
        plt.ylabel("Radial Distribution Function g(r)", fontsize=12)
        plt.title("Radial Distribution Functions by Atom Type", fontsize=14)
        plt.legend(frameon=True)
        plt.grid(True, alpha=0.3)

        # Set axis limits and ticks
        plt.xlim(left=0)
        plt.ylim(bottom=0)

        # Adjust layout to prevent label clipping
        plt.tight_layout()

        # Save figure with high DPI
        plt.savefig(f"data_out/rdf/rdf_all_types.pdf", dpi=300, bbox_inches="tight")
        plt.close()

        # Also save individual plots
        for label, rdf in rdf_bins.items():
            plt.figure(figsize=(8, 5))
            plt.plot(bins.bin_centers[label], rdf, linewidth=2, color="navy")

            plt.xlabel("Distance (Å)", fontsize=12)
            plt.ylabel("g(r)", fontsize=12)
            plt.title(f"RDF for Atom Type {atom_labels_to_name[label]}", fontsize=14)
            plt.grid(True, alpha=0.3)
            plt.xlim(left=0)
            plt.ylim(bottom=0)

            plt.tight_layout()
            plt.savefig(
                f"data_out/rdf/generate_rdf_{atom_labels_to_name[label]}.pdf",
                dpi=300,
                bbox_inches="tight",
            )
            plt.close()


def initialize_domain():
    """ """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Full domain with reservoirs
    box = [
        [0.0, 1.7576827931001799e02],
        [0.0, 1.7576827931001799e02],
        [-287.0, 237.0],
    ]

    # ignore reservoirs
    # Ignore water "reservoirs"
    box[2] = [-150, 175]

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
    Make sure the files match
    """

    len_1 = len(files_1)
    len_2 = len(files_2)

    files_1.sort()
    files_2.sort()

    if len_1 != len_2:
        raise ValueError("File lengths dont match")

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
            atom_name = atom_labels_to_name[label]
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

    if bridges:
        membrane_files = glob.glob(
            "/ocean/projects/cts200024p/rvickers/RV-P3/64x/03_90/perm_v2/membrane/*"
        )
        water_files = glob.glob(
            "/ocean/projects/cts200024p/rvickers/RV-P3/64x/03_90/perm_v2/pressure/*"
        )

        water_files.remove(
            "/ocean/projects/cts200024p/rvickers/RV-P3/64x/03_90/perm_v2/pressure/pressuredata.100000000.gz"
        )
        water_files.remove(
            "/ocean/projects/cts200024p/rvickers/RV-P3/64x/03_90/perm_v2/pressure/pressuredata.110005000.gz"
        )

        water_files.remove(
            "/ocean/projects/cts200024p/rvickers/RV-P3/64x/03_90/perm_v2/pressure/pressuredata.99995000.gz"
        )

    else:
        membrane_files = glob.glob("data/membrane_data/*")
        water_files = glob.glob("data/water_data/*")

    check_files(membrane_files, water_files)

    num_bins = 280

    max_radius = 2.8

    radii = {}
    for label in atom_labels_to_name:
        radii[label] = max_radius

    bins = pmmoto.domain_generation.rdf.generate_bins(radii, num_bins)

    sim_time = time.time()

    for n_file, (membrane_file, water_file) in enumerate(
            zip(membrane_files[2960:], water_files[2960:])
    ):

        # bins = pmmoto.domain_generation.rdf.generate_bins(radii, num_bins)

        if sd.rank == 0:
            iter_time = time.time()
            # print(f"Processed file {n_file} out of {len(membrane_files)}")

        membrane_positions, membrane_atom_type, _, timestep = (
            pmmoto.io.data_read.read_lammps_atoms(membrane_file, atom_id_charge_map)
        )

        water_positions, water_atom_type, _, timestep_2 = (
            pmmoto.io.data_read.read_lammps_atoms(water_file)
        )

        membrane_radii = gen_radii(membrane_atom_type, max_radius)
        water_radii = gen_radii(water_atom_type, max_radius)

        membrane = pmmoto.particles.initialize_atoms(
            sd,
            membrane_positions,
            membrane_radii,
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

        distance_bins = pmmoto.domain_generation.rdf.bin_distances(
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
            rdf = pmmoto.domain_generation.rdf.generate_rdf(distance_bins, bins)
            save_rdf(sd, bins, distance_bins, atom_labels_to_name, n_file)

    # Final save and generate plots
    rdf = pmmoto.domain_generation.rdf.generate_rdf(distance_bins, bins)
    generate_binned_distance_plots(sd, bins, distance_bins, atom_labels_to_name)
    generate_rdf_plots(sd, bins, rdf, atom_labels_to_name)
    save_rdf(sd, bins, distance_bins, atom_labels_to_name, n_file)


if __name__ == "__main__":
    bridges = True
    generate_rdf(bridges)
