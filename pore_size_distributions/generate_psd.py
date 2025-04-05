"""generates pore size distributions of membrane"""

import glob
from mpi4py import MPI
import numpy as np
import pmmoto
import os
import csv
import psd_helpers
import time


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
        [-100, 100],  # Truncated z-axis based on density profile results
    ]

    sd = pmmoto.initialize(
        voxels=(300, 300, 300),
        box=box,
        rank=rank,
        subdomains=(2, 2, 2),
        boundary_types=((2, 2), (2, 2), (2, 2)),
        verlet_domains=(20, 20, 20),
    )

    return sd


atom_type_element_map = {
    1: "C",
    3: "C",
    4: "O",
    5: "H",
    7: "N",
    8: "H",
    12: "O",
    14: "H",
}


def generate_psd(bridges):
    """
    Generate pore size density plots using PMMoTo.
    """

    start_time = time.time()
    sd = initialize_domain()

    mem_radii = {}
    for atom_label, element in atom_type_element_map.items():
        mem_radii[atom_label] = list(
            pmmoto.particles.uff_radius(atom_names=element).values()
        )[0]

    pore_size_radii = np.linspace(1, 10, 50)

    psd_counts = {}
    for radius in pore_size_radii:
        psd_counts[radius] = 0

    if bridges:
        membrane_files, _ = psd_helpers.get_bridges_files()

    else:
        membrane_files = glob.glob("data/membrane_data/*")

    # Just take [0] entry in membrane files to analyze one specific file
    for membrane_file in [membrane_files[0]]:

        membrane_positions, membrane_atom_type, _, _ = (
            pmmoto.io.data_read.read_lammps_atoms(membrane_file)
        )

        if sd.rank == 0:
            print("File Read")

        pm = pmmoto.domain_generation.gen_pm_atom_domain(
            sd, membrane_positions, mem_radii, membrane_atom_type
        )

        if sd.rank == 0:
            print("Domain Generated")

        psd = pmmoto.filters.porosimetry.pore_size_distribution(
            sd, pm, pore_size_radii, inlet=False
        )

        counts = pmmoto.core.utils.bin_image(sd, psd)

        for radius, _count in counts.items():
            if radius == 0:
                continue
            psd_counts[radius] += _count

        if sd.rank == 0:
            print(psd_counts)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Elapsed Simulation Time: {elapsed_time:.2f}")

    if sd.rank == 0:
        # Save psd_counts to csv output
        z_voxels = 300
        with open(
            f"data_out/psd/psd_counts_{z_voxels}.csv", "w", encoding="utf-8"
        ) as csvfile:
            for key, value in psd_counts.items():
                csvfile.write(f"{key}, {value}\n")

        # pmmoto.filters.porosimetry.plot_pore_size_distribution(
        #     "data_out/psd",
        #     psd_counts,
        #     plot_type="cdf",
        # )

    # pmmoto.io.output.save_img_data_parallel(
    #     "data_out/psd_images_small",
    #     sd,
    #     pm.img,
    #     additional_img={"psd": psd},
    # )


if __name__ == "__main__":
    bridges_run = False
    generate_psd(bridges_run)
