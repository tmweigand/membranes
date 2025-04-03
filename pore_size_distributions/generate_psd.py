"""generates pore size distributions of membrane"""

import glob
from mpi4py import MPI
import numpy as np
import pmmoto
import os
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

    sd = initialize_domain()

    membrane_files = glob.glob("data/membrane_data/*")
    membrane_file = membrane_files[0]

    membrane_positions, membrane_atom_type, _, _ = (
        pmmoto.io.data_read.read_lammps_atoms(membrane_file)
    )

    unique_atom_types = np.unique(membrane_atom_type)
    mem_radii = {}
    for type in unique_atom_types:
        element = atom_type_element_map[type]
        mem_radii[type] = list(
            pmmoto.particles.uff_radius(atom_names=element).values()
        )[0]

    # membrane = pmmoto.particles.initialize_atoms(
    #     sd,
    #     membrane_positions,
    #     mem_radii,
    #     membrane_atom_type,
    #     by_type=False,
    #     add_periodic=True,
    #     set_own=False,
    # )

    # atom_locations = membrane.return_coordinates()
    # atom_types =

    pm = pmmoto.domain_generation.gen_pm_atom_domain(
        sd, membrane_positions, mem_radii, membrane_atom_type
    )

    img = pmmoto.filters.porosimetry.pore_size_distribution(
        sd, pm, 0.5866666666666667, inlet=False, plot="pdf"
    )

    pmmoto.io.output.save_img_data_parallel(
        "data_out/psd_images",
        sd,
        pm.img,
        additional_img={"psd": img},
    )


generate_psd(None)
