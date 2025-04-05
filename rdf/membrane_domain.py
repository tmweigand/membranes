"""membrane_domain.py"""

import glob
import numpy as np
from mpi4py import MPI
import pmmoto
import matplotlib.pyplot as plt
import rdf_helpers
import time


def initialize_domain(voxel):
    """
    Initialize the membrane domain
    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    proc_size = comm.Get_size()

    if proc_size == 8:
        subdomains = (2, 2, 2)
    elif proc_size == 27:
        subdomains = (3, 3, 3)
    elif proc_size == 64:
        subdomains = (4, 4, 4)
    elif proc_size == 125:
        subdomains = (5, 5, 5)
    else:
        print("oops")
        subdomains = (0, 0, 0)

    # Full domain with reservoirs
    box = [
        [0.0, 176],
        [0.0, 176],
        [-100, 100],  # Ignore water reservoirs
    ]

    sd = pmmoto.initialize(
        voxels=(voxel, voxel, voxel),
        box=box,
        rank=rank,
        subdomains=subdomains,
        boundary_types=((2, 2), (2, 2), (0, 0)),
        verlet_domains=(20, 20, 20),
        inlet=((0, 0), (0, 0), (1, 0)),
        outlet=((0, 0), (0, 0), (0, 1)),
    )

    return sd


def generate_membrane_domain(subdomain, membrane_file):
    """
    Test for generating a radial distribution function from LAMMPS data
    """

    # membrane_atom_label_file = "rdf/atom_map.txt"
    # atom_labels_to_name = pmmoto.io.data_read.read_atom_map(membrane_atom_label_file)

    atom_folder = "rdf/bridges_results/bins_extended/"
    atom_map, rdf = pmmoto.io.data_read.read_binned_distances_rdf(atom_folder)

    bounded_rdf = {}
    for _id, _rdf in rdf.items():
        bounded_rdf[_id] = pmmoto.domain_generation.rdf.Bounded_RDF.from_rdf(
            _rdf, 1.0e-3
        )

    radii = {}
    for atom_id, _rdf in bounded_rdf.items():
        radii[atom_id] = _rdf.interpolate_radius_from_pmf(10.0)

    pm = pmmoto.domain_generation.gen_pm_atom_file(
        subdomain=subdomain, lammps_file=membrane_file, atom_radii=radii
    )

    cc, _ = pmmoto.filters.connected_components.connect_components(
        img=pm.img, subdomain=subdomain
    )

    connections = pmmoto.filters.connected_components.inlet_outlet_connections(
        subdomain=subdomain, labeled_img=cc
    )

    if subdomain.rank == 0:
        print(f"Connections:{connections}")
        # print(f"Simulation Time: {time.time() - start_time}")

    # connected_path = np.where(cc == 34, 1, 0)

    # water_path = pmmoto.filters.morphological_operators.dilate(sd, connected_path, 1.4)

    # pmmoto.io.output.save_img_data_parallel(
    #     "data_out/membrane_domain_2",
    #     sd,
    #     pm.img,
    #     additional_img={"water_path": water_path},
    # )


def profile_bridges():
    """
    Get some information for scaling
    """
    membrane_files, _ = rdf_helpers.get_bridges_files()
    membrane_file = membrane_files[0]
    voxels = np.arange(1500, 3500, 500)
    for voxel in voxels:
        start_time = time.time()
        subdomain = initialize_domain(voxel)
        generate_membrane_domain(subdomain, membrane_file)
        print(f"Elapsed Time for {voxel} is {time.time()-start_time}")


if __name__ == "__main__":

    profile_bridges()
    # membrane_files = glob.glob("data/membrane_data/*")
    # membrane_file = membrane_files[0]
    # subdomain = initialize_domain(1200)
    # generate_membrane_domain(subdomain, membrane_file)
