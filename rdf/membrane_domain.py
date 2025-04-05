"""membrane_domain.py"""

import glob
import numpy as np
from mpi4py import MPI
import pmmoto
import rdf_helpers
import time
from enum import Enum

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
proc_size = comm.Get_size()


def initialize_domain(voxels):
    """
    Initialize the membrane domain
    """

    if proc_size == 8:
        subdomains = (2, 2, 2)
    elif proc_size == 27:
        subdomains = (3, 3, 3)
    elif proc_size == 64:
        subdomains = (4, 4, 4)
    elif proc_size == 125:
        subdomains = (5, 5, 5)
    elif proc_size == 216:
        subdomains = (6, 6, 6)
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
        voxels=voxels,
        box=box,
        rank=rank,
        subdomains=subdomains,
        boundary_types=((2, 2), (2, 2), (0, 0)),
        verlet_domains=(20, 20, 20),
        inlet=((0, 0), (0, 0), (1, 0)),
        outlet=((0, 0), (0, 0), (0, 1)),
    )

    return sd


def generate_bounded_rdf():
    """
    Generate radii
    """
    atom_folder = "rdf/bridges_results/bins_extended/"
    atom_map, rdf = pmmoto.io.data_read.read_binned_distances_rdf(atom_folder)

    bounded_rdf = {}
    for _id, _rdf in rdf.items():
        bounded_rdf[_id] = pmmoto.domain_generation.rdf.Bounded_RDF.from_rdf(
            _rdf, 1.0e-3
        )

    return bounded_rdf


def determine_radii(bounded_rdf, pmf_value):
    """
    Collect the radii given a pmf cutoff
    """
    radii = {}
    for atom_id, _rdf in bounded_rdf.items():
        radii[atom_id] = _rdf.interpolate_radius_from_pmf(pmf_value)

    return radii


def generate_membrane_domain(pmf_value, subdomain, membrane_file):
    """
    Test for generating a radial distribution function from LAMMPS data
    """
    bounded_rdf = generate_bounded_rdf()

    radii = determine_radii(bounded_rdf, pmf_value)

    pmmoto.logger.info("Generating Porous Media.")
    pm = pmmoto.domain_generation.gen_pm_atom_file(
        subdomain=subdomain, lammps_file=membrane_file, atom_radii=radii, kd=False
    )

    pmmoto.logger.info("Connecting Components.")
    cc, _ = pmmoto.filters.connected_components.connect_components(
        img=pm.img, subdomain=subdomain
    )

    connections = pmmoto.filters.connected_components.inlet_outlet_connections(
        subdomain=subdomain, labeled_img=cc
    )

    if connections:
        pmmoto.logger.info(f"Connections found with {pmf_value}")
        return True
    else:
        pmmoto.logger.info(f"No Connections found with {pmf_value}")
        return False


def profile_bridges():
    """
    Get some information for scaling
    """
    membrane_files, _ = rdf_helpers.get_bridges_files()
    membrane_file = membrane_files[0]
    voxels = np.arange(3000, 6500, 500)
    for voxel in voxels:
        if rank == 0:
            print(
                f"Starting... {(voxel,voxel,voxel)}",
                flush=True,
            )
        start_time = time.time()
        subdomain = initialize_domain(voxel)
        generate_membrane_domain(15, subdomain, membrane_file)
        if rank == 0:
            print(f"Elapsed Time for {voxel} is {time.time()-start_time}", flush=True)
            print()
            print()


class Status(Enum):
    SUCCESS = "Success"
    LOWER_BOUND = "Lower bound is connected"
    UPPER_BOUND = "Upper bound is not connected"
    UNCONVERGED = "Value did not converge within maximum iterations"


def bisection_method(
    lower_bound,
    upper_bound,
    subdomain,
    membrane_file,
    tolerance=1e-4,
    max_iterations=100,
):

    # Ensure the initial bounds yield a connected patrh and a not connected path
    if not generate_membrane_domain(upper_bound, subdomain, membrane_file):
        pmmoto.logger.error("Upper Bound must have a Connection")
        return Status.UPPER_BOUND, None
    if generate_membrane_domain(lower_bound, subdomain, membrane_file):
        pmmoto.logger.error("Lower Bound must NOT have a Connection")
        return Status.LOWER_BOUND, None

    # Start the bisection loop
    iteration = 0

    while iteration < max_iterations:
        mid_bound = (lower_bound + upper_bound) / 2
        result = generate_membrane_domain(mid_bound, subdomain, membrane_file)

        change = np.abs(mid_bound - lower_bound)
        pmmoto.logger.info(f"Change of {change}")

        # If mid_point produces a True result, stop the bisection
        if result and change < tolerance:
            pmmoto.logger.info("Connected path with %e" % mid_bound)
            return Status.SUCCESS, mid_bound

        # Adjust the bounds based on the result at mid_bound
        if result:
            upper_bound = mid_bound
        else:
            lower_bound = mid_bound

        iteration += 1

    # If we exit the loop, we have failed to find a valid bound within the tolerance
    pmmoto.logger.error(
        f"Failed to find a valid bound within {max_iterations} iterations"
    )
    return Status.UNCONVERGED, None


# Function to write result and status to a file
def write_to_file(filename, file_number, processed_file, sim_status, result):
    with open(filename, "a") as f:
        f.write(
            f"{file_number}, File: {processed_file}, Status: {sim_status.name}, Result: {result}\n"
        )


if __name__ == "__main__":

    # Open the file for writing, clear previous content if needed
    if rank == 0:
        file_name = "connected_paths.out"
        with open(file_name, "w") as f:
            f.write("Connection Results\n")

    bridges = True

    # Grab Files
    if bridges:
        voxels_in = (3520, 3520, 4000)
        membrane_files, _ = rdf_helpers.get_bridges_files()
    else:
        voxels_in = (800, 800, 800)
        membrane_files = glob.glob("data/membrane_data/*")

    # Bounds for guesses.
    upper_pmf_data = 17.315
    lower_bound = 5

    sd = initialize_domain(voxels_in)
    for n_file, membrane_file in enumerate(membrane_files):
        status, pmf = bisection_method(lower_bound, upper_pmf_data, sd, membrane_file)
        if rank == 0:
            write_to_file(file_name, n_file, membrane_file, status, pmf)
