"""membrane_domain.py"""

import glob
import numpy as np
from mpi4py import MPI
import pmmoto
import rdf_helpers


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
proc_size = comm.Get_size()

logger = pmmoto.logger

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
    elif proc_size == 320:
        subdomains = (8, 8, 5)
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

    pm = pmmoto.domain_generation.gen_pm_atom_file(
        subdomain=subdomain,
        lammps_file=membrane_file,
        atom_radii=radii,
        type_map=atom_id_charge_map,
        add_periodic=True,
        kd=False,
    )

    cc = pmmoto.filters.connected_components.connect_components(
        img=pm.img, subdomain=subdomain, return_label_count=False
    ).astype(np.uint16)

    connections = pmmoto.filters.connected_components.inlet_outlet_connections(
        subdomain=subdomain, labeled_img=cc
    )

    connected = np.where(cc == connections[0], 1, 0).astype(np.uint8)

    conected_morph = pmmoto.filters.morphological_operators.dilate(
        subdomain, connected, 1.4
    )

    return conected_morph


# Function to write result and status to a file
def write_to_file(filename, file_number, processed_file, result):
    with open(filename, "a") as f:
        f.write(f"{file_number}, Last File: {processed_file},  Result: {result}\n")


def check_persistence(
    sd, connected_img, persistent_connected_img, n_persist, file_info=None
):
    """
    Check persistence of connected paths between membrane snapshots.

    Parameters
    ----------
    sd : SubDomain
        The subdomain object containing rank information
    connected_img : ndarray
        Current membrane connectivity image
    persistent_connected_img : ndarray
        Previous persistent connectivity image
    n_persist : int
        Current persistence count
    file_info : tuple, optional
        Tuple containing (file_number, membrane_file, output_filename)

    Returns
    -------
    tuple
        (new_persistent_img, new_n_persist)
    """
    if persistent_connected_img is None:
        return connected_img, n_persist

    # Check intersection of current and persistent paths
    persistent_connected_img = np.where(
        (persistent_connected_img == 1) & (connected_img == 1), 1, 0
    ).astype(np.uint8)

    # Find connected components
    cc = pmmoto.filters.connected_components.connect_components(
        img=persistent_connected_img, subdomain=sd, return_label_count=False
    ).astype(np.uint32)

    # Check inlet/outlet connections
    connections = pmmoto.filters.connected_components.inlet_outlet_connections(
        subdomain=sd, labeled_img=cc
    )

    if connections:
        # Path persists
        n_persist += 1
        logger.info("Persistence %i" % n_persist)
        return persistent_connected_img, n_persist
    else:
        # Path broken - log and reset
        logger.info("Length of that Persistence %i" % n_persist)
        if file_info and sd.rank == 0:
            n_file, membrane_file, file_name = file_info
            write_to_file(file_name, n_file, membrane_file, n_persist)
        return connected_img, 1


if __name__ == "__main__":

    bridges = True

    # Grab Files
    if bridges:
        # voxels_in = (3520, 3520, 4000)
        voxels_in = (3080, 3080, 4000)
        membrane_files, _ = rdf_helpers.get_bridges_files()
    else:
        voxels_in = (800, 800, 800)
        membrane_files = glob.glob("data/more_membrane_data/*")
        membrane_files.sort()

    # Open the file for writing, clear previous content if needed
    if rank == 0:
        file_name = "persistent_paths.out"
        with open(file_name, "w") as f:
            f.write("Persistence of Path - Time is Count \n")

    pmf = 5

    sd = initialize_domain(voxels_in)
    n_persist = 1
    persistent_connected_img = None

    for n_file, membrane_file in enumerate(membrane_files):
        connected_img = generate_membrane_domain(pmf, sd, membrane_file)

        file_info = (n_file, membrane_file, file_name) if rank == 0 else None
        persistent_connected_img, n_persist = check_persistence(
            sd, connected_img, persistent_connected_img, n_persist, file_info
        )
