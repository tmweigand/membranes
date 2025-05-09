"""persistent_connected_paths.py"""

import glob
import sys
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
        [-35, 65],  # Central dense region
        # [-100, 100],  # Ignore water reservoirs
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


def max_distance(pmf_value, subdomain, membrane_file):
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

    edt = pmmoto.filters.distance.edt(pm.img, subdomain=subdomain)
    max_edt = pmmoto.core.utils.determine_maximum(edt)

    cc = pmmoto.filters.connected_components.connect_components(
        img=pm.img, subdomain=subdomain, return_label_count=False
    )

    connections = pmmoto.filters.connected_components.inlet_outlet_connections(
        subdomain=subdomain, labeled_img=cc
    )

    if connections == []:
        max_edt_connected = None
    else:
        connected = np.where(cc == connections[0], 1, 0).astype(np.uint8)
        edt = pmmoto.filters.distance.edt(connected, subdomain=subdomain)
        max_edt_connected = pmmoto.core.utils.determine_maximum(edt)

    return max_edt, max_edt_connected


# Function to write result and status to a file
def write_to_file(filename, file_number, processed_file, result, result2):
    with open(filename, "a") as f:
        f.write(
            f"{file_number}, Last File: {processed_file},  All: {result} Connected: {result2} \n"
        )


if __name__ == "__main__":

    bridges = True

    # Grab Files
    if bridges:
        # voxels_in = (3520, 3520, 4000)
        voxels_in = (3080, 3080, 3500)
        membrane_files, _ = rdf_helpers.get_bridges_files()
    else:
        voxels_in = (1200, 1200, 1200)
        membrane_files = glob.glob("data/more_membrane_data/*")
        membrane_files.sort()

    # Use command-line argument for pmf
    if len(sys.argv) < 2:
        raise ValueError("Need to specify the pmf")
    pmf = float(sys.argv[1])

    # pmf = 1.4  # average
    # pmf = 3.61  # max observed min G for connected path
    # pmf = 17.315  # Maximum G from data

    # Open the file for writing, clear previous content if needed
    if rank == 0:
        file_name = f"max_distance_{pmf}.out"
        with open(file_name, "w") as f:
            f.write("In Angstroms - From water center \n")

    sd = initialize_domain(voxels_in)
    for n_file, membrane_file in enumerate(membrane_files):
        max_edt, max_edt_connected = max_distance(
            pmf_value=pmf, subdomain=sd, membrane_file=membrane_file
        )
        if rank == 0:
            write_to_file(file_name, n_file, membrane_file, max_edt, max_edt_connected)
