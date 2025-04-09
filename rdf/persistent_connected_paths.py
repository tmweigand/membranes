"""membrane_domain.py"""

import glob
import numpy as np
from mpi4py import MPI
import pmmoto
import rdf_helpers
import time
from enum import Enum
import vtk
from vtk.util import numpy_support

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


def write_pvd_file(file_name, num_ranks):

    pvd_filename = file_name + "/" + file_name.split("/")[-1] + ".pvd"

    with open(pvd_filename, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
        f.write("  <Collection>\n")
        for rank in range(num_ranks):
            file_proc = file_name.split("/")[-1] + "_" + str(rank) + ".vtp"
            f.write(f'    <DataSet timestep="0" part="{rank}" file="{file_proc}"/>\n')
        f.write("  </Collection>\n")
        f.write("</VTKFile>\n")


def create_surface(file_name, data, subdomain):

    threshold = 0.5  # Isosurface level

    # Convert NumPy array to VTK image data
    vtk_data = vtk.vtkImageData()
    vtk_data.SetDimensions(data.shape)
    vtk_data.SetOrigin(subdomain.get_origin())
    vtk_data.SetSpacing(subdomain.domain.resolution)  # Adjust spacing as needed

    # Convert NumPy array to VTK-compatible format
    flat_data = data.ravel(order="F")  # VTK expects Fortran-ordering
    vtk_array = numpy_support.numpy_to_vtk(
        flat_data, deep=True, array_type=vtk.VTK_FLOAT
    )
    vtk_data.GetPointData().SetScalars(vtk_array)

    # Apply Marching Cubes algorithm
    mc = vtk.vtkMarchingCubes()
    mc.SetInputData(vtk_data)
    mc.SetValue(0, threshold)  # Set isosurface value
    mc.Update()

    # Apply Laplacian smoothing
    smoother = vtk.vtkSmoothPolyDataFilter()
    smoother.SetInputConnection(mc.GetOutputPort())
    smoother.SetNumberOfIterations(60)  # More iterations = smoother mesh
    smoother.SetRelaxationFactor(0.1)  # Lower values = conservative smoothing
    smoother.FeatureEdgeSmoothingOff()  # Turn off edge preservation for a softer result
    smoother.Update()

    pmmoto.io.io_utils.check_file_path(file_name)
    file_proc = (
        file_name + "/" + file_name.split("/")[-1] + "_" + str(subdomain.rank) + ".vtp"
    )

    # Save the mesh as a VTP (VTK PolyData) file
    vtp_writer = vtk.vtkXMLPolyDataWriter()
    vtp_writer.SetFileName(file_proc)
    vtp_writer.SetInputConnection(smoother.GetOutputPort())
    vtp_writer.Write()

    if subdomain.rank == 0:
        write_pvd_file(file_name=file_name, num_ranks=subdomain.domain.num_subdomains)


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

    pm = pmmoto.domain_generation.gen_pm_atom_file(
        subdomain=subdomain,
        lammps_file=membrane_file,
        atom_radii=radii,
        type_map=atom_id_charge_map,
        add_periodic=True,
        kd=False,
    )

    cc, _ = pmmoto.filters.connected_components.connect_components(
        img=pm.img, subdomain=subdomain
    )

    connections = pmmoto.filters.connected_components.inlet_outlet_connections(
        subdomain=subdomain, labeled_img=cc
    )

    connected = np.where(cc == connections[0], 1, 0)

    conected_morph = pmmoto.filters.morphological_operators.dilate(
        subdomain, connected, 1.4
    )

    return conected_morph


# Function to write result and status to a file
def write_to_file(filename, file_number, processed_file, result):
    with open(filename, "a") as f:
        f.write(f"{file_number}, Last File: {processed_file},  Result: {result}\n")


if __name__ == "__main__":

    bridges = True

    # Grab Files
    if bridges:
        voxels_in = (3520, 3520, 4000)
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

    for n_file, membrane_file in enumerate(membrane_files):
        connected_img = generate_membrane_domain(pmf, sd, membrane_file)

        if n_file == 0 and n_persist == 1:
            persistent_connected_img = connected_img
        else:
            persistent_connected_img = np.where(
                (persistent_connected_img == 1) & (connected_img == 1), 1, 0
            )

            cc, _ = pmmoto.filters.connected_components.connect_components(
                img=persistent_connected_img, subdomain=sd
            )

            connections = pmmoto.filters.connected_components.inlet_outlet_connections(
                subdomain=sd, labeled_img=cc
            )

            if connections:
                n_persist += 1
                logger.info("Persistence %i" % n_persist)
            else:
                logger.info("Length fo that Persistence %i" % n_persist)
                if sd.rank == 0:
                    write_to_file(file_name, n_file, membrane_file, n_persist)

                n_persist = 1
                persistent_connected_img = connected_img
