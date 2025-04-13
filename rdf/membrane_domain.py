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

logger = pmmoto.logger

import vtk
from vtk.util import numpy_support

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
        [-35, 65],
        # [-100, 100],
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


def determine_uff_radii():
    """
    Collect the radii given a pmf cutoff
    """
    atom_folder = "rdf/bridges_results/bins_extended/"
    atom_map, _ = pmmoto.io.data_read.read_binned_distances_rdf(atom_folder)
    radii = {}
    for atom_id, atom_data in atom_map.items():
        radii[atom_id] = (
            list(pmmoto.particles.uff_radius(atom_names=atom_data["element"]).values())[
                0
            ]
            + 1.4
        )

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
        kd=False,
        add_periodic=True,
    )

    # pm_morph = pmmoto.filters.morphological_operators.dilate(subdomain, pm.img, 1.4)

    cc, _ = pmmoto.filters.connected_components.connect_components(
        img=pm.img, subdomain=subdomain
    )

    connections = pmmoto.filters.connected_components.inlet_outlet_connections(
        subdomain=subdomain, labeled_img=cc
    )

    if connections:
        logger.info(f"Connections found! {connections}")

    else:
        logger.info(f"No Connections found.")
        return

    connected = np.where(cc == connections[0], 1, 0)

    _morph = pmmoto.filters.morphological_operators.dilate(subdomain, connected, 1.4)

    # _edt = pmmoto.filters.distance.edt(_morph.astype(np.uint8), subdomain)

    # create_surface("data_out/pm_morph_uff", pm_morph, subdomain)

    create_surface("data_out/connected_morph_check", _morph, subdomain)

    # return _morph

    # pmmoto.io.output.save_img_data_parallel(
    #     "data_out/membrane_domain",
    #     subdomain,
    #     pm_morph,  # additional_img={"edt": _edt}
    # )


def compare_radii(pmf_value, subdomain, membrane_file):
    """
    Generate plots for comparing approaches
    """
    bounded_rdf = generate_bounded_rdf()

    rdf_radii = determine_radii(bounded_rdf, pmf_value)
    uff_radii = determine_uff_radii()

    rdf_pm = pmmoto.domain_generation.gen_pm_atom_file(
        subdomain=subdomain,
        lammps_file=membrane_file,
        atom_radii=rdf_radii,
        type_map=atom_id_charge_map,
        kd=False,
        add_periodic=True,
    )

    rdf_pm.img = pmmoto.filters.morphological_operators.dilate(
        subdomain, rdf_pm.img, 1.4
    )

    uff_pm = pmmoto.domain_generation.gen_pm_atom_file(
        subdomain=subdomain,
        lammps_file=membrane_file,
        atom_radii=uff_radii,
        type_map=atom_id_charge_map,
        kd=False,
        add_periodic=True,
    )

    uff_pm.img = pmmoto.filters.morphological_operators.dilate(
        subdomain, uff_pm.img, 1.4
    )

    mask = np.where((uff_pm.img == 0) & (rdf_pm.img == 1), 2, uff_pm.img)

    pmmoto.io.output.save_img_data_parallel(
        "data_out/membrane_domain_compare",
        subdomain,
        rdf_pm.img,
        additional_img={"uff": uff_pm.img, "mask": mask},
    )


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


def save_water_locations(subdomain, water_file, img=None):
    """
    Save the locations of water molecules
    """
    water_positions, water_atom_type, _, _ = pmmoto.io.data_read.read_lammps_atoms(
        water_file
    )

    water_radii = {15: 1.4, 16: 1.4}

    water = pmmoto.particles.initialize_atoms(
        subdomain,
        water_positions,
        water_radii,
        water_atom_type,
        by_type=True,
        trim_within=True,
    )

    hydrogen = water.return_list(16)
    oxygen = water.return_list(15)

    save_based_on_img("data_out/water_oxygen", subdomain, img, oxygen)
    save_based_on_img("data_out/water_hydrogen", subdomain, img, hydrogen)


def save_based_on_img(file_name, subdomain, img, particle):

    particle_array = particle.return_np_array()

    particle_out = np.zeros_like(particle_array)
    particle_out[:, 0:3] = particle_array[:, 0:3]

    if img is not None:
        for n, atom in enumerate(particle_array):
            index = subdomain.get_img_index(atom[0:3])
            if index and img[index] == 1:
                particle_out[n, 3] = 1

    particle_out = particle_out[particle_out[:, 3] == 1]

    pmmoto.io.output.save_particle_data(file_name, subdomain, particle_out)


if __name__ == "__main__":

    bridges = False

    # Grab Files
    if bridges:
        voxels_in = (3520, 3520, 4000)
        membrane_files, _ = rdf_helpers.get_bridges_files()
    else:
        voxels_in = (800, 800, 800)
        membrane_files = glob.glob("data/membrane_data/membranedata.100020000")
        water_files = glob.glob("data/water_data/pressuredata.100020000")

    sd = initialize_domain(voxels_in)
    _water_file = water_files[0]
    _membrane_file = membrane_files[0]

    upper_pmf_data = 3.507219619750977

    # Determine conncetions and PLots with Water Locations
    generate_membrane_domain(upper_pmf_data, sd, _membrane_file)
    # connect_water = generate_membrane_domain(upper_pmf_data, sd, _membrane_file)
    # save_water_locations(sd, _water_file, connect_water)

    # Plot for comparing radii from RDF and UFF
    # compare_radii(upper_pmf_data, sd, _membrane_file)
