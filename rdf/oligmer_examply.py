"""oligomer_example.py"""

import pmmoto
import numpy as np
import vtk
from vtk.util import numpy_support


import matplotlib.pyplot as plt


def create_surface(data):

    # data = np.random.rand(30, 30, 30).astype(np.float32)  # Ensure it's float32
    threshold = 0.5  # Isosurface level

    # Convert NumPy array to VTK image data
    vtk_data = vtk.vtkImageData()
    vtk_data.SetDimensions(data.shape)
    vtk_data.SetSpacing(1.0, 1.0, 1.0)  # Adjust spacing as needed

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

    # Save the mesh as an STL file
    stl_writer = vtk.vtkSTLWriter()
    stl_writer.SetFileName("output_mesh.stl")
    stl_writer.SetInputConnection(smoother.GetOutputPort())
    stl_writer.Write()

    print("Marching Cubes mesh saved as 'output_mesh.stl'.")


def oligomer_uff():
    """
    Generate nice plots to deomonstrate this work on an oligomer
    """

    o_coords_files = "domain_in/oligomer_coordinates.txt"
    data = np.genfromtxt(o_coords_files, delimiter=",", dtype=str, skip_header=1)
    atom_names = data[:, 0]  # Atom names

    # Differentiate atom names but location in polyamide
    atom_labels = data[:, 1]

    uff_radii = pmmoto.particles.uff_radius(atom_names=atom_names)
    atom_types = pmmoto.particles.convert_atoms_elements_to_ids(atom_names)

    print(uff_radii)

    coordinates = data[:, 2:].astype(float)  # Convert the rest to float

    box = ((-10, 30), (-10, 30), (-10, 30))

    sd = pmmoto.initialize(voxels=(500, 500, 500), box=box)

    pm = pmmoto.domain_generation.gen_pm_atom_domain(
        subdomain=sd,
        atom_locations=coordinates,
        atom_radii=uff_radii,
        atom_types=atom_types,
    )

    pores = pmmoto.filters.morphological_operators.erode(sd, pm.img, 2.8 / 2)
    pores2 = pmmoto.filters.morphological_operators.dilate(sd, pores, 2.8 / 2)

    create_surface(pores2)

    pmmoto.io.output.save_img_data_parallel(
        "data_out/oligomer_van_der_waals",
        sd,
        pm.img,
        # additional_img={"erode": pores, "dilate": pores2},
    )


def oligomer_rdf():
    """
    Generate nice plots to deomonstrate this work on an oligomer
    """

    o_coords_files = "domain_in/oligomer_coordinates.txt"
    data = np.genfromtxt(o_coords_files, delimiter=",", dtype=str, skip_header=1)
    oligomer_elements = data[:, 0]  # Atom names

    # Differentiate atom names but location in polyamide
    atom_species_name = data[:, 1]
    unique_species = np.unique(atom_species_name)

    # atom_folder = "data_out/rdf_data/"
    atom_folder = "rdf/bridges_results/bins/"
    atom_map, rdf = pmmoto.io.data_read.read_binned_distances_rdf(atom_folder)

    bounded_rdf = {}
    for _id, _rdf in rdf.items():
        bounded_rdf[_id] = pmmoto.domain_generation.rdf.Bounded_RDF.from_rdf(
            _rdf, 1.0e-3
        )

    # # Convert atom_map to lookup based on species name
    # atom_by_label = {
    #     value["label"]: {"element": value["element"], "id": key}
    #     for key, value in atom_map.items()
    # }

    # rdf_radii = {}
    # for label in unique_species:
    #     if label == "H":
    #         id = 9  # BH1
    #     elif label == "O":
    #         id = 8  # O_CONH
    #     else:
    #         id = atom_by_label[label]["id"]

    #     rdf_radii[label] = rdf[label].interpolate_radius_from_pmf(5)

    # # atom_types = pmmoto.particles.convert_atoms_elements_to_ids(atom_names)

    # coordinates = data[:, 2:].astype(float)  # Convert the rest to float

    # box = ((-10, 30), (-10, 30), (-10, 30))

    # sd = pmmoto.initialize(voxels=(500, 500, 500), box=box)

    # # rdf_radii = {}
    # # for key, _rdf in rdf.items():
    # #     rdf_radii[_rdf.atom_id] = _rdf.r_from_G(5000)

    # # print(rdf_radii)

    # # atom_types = []
    # # for

    # # print(rdf_radii)

    # pm = pmmoto.domain_generation.gen_pm_atom_domain(
    #     subdomain=sd,
    #     atom_locations=coordinates,
    #     atom_radii=rdf_radii,
    #     atom_types=atom_types,
    # )

    # pores = pmmoto.filters.morphological_operators.erode(sd, pm.img, 2.8 / 2)
    # pores2 = pmmoto.filters.morphological_operators.dilate(sd, pores, 2.8 / 2)

    # create_surface(pores2)

    # pmmoto.io.output.save_img_data_parallel(
    #     "data_out/oligomer_rdf",
    #     sd,
    #     pm.img,
    #     # additional_img={"erode": pores, "dilate": pores2},
    # )


if __name__ == "__main__":
    # oligomer_uff()
    oligomer_rdf()
