"""oligomer_example.py"""

import pmmoto
import numpy as np
import vtk
from vtk.util import numpy_support


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


def oligomer():
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


if __name__ == "__main__":
    oligomer()
