"""oligomer_example.py"""

import pmmoto
import numpy as np


def oligomer():
    """
    Generate nice plots to deomonstrate this work on an oligomer
    """

    o_coords_files = "domain_in/oligomer_coordinates.txt"

    data = np.genfromtxt(o_coords_files, delimiter=",", dtype=str, skip_header=1)

    _atom_ids = data[:, 0]  # First column

    spheres = data[:, 1:].astype(float)  # Convert the rest to float
    # radii = data[:, 4].astype(float)

    box = ((-10, 30), (-10, 30), (10, 30))

    sd = pmmoto.initialize(voxels=(500, 500, 500), box=box)

    # count = 0
    # atom_radii = {}
    # for _id in _atom_ids:
    #     if _id[0] == "H":
    #         atom_radii[count] = 1.5
    #     else:
    #         atom_radii[count] = 2.75
    #     count += 1

    # count = 0
    # atom_ids = {}
    # for atom in _atom_ids:
    #     if atom not in atom_ids:
    #         atom_ids[count] = count
    #         count += 1

    pm = pmmoto.domain_generation.gen_pm_spheres_domain(sd, spheres)

    pores = pmmoto.filters.morphological_operators.erode(sd, pm.img, 2.8)
    pores2 = pmmoto.filters.morphological_operators.dilate(sd, pores, 2.8)

    pmmoto.io.output.save_img_data_parallel(
        "data_out/oligomer_van_der_waals",
        sd,
        pm.img,
        additional_img={"erode": pores, "dilate": pores2},
    )


if __name__ == "__main__":
    oligomer()
