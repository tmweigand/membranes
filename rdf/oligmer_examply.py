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

    coordinates = data[:, 1:].astype(float)  # Convert the rest to float

    box = ((-10, 30), (-10, 30), (10, 30))

    sd = pmmoto.initialize(voxels=(500, 500, 500), box=box)

    count = 0
    atom_radii = {}
    for _id in _atom_ids:
        if _id[0] == "H":
            atom_radii[count] = 1.5
        else:
            atom_radii[count] = 2.75
        count += 1

    count = 0
    atom_ids = {}
    for atom in _atom_ids:
        if atom not in atom_ids:
            atom_ids[count] = count
            count += 1

    # atoms = pmmoto.domain_generation.particles.initialize_atoms(
    #     sd, coordinates, atom_radii, atom_ids, by_type=False
    # )

    pm = pmmoto.domain_generation.gen_pm_atom_domain(
        sd, coordinates, atom_radii, atom_ids
    )

    pmmoto.io.output.save_img_data_parallel(
        "data_out/oligomer_equilibrium_relaxed", sd, pm.img
    )


if __name__ == "__main__":
    oligomer()
