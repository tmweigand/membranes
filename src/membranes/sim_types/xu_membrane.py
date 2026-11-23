import numpy as np


class XuMembrane:

    oxygen_type = 10

    box = [[0, 48.102], [0, 48.102], [0, 286.642]]

    region_bounds = {
        "feed": (66.0, 87.0),
        "membrane": (87.0, 190.0),
        "permeate": (190.0, 250.0),
    }

    cutoffs = [3.0, 8.0]  # In  angstroms
    cutoff_step = 0.2
    initial_iteration = 100005000
    simulation_timestep = 1e-6
    data_columns = {
        "id": int,
        "mol": float,
        "type": np.uint8,
        # "mass": float,
        "q": float,
        "x": float,
        "y": float,
        "z": float,
        "xu": float,
        "yu": float,
        "zu": float,
        # "vx": float,
        # "vy": float,
        # "vz": float,
        # "v_peratompress": float,
        # "c_peratomvol[1]": float,
    }
