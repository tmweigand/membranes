import numpy as np


class VickersMembrane:

    oxygen_type = 15

    box = [[0, 175.8], [0, 175.8], [-287.0, 237.0]]

    region_bounds = {
        "feed": (-287.0, -100.0),
        "membrane": (-100.0, 100.0),
        "permeate": (100.0, 237.0),
    }

    cutoffs = [3.0, 6.0]  # In  angstroms
    cutoff_step = 0.2
    initial_iteration = 71045000
    simulation_timestep = 2e-6
    data_columns = {
        "id": int,
        "mol": float,
        "type": np.uint8,
        "x": float,
        "y": float,
        "z": float,
    }
