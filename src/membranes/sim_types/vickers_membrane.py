import numpy as np
from ..regions import make_region_bounds


class VickersMembrane:

    input_type = None

    oxygen_type = 15

    box = [[0, 175.8], [0, 175.8], [-287.0, 237.0]]

    membrane_min = -100
    membrane_max = 100

    eps = 0  # boundary cushion

    membrane_bounds = make_region_bounds(
        membrane_min, membrane_max, n_bins=3, name="mem"
    )

    region_bounds = {
        "feed": (box[2][0] + eps, membrane_min - eps),
        **membrane_bounds,
        "permeate": (membrane_max + eps, box[2][1] - eps),
    }

    cutoffs = np.linspace(3.0, 4.0, 10)
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
