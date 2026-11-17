import numpy as np


class TestSim:

    input_mode = "frame_files"

    oxygen_type = 1

    box = [[0, 20], [0, 20], [0, 20]]

    region_bounds = {
        "feed": (0, 10),
        "membrane": (10, 15),
        "permeate": (15, 20),
    }

    cutoffs = [3.0, 8.0]  # In  angstroms
    cutoff_step = 0.2

    data_columns = {
        "id": int,
        "type": np.uint8,
        "x": float,
        "y": float,
        "z": float,
    }
