import numpy as np


class VickersMembrane:

    input_mode = "frame_files"

    oxygen_type = 15

    box = [[0, 100], [0, 100], [0, 100]]

    region_bounds = {
        "feed": (0, 0),
        "membrane": (0, 0),
        "permeate": (0, 0),
    }

    cutoffs = [3.0, 8.0]  # In  angstroms
    cutoff_step = 0.2

    data_columns = {
        "id": int,
        "mol": float,
        "type": np.uint8,
        # "mass": float,
        # "q": float,
        # "x": float,
        # "y": float,
        # "z": float,
        "xu": float,
        "yu": float,
        "zu": float,
        # "vx": float,
        # "vy": float,
        # "vz": float,
        # "v_peratompress": float,
        # "c_peratomvol[1]": float,
    }
