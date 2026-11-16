class XuMembrane:

    input_mode = "trajectory"

    oxygen_type = 15

    box = [[0, 100], [0, 100], [0, 100]]

    region_bounds = {
        "feed": (70.0, 80.0),
        "membrane": (95.0, 170.0),
        "permeate": (195.0, 240.0),
    }

    cutoffs = [3.0, 8.0]  # In  angstroms
    cutoff_step = 0.2
