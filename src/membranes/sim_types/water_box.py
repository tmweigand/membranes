import numpy as np


class WaterBox:

    input_mode = "frame_files"

    oxygen_type = 1
    box = [[0, 100], [0, 100], [0, 100]]
    region_bounds = {"whole_box": (-np.inf, np.inf)}

    cutoffs = [1.0, 8.0]  # In  angstroms
    cutoff_step = 0.2
