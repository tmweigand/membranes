import numpy as np

# Input mode: "trajectory" or "frame_files"
input_mode = "frame_files"

# Individual frame files (PMMoTo)
use_helpers_get_files = False  ## Set to True to use the helper function to get frame files, False to specify path manually
machine = "mac"
file_type = "gzip"

frame_files_pattern = "/hickory/users/r/g/rghanta/water_data/Xu/CVFF_9_SS/*"

# Water oxygen atom TYPE
oxygen_type = 10

# Coordination cutoffs (nm)
cutoffs_nm = np.arange(30, 81, 2, dtype=np.int64) / 100.0

# Regions based on oxygen z coordinates (Angstrom)
region_bounds = {
    "feed": (70, 80),
    "membrane": (95.0, 170.0),
    "permeate": (195.0, 240.0),
}

# Box dimensions (Angstrom)
frame_files_box = np.array(
    [48.0473, 48.0473, 286.642, 90.0, 90.0, 90.0], dtype=np.float64
)

# PMMoTo reader parameters
initial_iteration = 71045000
simulation_timestep = 2e-6
unwrapped_coordinates = False

# Output directory
output_dir = "membrane_MD/Analysis/coordination_number/Xu/CVFF_9_SS/"
