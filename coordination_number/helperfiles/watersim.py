import numpy as np

# Input mode: "trajectory" or "frame_files"
input_mode = "trajectory"

# Trajectory file (MDAnalysis)
trajectory_file = (
    "/hickory/users/r/g/rghanta/water_data/tip3p_orig/100bar/traj_100bar.lammpstrj"
)

# Water oxygen atom TYPE
oxygen_type = 1

# Coordination cutoffs (nm)
cutoffs_nm = np.arange(30, 81, 2, dtype=np.int64) / 100.0

# Regions
region_bounds = {"whole_box": (-np.inf, np.inf)}

# Output directory
output_dir = "water_MD/Analysis/coordination_number/tip3p_orig/100bar/"
