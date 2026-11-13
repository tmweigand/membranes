import numpy as np


# Region definitions (Angstrom)
def get_region_bounds(simtype):

    if simtype == "watersim":
        return {"whole_box": (-np.inf, np.inf)}

    elif simtype == "membranesim":
        return {
            "feed": (70.0, 80.0),
            "membrane": (95.0, 170.0),
            "permeate": (195.0, 240.0),
        }

    else:
        raise ValueError(f"Unsupported simulation type: {simtype}")


initial_iteration = 71045000
simulation_timestep = 2e-6
unwrapped_coordinates = False

# Individual frame file pattern
file_pattern = "*"

# Existing cutoff increment: 0.2 Angstrom
cutoff_step_A = 0.2
