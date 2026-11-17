import os
import glob

from .sim_types.water_box import WaterBox
from .sim_types.vickers_membrane import VickersMembrane
from .sim_types.xu_membrane import XuMembrane
from .sim_types.test_sim import TestSim

SIM_TYPES = {
    "water_box": WaterBox,
    "vickers_membrane": VickersMembrane,
    "xu_membrane": XuMembrane,
    "test_sim": TestSim,
}


class InputConfig:
    """Configuration class that store information for various simulations"""

    def __init__(self, in_dir, out_dir, sim_type):
        self.in_dir = in_dir
        self.out_dir = out_dir
        self.sim_type = sim_type
        self.sim = SIM_TYPES[sim_type]()

        self.box = self.sim.box

        if self.sim.cutoffs[0] <= 0 or self.sim.cutoffs[1] < self.sim.cutoffs[0]:
            raise ValueError("Invalid coordination cutoff range.")

        self.initial_iteration = 71045000
        self.simulation_timestep = 2e-6
        self.unwrapped_coordinates = False

        # Individual frame file pattern
        self.file_pattern = "*"

        # Sim Type Checks
        if self.sim.input_mode not in {"trajectory", "frame_files"}:
            raise ValueError("input_mode must be 'trajectory' or 'frame_files'.")

        if not self.sim.region_bounds:
            raise ValueError("region_bounds must contain at least one region.")

    def make_output_directory(self):
        """Make the output directory"""
        os.makedirs(self.out_dir, exist_ok=True)

    def get_input_files(self):
        """Return a list of files from the specified input directory

        TODO: Add sorting option that is simulation specific
        """
        return glob.glob(os.path.join(os.path.abspath(self.in_dir), "*"))
