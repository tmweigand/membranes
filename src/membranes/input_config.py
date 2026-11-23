import os
import glob
import re

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

        self.initial_iteration = self.sim.initial_iteration
        self.simulation_timestep = self.sim.simulation_timestep

        # Individual frame file pattern
        self.file_pattern = "*"

    # def get_input_files(self):
    #     """Return a list of files from the specified input directory

    #     TODO: Add sorting option that is simulation specific
    #     """

    #     return glob.glob(os.path.join(os.path.abspath(self.in_dir), "*"))

    def get_input_files(self):
        """Return numerically sorted files from the input directory."""

        files = glob.glob(os.path.join(os.path.abspath(self.in_dir), "*"))

        return sorted(
            files, key=lambda f: int(re.search(r"\d+", os.path.basename(f)).group())
        )
