import os
import glob
import re

from .sim_types.water_box import WaterBox
from .sim_types.vickers_membrane import VickersMembrane
from .sim_types.xu_membrane import CVFF_9_SS, CVFF_20_SS, GAFF_9_SS, OPLS_9_SS
from .sim_types.test_sim import TestSim

SYSTEMS = {
    "water_box": WaterBox,
    "vickers_membrane": VickersMembrane,
    "xu_cvff_9_ss": CVFF_9_SS,
    "xu_cvff_20_ss": CVFF_20_SS,
    "xu_gaff_9_ss": GAFF_9_SS,
    "xu_opls_9_ss": OPLS_9_SS,
    "test_sim": TestSim,
}


class InputConfig:
    """Configuration class that store information for various simulations"""

    def __init__(self, in_dir, out_dir, sim_type):
        self.in_dir = in_dir
        self.out_dir = out_dir
        self.sim_type = sim_type
        self.sim = SYSTEMS[sim_type]()
        self.box = self.sim.box

    def get_input_files(self):
        """Return [file] if in_dir is a file, else the numerically sorted files in it."""
        path = os.path.abspath(self.in_dir)

        if os.path.isfile(path):
            return [path]
        if not os.path.isdir(path):
            raise FileNotFoundError(f"Input path not found: {path}")

        files = [f for f in glob.glob(os.path.join(path, "*")) if os.path.isfile(f)]

        def sort_key(f):
            numbers = re.findall(r"\d+", os.path.basename(f))
            return (0, int(numbers[-1]), f) if numbers else (1, 0, f)

        return sorted(files, key=sort_key)
