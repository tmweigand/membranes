import pmmoto
import numpy as np
from ..regions import make_region_bounds


class XuMembrane:

    oxygen_type = 10

    cutoffs = np.linspace(3.0, 4.0, 10)

    initial_iteration = 100005000
    simulation_timestep = 1e-6
    data_columns = {
        "id": int,
        "mol": float,
        "type": np.uint8,
        "x": float,
        "y": float,
        "z": float,
    }

    input_type = "LAMMPS.dat"


class CVFF_9_SS(XuMembrane):
    """
    Xu et al 2026 simulations results

    CFVV forcefield
    9 nm long domain
    SS self-supported restriction


    Box dimension in file are misleading.
    The (min,max) of all oxygen atom positions is:
        (0.00314,48.03018)
        (0.006,48.04603)
        (65.80597,247.00644)

    """

    box = [[0, 48.047298], [0, 48.047298], [65.80597, 247.00644]]

    membrane_min = 100  # Guess
    membrane_max = 200  # Guess

    eps = 0  # boundary cushion

    membrane_bounds = make_region_bounds(
        membrane_min, membrane_max, n_bins=2, name="mem"
    )

    region_bounds = {
        "feed": (box[2][0] + eps, membrane_min - eps),
        **membrane_bounds,
        "permeate": (membrane_max + eps, box[2][1] - eps),
    }


class CVFF_20_SS(XuMembrane):
    """
    Xu et al 2026 simulations results

    CFVV forcefield
    20 nm long domain
    SS self-supported restriction

    Box dimension in file are misleading.
    The (min,max) of all oxygen atom positions is:
        (0.00062,48.10156)
        (0.00062,48.10157)
        (69.37631,394.98597)

    """

    box = [[0, 48.102001], [0, 48.102001], [69.37631, 394.98597]]

    membrane_min = 200  # Guess
    membrane_max = 320  # Guess

    eps = 0  # boundary cushion

    membrane_bounds = make_region_bounds(
        membrane_min, membrane_max, n_bins=2, name="mem"
    )

    region_bounds = {
        "feed": (box[2][0] + eps, membrane_min - eps),
        **membrane_bounds,
        "permeate": (membrane_max + eps, box[2][1] - eps),
    }


class GAFF_9_SS(XuMembrane):
    """
    Xu et al 2026 simulations results

    GAFF forcefield
    9 nm long domain
    SS self-supported restriction

    Box dimension in file are misleading.
    The (min,max) of all oxygen atom positions is:
        (0.00585,47.35714)
        (0.30286,47.9865)
        (125.30069,203.47569)
    """

    box = [[0, 48.047298431396484], [0, 48.047298431396484], [125.30069, 203.47569]]

    membrane_min = 140  # Guess
    membrane_max = 180  # Guess

    eps = 0  # boundary cushion

    membrane_bounds = make_region_bounds(
        membrane_min, membrane_max, n_bins=2, name="mem"
    )

    region_bounds = {
        "feed": (box[2][0] + eps, membrane_min - eps),
        **membrane_bounds,
        "permeate": (membrane_max + eps, box[2][1] - eps),
    }


class OPLS_9_SS(XuMembrane):
    """
    Xu et al 2026 simulations results

    OPLS forcefield
    9 nm long domain
    SS self-supported restriction

    Box dimension in file are misleading.
    The (min,max) of all oxygen atom positions is:
        (0.53371,48.02555)
        (0.32448,48.03909)
        (118.62378,201.37071)
    """

    box = [[0, 48.0469970703125], [0, 48.0469970703125], [118.62378, 201.37071]]

    membrane_min = 130  # Guess
    membrane_max = 180  # Guess

    eps = 0  # boundary cushion

    membrane_bounds = make_region_bounds(
        membrane_min, membrane_max, n_bins=2, name="mem"
    )

    region_bounds = {
        "feed": (box[2][0] + eps, membrane_min - eps),
        **membrane_bounds,
        "permeate": (membrane_max + eps, box[2][1] - eps),
    }
