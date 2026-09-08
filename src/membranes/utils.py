"""src/membranes/utils.py"""

import shutil
import numpy as np

__all__ = ["check_dependencies", "random_int"]

REQUIRED_TOOLS = ["antechamber", "parmchk2", "tleap", "obabel", "lmp"]


def check_dependencies():
    missing = [t for t in REQUIRED_TOOLS if shutil.which(t) is None]
    if missing:
        raise EnvironmentError(
            f"Missing required tools: {missing}\n"
            "Install with: conda install -c conda-forge ambertools openbabel lammps"
        )


_rngs: dict[int, np.random.Generator] = {}


def random_int(seed: int | None) -> int:
    if seed is None:
        return int(np.random.randint(1, 32768))

    if seed not in _rngs:
        _rngs[seed] = np.random.default_rng(seed)

    return int(_rngs[seed].integers(1, 32768))
