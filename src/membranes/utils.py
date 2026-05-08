"""src/membranes/utils.py"""

import shutil

REQUIRED_TOOLS = ["antechamber", "parmchk2", "tleap", "obabel", "lmp"]


def check_dependencies():
    missing = [t for t in REQUIRED_TOOLS if shutil.which(t) is None]
    if missing:
        raise EnvironmentError(
            f"Missing required tools: {missing}\n"
            "Install with: conda install -c conda-forge ambertools openbabel lammps"
        )
