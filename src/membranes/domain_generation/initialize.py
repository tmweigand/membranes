"""initialize.py"""

from membranes.domain_generation.lammps_init import LAMMPSInitialize
from membranes.domain_generation.style_properties import StyleProperties
from membranes.domain_generation.system_properties import SystemProperties

__all__ = ["initialize_simulation"]


def initialize_simulation(in_dir, log_file):
    """Initialize simulation"""
    wrapper = LAMMPSInitialize(
        log_file=log_file,
        units="real",
        atom_style="full",
    )

    if in_dir == "rv":
        StyleProperties(
            lmp=wrapper.lmp,
            dihedral_style="charmm",
            improper_style="harmonic",
            special_bonds="charmm",
        )
    else:
        StyleProperties(lmp=wrapper.lmp)

    SystemProperties(
        wrapper.lmp, dimension=3, n_atom_types=30, boundary=("p", "p", "p")
    )
    return wrapper
