"""test_lammps_wrapper.py"""

import pytest
import membranes


def test_import_lammps():
    try:
        from lammps import lammps
    except ImportError as e:
        pytest.fail(f"Failed to import lammps: {e}")


def test_lammps_instance():
    from lammps import lammps

    args = ["-log", "none"]
    lmp = lammps(cmdargs=args)
    lmp.close()


def test_lammps_init():
    lammps_wrapper = membranes.domain_generation.lammps_init.LAMMPSInitialize(
        log_file="test_log", units="real", atom_style="full"
    )

    lammps_wrapper.lmp.close()
