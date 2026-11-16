import os
import glob
import re
import pmmoto
import numpy as np
import MDAnalysis as mda


from ..box import validate_box

__all__ = ["safe_label", "read_lammps_file"]


def safe_label(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("_")


def get_timestep(filename):
    numbers = re.findall(r"\d+", os.path.basename(filename))

    if not numbers:
        raise ValueError(f"Could not identify timestep in {filename}")

    return int(numbers[-1])


def _check_input_path(file):
    """Ensure file exists before attempting to read"""
    if not os.path.isfile(file):
        raise FileNotFoundError(f"File not found: {file}")


def read_trajectory_frames(file, config):

    universe = mda.Universe(file, format="LAMMPSDUMP")
    all_types = np.asarray(universe.atoms.types, dtype=int)

    for ts in universe.trajectory:
        timestep = ts.data.get("step", ts.frame)
        data = {
            "ids": universe.atoms.ids,
            "positions": universe.atoms.positions,
            "types": universe.atoms.types,
            "box": ts.dimensions,
        }

    return data


def get_frame_files(config):
    search_pattern = os.path.join(config.input_path, config.file_pattern)
    frame_files = glob.glob(search_pattern)
    if not frame_files:
        raise FileNotFoundError(f"No frame files found matching: {search_pattern}")
    frame_files.sort(key=get_timestep)

    return frame_files


def _group_coordinates(atom_data):
    """Replace x/y/z and xu/yu/zu columns with Nx3 arrays."""
    grouped = dict(atom_data)

    for name, keys in (
        ("coords", ("x", "y", "z")),
        ("coords_unwrapped", ("xu", "yu", "zu")),
    ):
        if all(k in atom_data for k in keys):
            grouped[name] = np.column_stack([grouped.pop(k) for k in keys])

    return grouped


def _select_atoms_by_type(atom_data, atom_type):
    """Return only the atoms whose type equals atom_type."""
    types = atom_data["type"] if "type" in atom_data else atom_data["types"]
    mask = np.asarray(types) == atom_type
    n_atoms = len(mask)

    return {
        key: (
            value[mask]
            if isinstance(value, np.ndarray) and value.shape[:1] == (n_atoms,)
            else value
        )
        for key, value in atom_data.items()
    }


def _read_frame_files(file, config, extract_atom):

    _check_input_path(file)

    atom_data, domain, timestep = pmmoto.io.data_read.read_lammps_atoms(
        input_file=file,
        data_columns=config.sim.data_columns,
        initial_iteration=config.initial_iteration,
        simulation_timestep=config.simulation_timestep,
    )

    if extract_atom is not None:
        atom_data = _select_atoms_by_type(atom_data, extract_atom)

    atom_data = _group_coordinates(atom_data)

    return atom_data, domain, timestep


def read_lammps_file(file, config, extract_atom=None):
    if config.sim.input_mode == "trajectory":
        atom_data = read_trajectory_frames(file, config)
    elif config.sim.input_mode == "frame_files":
        atom_data, domain, timestep = _read_frame_files(file, config, extract_atom)

    return atom_data, domain, timestep
