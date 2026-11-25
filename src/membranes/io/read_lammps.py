import os
import glob
import re
import pmmoto
import numpy as np

__all__ = ["safe_label", "read_lammps_file"]

# Columns of the Atoms section (before any image flags), per LAMMPS atom_style.
# A .data file is structured output, so its columns are fixed by the atom_style
# and are independent of config.sim.data_columns (which describes dump files).
ATOM_STYLE_COLUMNS = {
    "atomic": ["id", "type", "x", "y", "z"],
    "charge": ["id", "type", "charge", "x", "y", "z"],
    "molecular": ["id", "mol", "type", "x", "y", "z"],
    "bond": ["id", "mol", "type", "x", "y", "z"],
    "angle": ["id", "mol", "type", "x", "y", "z"],
    "full": ["id", "mol", "type", "charge", "x", "y", "z"],
}

# Columns that should be stored as integers.
INT_COLUMNS = {"id", "mol", "type"}


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


def get_frame_files(config):
    search_pattern = os.path.join(config.input_path, config.file_pattern)
    frame_files = glob.glob(search_pattern)
    if not frame_files:
        raise FileNotFoundError(f"No frame files found matching: {search_pattern}")
    frame_files.sort(key=get_timestep)

    return frame_files


def _group_coordinates(atom_data):
    """Replace x/y/z columns with Nx3 arrays."""
    grouped = dict(atom_data)

    for name, keys in (("positions", ("x", "y", "z")),):
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
        initial_iteration=config.sim.initial_iteration,
        simulation_timestep=config.sim.simulation_timestep,
    )

    if extract_atom is not None:
        atom_data = _select_atoms_by_type(atom_data, extract_atom)

    atom_data = _group_coordinates(atom_data)

    return atom_data, domain, timestep


def _parse_lammps_data(file, atom_style=None):
    """Parse the box and Atoms section of a LAMMPS .data file.

    Returns (atom_data, domain): atom_data maps column name -> 1D array sorted
    by atom id; domain is a (3, 2) array [[xlo, xhi], [ylo, yhi], [zlo, zhi]].

    The columns come from the atom_style written by LAMMPS on the Atoms line
    (e.g. "Atoms # full"). If the file has no such comment, pass atom_style.
    Image flags and tilt factors are ignored.
    """
    with open(file, encoding="utf-8") as f:
        raw_lines = f.read().splitlines()
    lines = [line.split("#")[0].strip() for line in raw_lines]

    # Header: atom count and box bounds
    n_atoms = None
    bounds = {}
    atoms_start = None
    style_in_file = None
    for k, line in enumerate(lines):
        parts = line.split()
        if len(parts) == 2 and parts[1] == "atoms":
            n_atoms = int(parts[0])
        elif len(parts) == 4 and parts[2:] in (
            ["xlo", "xhi"],
            ["ylo", "yhi"],
            ["zlo", "zhi"],
        ):
            bounds[parts[2][0]] = (float(parts[0]), float(parts[1]))
        elif line == "Atoms":
            atoms_start = k + 1
            comment = raw_lines[k].partition("#")[2].split()
            style_in_file = comment[0] if comment else None
            break

    if atoms_start is None:
        raise ValueError(f"No 'Atoms' section found in {file}")
    if set(bounds) != {"x", "y", "z"}:
        raise ValueError(f"Could not read box bounds from {file}")

    # Explicit atom_style wins, then the style written in the file
    atom_style = atom_style or style_in_file
    if atom_style is None:
        raise ValueError(
            f"No atom_style found on the 'Atoms' line of {file}. "
            f"Set atom_style in the config. Choose from {sorted(ATOM_STYLE_COLUMNS)}"
        )
    if atom_style not in ATOM_STYLE_COLUMNS:
        raise ValueError(
            f"Unsupported atom_style '{atom_style}'. "
            f"Choose from {sorted(ATOM_STYLE_COLUMNS)}"
        )
    columns = ATOM_STYLE_COLUMNS[atom_style]
    ncols = len(columns)

    # Atoms section: skip blank lines after the header, stop at the next blank
    rows = []
    for line in lines[atoms_start:]:
        if not line:
            if rows:
                break
            continue
        parts = line.split()
        # ncols fields, or ncols + 3 when image flags (nx ny nz) are present
        if len(parts) not in (ncols, ncols + 3):
            raise ValueError(
                f"Atoms line has {len(parts)} fields but atom_style "
                f"'{atom_style}' expects {ncols} (or {ncols + 3} with image "
                f"flags): {line!r}"
            )
        rows.append(parts[:ncols])

    if n_atoms is not None and len(rows) != n_atoms:
        raise ValueError(f"Header says {n_atoms} atoms, read {len(rows)} in {file}")

    data = np.array(rows, dtype=float)
    data = data[np.argsort(data[:, 0], kind="stable")]  # id is always column 0

    atom_data = {
        name: data[:, i].astype(np.int64) if name in INT_COLUMNS else data[:, i]
        for i, name in enumerate(columns)
    }
    domain = np.array([bounds["x"], bounds["y"], bounds["z"]])

    return atom_data, domain


def _read_lammps_data_files(file, config, extract_atom):
    _check_input_path(file)

    # Optional override; normally the style is read from the file itself
    atom_style = getattr(config.sim, "atom_style", None)
    atom_data, domain = _parse_lammps_data(file, atom_style)

    try:
        timestep = get_timestep(file)
    except ValueError:  # e.g. "system.data" has no number in the name
        timestep = getattr(config.sim, "initial_iteration", 0)

    if extract_atom is not None:
        atom_data = _select_atoms_by_type(atom_data, extract_atom)

    atom_data = _group_coordinates(atom_data)

    return atom_data, domain, timestep


def read_lammps_file(file, config, extract_atom=None):
    if config.sim.input_type == "LAMMPS.dat":
        atom_data, domain, timestep = _read_lammps_data_files(
            file, config, extract_atom
        )
    else:
        atom_data, domain, timestep = _read_frame_files(file, config, extract_atom)

    return atom_data, domain, timestep
