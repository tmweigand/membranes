import numpy as np

from helpers.box import validate_box


def prepare_frame(frame_data):

    atom_ids = np.asarray(frame_data["ids"])
    positions = np.asarray(frame_data["positions"], dtype=np.float64)
    types = np.asarray(frame_data["types"])
    box = validate_box(
        frame_data["box"], frame_data.get("frame_label", "simulation frame")
    )

    if len(atom_ids) != len(positions):
        raise ValueError("Number of atom IDs does not match number of positions.")

    if len(atom_ids) != len(types):
        raise ValueError("Number of atom IDs does not match number of atom types.")

    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("Atom positions must have shape (N, 3).")

    return {"ids": atom_ids, "positions": positions, "types": types, "box": box}


def select_atoms_by_type(frame_data, atom_type):

    mask = frame_data["types"].astype(int) == int(atom_type)

    return {
        "ids": frame_data["ids"][mask],
        "positions": frame_data["positions"][mask],
        "types": frame_data["types"][mask],
    }
