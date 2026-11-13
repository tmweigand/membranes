import os

import numpy as np
import MDAnalysis as mda


def read_trajectory_frames(config, stats, analysis_function, process_frame):
    if not os.path.isfile(config.input_path):
        raise FileNotFoundError(f"Trajectory file not found: {config.input_path}")

    input_description = os.path.abspath(config.input_path)

    universe = mda.Universe(config.input_path, format="LAMMPSDUMP")

    # print("Trajectory loaded.")
    # print("Number of atoms:", len(universe.atoms))
    # print("Number of frames:", len(universe.trajectory))

    all_types = np.asarray(universe.atoms.types, dtype=int)

    # print("Available atom types:", np.unique(all_types))

    if not np.any(all_types == int(config.oxygen_type)):

        raise ValueError(
            f"No atoms with type {config.oxygen_type} found. "
            f"Available types: {np.unique(all_types)}"
        )

    frame_count = 0

    for ts in universe.trajectory:
        timestep = ts.data.get("step", ts.frame)
        frame_data = {
            "ids": universe.atoms.ids,
            "positions": universe.atoms.positions,
            "types": universe.atoms.types,
            "box": ts.dimensions,
        }
        process_frame(frame_data, ts.frame, timestep, config, stats, analysis_function)

        frame_count += 1

    return frame_count, input_description
