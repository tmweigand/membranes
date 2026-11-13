import os
import glob
import re
import pmmoto


def get_timestep(filename):
    numbers = re.findall(r"\d+", os.path.basename(filename))

    if not numbers:
        raise ValueError(f"Could not identify timestep in {filename}")

    return int(numbers[-1])


def get_frame_files(config):
    search_pattern = os.path.join(config.input_path, config.file_pattern)
    frame_files = glob.glob(search_pattern)
    if not frame_files:
        raise FileNotFoundError(f"No frame files found matching: {search_pattern}")
    frame_files.sort(key=get_timestep)

    return frame_files


def read_frame_files(config, stats, analysis_function, process_frame):
    frame_files = get_frame_files(config)
    frame_count = 0
    for frame_index, path in enumerate(frame_files):
        atom_ids, positions, types, masses, domain, timestep = (
            pmmoto.io.data_read.read_lammps_atoms(
                path,
                unwrapped_coordinates=config.unwrapped_coordinates,
                initial_iteration=config.initial_iteration,
                simulation_timestep=config.simulation_timestep,
            )
        )
        frame_data = {
            "ids": atom_ids,
            "positions": positions,
            "types": types,
            "box": config.frame_files_box,
        }

        process_frame(
            frame_data, frame_index, timestep, config, stats, analysis_function
        )

        frame_count += 1
    input_description = config.input_path

    return frame_count, input_description


def safe_label(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("_")
