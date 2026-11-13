from helpers.trajectory import read_trajectory_frames
from helpers.frame_files import read_frame_files


def get_simulation_frames(input_mode, config, stats, analysis_function, process_frame):

    if input_mode == "trajectory":

        return read_trajectory_frames(config, stats, analysis_function, process_frame)

    elif input_mode == "frame_files":

        return read_frame_files(config, stats, analysis_function, process_frame)

    else:
        raise ValueError(f"Unsupported input mode: {input_mode}")
