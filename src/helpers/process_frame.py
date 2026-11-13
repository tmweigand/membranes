from helpers.atom_data import prepare_frame


def process_frame(frame_data, frame_index, timestep, config, stats, analysis_function):
    prepared_data = prepare_frame(frame_data)

    ### This function needs to be updated/modified based on the type of analysis being performed.
    analysis_function(prepared_data, config, stats, frame_index, timestep)
