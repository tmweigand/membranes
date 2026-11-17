import membranes


def test_init_coordination_number():

    input = "test_data/"
    output = "tests_out"
    simtype = "test_sim"

    config = membranes.input_config.InputConfig(
        input,
        output,
        simtype,
    )

    analysis = membranes.analysis.CoordinationNumber(config)
    for input_file in config.get_input_files():
        oxygen_data, _, _ = membranes.io.read_lammps_file(
            file=input_file,
            config=config,
            extract_atom=analysis.oxygen_type,
        )

        n_pairs = analysis.process_frame(oxygen_data)

        print(n_pairs)
        # NEED to ASSERT something
