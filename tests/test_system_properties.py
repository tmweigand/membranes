import membranes


def test_style_properties_init():
    lammps_wrapper = membranes.domain_generation.lammps_init.LAMMPSInitialize(
        log_file="test_log", units="real", atom_style="full"
    )

    system = membranes.domain_generation.system_properties.SystemProperties(
        lammps_wrapper.lmp, dimension=3, n_atom_types=2, boundary=("p", "p", "p")
    )

    box_length = [[-1, 1], [-1, 1], [-1, 1]]

    system.set_box(box_length=box_length, box_name="simbox")

    lmp_box = system.lmp.extract_box()

    # Lower
    assert lmp_box[0] == [-1.0, -1.0, -1.0]

    # Upper
    assert lmp_box[1] == [1.0, 1.0, 1.0]

    # Periodic?
    assert lmp_box[5] == [1, 1, 1]
