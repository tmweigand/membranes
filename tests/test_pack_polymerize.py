import membranes


# def test_add_molecule():
#     style_props = membranes.domain_generation.system.StyleProperties()
#     polymerizer = membranes.domain_generation.pack_polymerize.Polymerizer(
#         style_props.lammps
#     )

#     # Initialize simulation box
#     # polymerizer.lammps.cmd.region("box", "block", 0, 1, 0, 1, 0, 1)
#     # polymerizer.lammps.cmd.create_box(50, "box")

#     polymerizer.read_tmc("domain_in/TMC_converted.lmps")
#     polymerizer.add_molecule("mpd", "domain_in/MPD.mol")


# def test_read_tmc():
#     style_props = membranes.domain_generation.system.StyleProperties()
#     polymerizer = membranes.domain_generation.pack_polymerize.Polymerizer(
#         style_props.lammps
#     )

#     polymerizer.read_tmc("domain_in/TMC_converted.lmps")

#     # print(box)
