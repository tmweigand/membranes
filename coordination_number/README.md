Instructions to perform Coordination number analysis 

1) Enter type of simulation with a flag --simtype "membranesim" or "watersim". Then give the input information in src/membranes/input_type/ accordingly.
2) For "membranesim": if input is old bridges data where each frame is a separate file, then set "use_helpers_get_files" to "True", else set "use_helpers_get_files" to "False" and give path to the files location in "frame_files_pattern".
3) For "watersim": give the location of the .lammpstrj file in "trajectory_file".
4) In both the cases, give the location where the output should be written in "output_dir".
5) Give molecule ID of oxygens in "oxygen_type" depending on the --simtype.
6) If region based analysis needs to be performed: give the region definition and bounds in "region_bounds".
7) "cutoff_nm" range should be defined to do a sweep of all the cutoffs for coordination number analysis.
8) For input_type == membranesim: it uses pmmoto lammps reader to read the input files, so set initial_iteration = 0, simulation_timestep = 1e-6, unwrapped_coordinates = False. Also, give the box dimensions in "frame_files_box".

