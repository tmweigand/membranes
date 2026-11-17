# Instructions to perform Coordination number analysis 

## Example command:

Water simulation using a trajectory file:

python coordination_number/coordination_number.py --simtype watersim --input_mode trajectory --input path/to/traj.lammpstrj --output path/to/output/ --oxygen-type 1

Membrane simulation using individual frame files:

python coordination_number/coordination_number.py --simtype membranesim --input_mode frame_files --input path/to/frame_files/ --output path/to/output/ --oxygen-type 10 --box "x" "y" "z" 90 90 90

Details of the arguments:

1. Enter the type of simulation using `--simtype "membranesim"` or `--simtype "watersim"`.
2. Select the input format using `--input_mode "trajectory"` or `--input_mode "frame_files"`.
   - For `trajectory`, provide the location of the `.lammpstrj` file using `--input`.
   - For `frame_files`, provide the directory containing the individual frame files using `--input`.
3. Provide the output directory using `--output`.
4. Specify the oxygen atom type using `--oxygen-type`, depending on the simulation type.
5. If using `frame_files` input mode, provide the simulation box dimensions using `--box` (in Angstrom, followed by the three box angles in degrees).
6. If region-based analysis needs to be performed, the region definitions and bounds are selected automatically based on `--simtype` in `input_type/config.py`.
   - For `watersim`, the entire simulation box (`whole_box`) is considered.
   - For `membranesim`, the feed, membrane, and permeate regions are considered. Their bounds can be modified in `input_type/config.py`.
7. Specify the coordination cutoff range using `--cutoff-min` and `--cutoff-max` in Angstrom. The cutoff increment is set to 0.2 Angstrom in `input_type/config.py`. The cutoff values are converted internally to nm for the analysis.
8. For `frame_files` input mode, the PMMoTo LAMMPS reader is used. Its parameters (`initial_iteration`, `simulation_timestep`, and `unwrapped_coordinates`) are specified in `input_type/config.py`.
9. The default cutoff ranges are 1–8 Å for `watersim` and 3–8 Å for `membranesim`. These can be overridden using `--cutoff-min` and `--cutoff-max`.
