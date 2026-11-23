# Instructions to perform Coordination number analysis 

## Example command:

python coordination_number/coordination_number.py --simtype vickers_membrane --input data/water_data/ --output data_out/vickers_membrane/ 

Details of the arguments:

1. Enter the type of simulation using `--simtype "vickers_membrane"` or  `"water_box"` or `"xu_membrane"`.
2. Provide the output directory using `--output`.
3. Specify the oxygen atom type using `--oxygen-type` in `sim_types`, depending on the simulation type.
4. If region-based analysis needs to be performed, the region definitions and bounds are selected automatically based on `--simtype`.
   - For `water_box`, the entire simulation box (`whole_box`) is considered.
   - For `vickers_membrane`, the feed, membrane, and permeate regions are considered. Their bounds can be modified in `sim_types`.