"""coordination_number.py"""

import argparse
import numpy as np
import membranes
import time

from membranes.input_config import SYSTEMS


def validate_frame(atom_data, input_file, atom_type):
    """Raise ValueError if a frame is unusable."""
    if atom_data["positions"].size == 0:
        raise ValueError(
            f"No oxygen atoms found in frame {input_file}. "
            f"Check oxygen_type={atom_type}."
        )
    ids = atom_data["id"]
    if len(np.unique(ids)) != len(ids):
        raise ValueError(f"Duplicate oxygen atom IDs in frame {input_file}.")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Calculate water oxygen coordination numbers."
    )
    parser.add_argument(
        "--system",
        required=True,
        choices=list(SYSTEMS),
        help="Which system to analyze",
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to the input simulation folder or trajectory file",
    )
    parser.add_argument("--output", required=True, help="Path to the output directory")

    return parser.parse_args()


def main():

    args = parse_arguments()

    # Set input and output paths from command-line arguments
    config = membranes.input_config.InputConfig(
        args.input,
        args.output,
        args.system,
    )

    region_classifier = membranes.RegionClassifier(
        bounds=config.sim.region_bounds, box=membranes.Box(domain=config.sim.box)
    )

    stats = membranes.analysis.CoordinationStats(
        cutoffs=config.sim.cutoffs, regions=region_classifier
    )

    writer = membranes.analysis.CoordinationWriter(
        in_dir=args.input,
        out_dir=args.output,
        stats=stats,
    )

    for input_file in config.get_input_files():

        oxygen_data, _, timestep = membranes.io.read_lammps_file(
            file=input_file,
            config=config,
            extract_atom=config.sim.oxygen_type,
        )

        validate_frame(oxygen_data, input_file, config.sim.oxygen_type)

        stats.update(oxygen_data["positions"], oxygen_data["id"])

    writer.plots()


if __name__ == "__main__":
    main()
