"""coordination_number.py"""

import argparse
import numpy as np
import membranes
import time

ANGSTROM_PER_NM = 10.0


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
        "--simtype",
        type=str,
        required=True,
        choices=["water_box", "vickers_membrane", "xu_membrane"],
        help="Choose simulation type: water_box or vickers_membrane or xu_membrane",
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
        args.simtype,
    )
    cutoff = np.arange(
        config.sim.cutoffs[0],
        config.sim.cutoffs[1] + config.sim.cutoff_step / 2,
        config.sim.cutoff_step,
    )

    region_classifier = membranes.RegionClassifier(
        bounds=config.sim.region_bounds, box=membranes.Box(domain=config.sim.box)
    )
    stats = membranes.analysis.CoordinationStats(
        cutoffs=cutoff, regions=region_classifier
    )
    for input_file in config.get_input_files():
        start = time.perf_counter()

        oxygen_data, _, timestep = membranes.io.read_lammps_file(
            file=input_file,
            config=config,
            extract_atom=config.sim.oxygen_type,
        )
        read_time = time.perf_counter() - start

        validate_frame(oxygen_data, input_file, config.sim.oxygen_type)

        start = time.perf_counter()

        stats.update(oxygen_data["positions"], oxygen_data["id"])
        analysis_time = time.perf_counter() - start

        print(
            f"Frame {stats.counter} | "
            f"Read: {read_time:.3f} s | "
            f"Analysis + stats: {analysis_time:.3f} s"
        )

        print(
            f"Frame {stats.counter} | timestep {timestep} | "
            f"oxygen atoms {len(oxygen_data['positions'])} | "
        )
    membranes.analysis.CoordinationWriter(
        in_dir=args.input,
        out_dir=args.output,
        frames_processed=stats.counter,
        stats=stats,
    ).write_outputs()


if __name__ == "__main__":
    main()
