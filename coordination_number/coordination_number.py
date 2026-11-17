"""coordination_number.py"""

import argparse
import numpy as np
import membranes

ANGSTROM_PER_NM = 10.0


def validate_frame(atom_data, input_file, atom_type):
    """Raise ValueError if a frame is unusable."""
    if atom_data["coords_unwrapped"].size == 0:
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

    analysis = membranes.analysis.CoordinationNumber(config)

    for input_file in config.get_input_files():
        oxygen_data, _, timestep = membranes.io.read_lammps_file(
            file=input_file,
            config=config,
            extract_atom=analysis.oxygen_type,
        )
        validate_frame(oxygen_data, input_file, analysis.oxygen_type)

        n_pairs = analysis.process_frame(oxygen_data)

        print(
            f"Frame {analysis.frames_processed} | timestep {timestep} | "
            f"oxygen atoms {len(oxygen_data['coords_unwrapped'])} | "
            f"pairs at {analysis.max_cutoff_angstrom / ANGSTROM_PER_NM:.2f} nm: "
            f"{n_pairs}"
        )

    # Summary statistics and output files
    membranes.analysis.CoordinationWriter(
        config,
        analysis.stats,
        analysis.region_bounds,
        analysis.frames_processed,
    ).write_all()


if __name__ == "__main__":
    main()
