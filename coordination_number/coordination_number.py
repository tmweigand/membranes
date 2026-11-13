"""
water_coordination.py

Calculate water-oxygen coordination numbers for either:
  1. A multi-frame LAMMPS dump read by MDAnalysis, or
  2. A folder of single-frame PMMoTo-compatible files.

Outputs (under output_dir):
  - coordination_summary.dat
      One row per cutoff and region: observations, mean CN, standard deviation.
  - coordination_distribution_<cutoff>_<region>.dat
      Coordination number, frequency, and probability for each cutoff/region.
  - README.txt
      Records the input configuration and describes the outputs.

Distances in the input coordinates and box lengths are assumed to be Angstrom.
"""

import os
from collections import Counter

import argparse
import importlib

import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import self_capped_distance

from helpers.process_frame import process_frame
from helpers.atom_data import select_atoms_by_type
from helpers.regions import classify_regions
from helpers.frame_files import safe_label
from helpers.input_handler import get_simulation_frames
from input_type import config


def calculate_oxygen_pairs(oxygen_positions, box, max_cutoff_nm):
    pairs, distances = self_capped_distance(
        oxygen_positions,
        max_cutoff=float(max_cutoff_nm) * 10.0,
        box=box,
        return_distances=True,
    )
    return pairs, distances


def coordination_numbers(pair_indices, pair_distances, cutoff_nm, oxygen_count):
    cutoff_angstrom = np.nextafter(float(cutoff_nm) * 10.0, np.inf)
    selected_pairs = pair_indices[pair_distances <= cutoff_angstrom]
    coordination = np.zeros(oxygen_count, dtype=np.int64)

    if len(selected_pairs):
        np.add.at(coordination, selected_pairs[:, 0], 1)
        np.add.at(coordination, selected_pairs[:, 1], 1)

    return coordination


# This function is called for each frame of the simulation and goes to get_simulation_frames function, and updates the statistics for each cutoff and region
def analyze_coordination_frame(frame_data, config, stats, frame_index, timestep):
    # Select oxygen atoms using the configured atom TYPE.
    oxygen_data = select_atoms_by_type(frame_data, config.oxygen_type)
    oxygen_ids = oxygen_data["ids"]
    oxygen_positions = oxygen_data["positions"]

    if oxygen_positions.size == 0:
        raise ValueError(
            f"No oxygen atoms found in frame {frame_index}. "
            f"Check oxygen_type={config.oxygen_type}."
        )

    if len(np.unique(oxygen_ids)) != len(oxygen_ids):
        raise ValueError(f"Duplicate oxygen atom IDs in frame {frame_index}.")

    # Classify oxygen atoms into regions.
    regions = classify_regions(oxygen_positions[:, 2], config.region_bounds)

    # Calculate oxygen-oxygen pairs only once per frame.
    pairs, distances = calculate_oxygen_pairs(
        oxygen_positions, frame_data["box"], max(config.cutoffs_nm)
    )

    print(
        f"Frame {frame_index} | timestep {timestep} | "
        f"oxygen atoms {len(oxygen_positions)} | "
        f"pairs at {max(config.cutoffs_nm):.2f} nm: {len(pairs)}"
    )

    # Calculate coordination numbers at each cutoff.
    for cutoff_nm in config.cutoffs_nm:
        cn = coordination_numbers(pairs, distances, cutoff_nm, len(oxygen_positions))
        # Update statistics separately for each region.
        for region in config.region_bounds:
            update_stats(stats[float(cutoff_nm)][region], cn[regions == region])


def update_stats(stats, values):
    values = np.asarray(values, dtype=np.int64)
    if values.size == 0:
        return
    stats["sum"] += int(np.sum(values, dtype=np.int64))
    stats["sum_sq"] += float(np.sum(values.astype(np.float64) ** 2))
    stats["count"] += int(values.size)
    stats["frequency"].update(values.tolist())


def initialize_stats(config):
    stats = {
        float(cutoff): {
            region: {"sum": 0, "sum_sq": 0.0, "count": 0, "frequency": Counter()}
            for region in config.region_bounds
        }
        for cutoff in config.cutoffs_nm
    }
    return stats


def write_summary(stats, config):
    summary_path = os.path.join(config.output_dir, "coordination_summary.dat")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("# cutoff_nm mean_coordination_number std_coordination_number region\n")
        for cutoff_nm in config.cutoffs_nm:
            for region in config.region_bounds:
                item = stats[float(cutoff_nm)][region]
                count = item["count"]
                if count:
                    mean = item["sum"] / count
                    variance = max(item["sum_sq"] / count - mean**2, 0.0)
                    std = np.sqrt(variance)
                else:
                    mean = np.nan
                    std = np.nan
                f.write(f"{cutoff_nm:.2f} " f"{mean:.8f} " f"{std:.8f} " f"{region}\n")

                print(
                    f"Cutoff {cutoff_nm:.2f} nm | {region} | "
                    f"N={count} | mean CN={mean:.4f} | "
                    f"std={std:.4f}"
                )
    return summary_path


def write_distributions(stats, config):
    for cutoff_nm in config.cutoffs_nm:
        cutoff_label = f"{cutoff_nm:.2f}nm"
        for region in config.region_bounds:
            frequency = stats[float(cutoff_nm)][region]["frequency"]
            region_label = safe_label(region)

            distribution_path = os.path.join(
                config.output_dir,
                f"coordination_distribution_" f"{cutoff_label}_{region_label}.dat",
            )
            total = sum(frequency.values())
            with open(distribution_path, "w", encoding="utf-8") as f:
                for cn_value in sorted(frequency):
                    count = frequency[cn_value]
                    probability = count / total if total else np.nan
                    f.write(f"{cn_value:d} " f"{count:d} " f"{probability:.12g}\n")


def write_readme(input_description, frame_count, config):
    path = os.path.join(config.output_dir, "OUTPUT_README.txt")

    with open(path, "w", encoding="utf-8") as f:
        f.write("Water oxygen coordination analysis\n")
        f.write("==================================\n\n")
        f.write(f"Input data: {input_description}\n")
        f.write(f"Frames processed: {frame_count}\n")
        f.write("Coordinate and box-length units: Angstrom\n")
        f.write("Cutoff units: nm (converted to Angstrom for distance calculations)\n")
        f.write(
            f"Cutoffs (nm): " f"{', '.join(f'{x:.2f}' for x in config.cutoffs_nm)}\n\n"
        )
        f.write("Region definitions (z_min, z_max), Angstrom:\n")
        for region, limits in config.region_bounds.items():
            f.write(f"  {region}: {limits}\n")
        f.write("\nOutputs\n" "-------\n")
        f.write(
            "coordination_summary.dat: mean and standard deviation for each cutoff and region.\n"
        )
        f.write(
            "coordination_distribution_<cutoff>_<region>.dat:coordination number, frequency and probability.\n"
        )
        f.write(
            f"\nAll outputs are stored in: " f"{os.path.abspath(config.output_dir)}\n"
        )


def write_outputs(stats, input_description, frame_count, config):
    summary_path = write_summary(stats, config)
    write_distributions(stats, config)
    write_readme(input_description, frame_count, config)
    print("\nCoordination analysis completed.")
    print("Frames processed:", frame_count)
    print("Summary:", os.path.abspath(summary_path))
    print(
        "README:", os.path.abspath(os.path.join(config.output_dir, "OUTPUT_README.txt"))
    )
    print("All output files:", os.path.abspath(config.output_dir))


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Calculate water oxygen coordination numbers."
    )
    parser.add_argument(
        "--simtype",
        type=str,
        required=True,
        choices=["watersim", "membranesim"],
        help="Choose simulation type: watersim or membranesim",
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to the input simulation folder or trajectory file",
    )
    parser.add_argument("--output", required=True, help="Path to the output directory")
    parser.add_argument(
        "--input_mode",
        choices=["trajectory", "frame_files"],
        required=True,
        help="Input mode: trajectory file or individual frame files",
    )
    parser.add_argument(
        "--oxygen-type", type=int, required=True, help="Atom type of water oxygen"
    )

    parser.add_argument(
        "--box",
        type=float,
        nargs=6,
        default=None,
        metavar=("LX", "LY", "LZ", "ALPHA", "BETA", "GAMMA"),
        help="Box dimensions in Angstrom and angles in degrees",
    )

    parser.add_argument(
        "--cutoff-min",
        type=float,
        default=None,
        help="Minimum coordination cutoff in Angstrom",
    )

    parser.add_argument(
        "--cutoff-max",
        type=float,
        default=None,
        help="Maximum coordination cutoff in Angstrom",
    )
    return parser.parse_args()


# MAIN


def main():

    args = parse_arguments()

    # Set input and output paths from command-line arguments
    config.input_path = args.input
    config.output_dir = args.output

    # Simulation type and regions
    config.simtype = args.simtype
    config.region_bounds = config.get_region_bounds(args.simtype)

    # Oxygen atom type
    config.oxygen_type = args.oxygen_type

    # Box dimensions
    config.frame_files_box = (
        np.asarray(args.box, dtype=np.float64) if args.box is not None else None
    )

    # Default cutoff ranges preserve your previous settings
    if args.simtype == "watersim":
        cutoff_min_A = 1.0 if args.cutoff_min is None else args.cutoff_min
    else:
        cutoff_min_A = 3.0 if args.cutoff_min is None else args.cutoff_min

    cutoff_max_A = 8.0 if args.cutoff_max is None else args.cutoff_max

    if cutoff_min_A <= 0 or cutoff_max_A < cutoff_min_A:
        raise ValueError("Invalid coordination cutoff range.")

    # Convert Angstrom to nm
    config.cutoffs_nm = (
        np.arange(
            cutoff_min_A, cutoff_max_A + config.cutoff_step_A / 2, config.cutoff_step_A
        )
        / 10.0
    )

    # Validate box requirement for frame-file simulations
    if args.input_mode == "frame_files" and config.frame_files_box is None:
        raise ValueError("--box must be provided when using frame_files input mode.")

    os.makedirs(config.output_dir, exist_ok=True)
    if args.input_mode not in {"trajectory", "frame_files"}:
        raise ValueError("input_mode must be 'trajectory' or 'frame_files'.")

    if not config.region_bounds:
        raise ValueError("region_bounds must contain at least one region.")

    # Each cutoff has aggregate stats and a CN histogram for each region.
    stats = initialize_stats(config)

    # Read frames and process them, updating stats.
    frame_count, input_description = get_simulation_frames(
        args.input_mode, config, stats, analyze_coordination_frame, process_frame
    )

    # Summary statistics and output files
    write_outputs(stats, input_description, frame_count, config)


if __name__ == "__main__":
    main()
