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
import glob
import re
from collections import Counter

import argparse
import importlib

import numpy as np
import pmmoto
import MDAnalysis as mda
from MDAnalysis.lib.distances import self_capped_distance

from helpers.box import validate_box
from helpers.frame_files import get_frame_files, safe_label
from helpers.regions import classify_regions


def coordination_numbers(pair_indices, pair_distances, cutoff_nm, oxygen_count):
    cutoff_angstrom = np.nextafter(float(cutoff_nm) * 10.0, np.inf)
    selected_pairs = pair_indices[pair_distances <= cutoff_angstrom]
    coordination = np.zeros(oxygen_count, dtype=np.int64)

    if len(selected_pairs):
        np.add.at(coordination, selected_pairs[:, 0], 1)
        np.add.at(coordination, selected_pairs[:, 1], 1)

    return coordination


def process_frame(
    frame_index,
    timestep,
    atom_ids,
    positions,
    types,
    box,
    oxygen_type,
    region_bounds,
    cutoffs_nm,
    stats,
):
    """Process one frame and update coordination statistics."""

    atom_ids = np.asarray(atom_ids)
    positions = np.asarray(positions, dtype=np.float64)
    types = np.asarray(types)

    oxygen_mask = types.astype(int) == int(oxygen_type)

    oxygen_ids = atom_ids[oxygen_mask]
    oxygen_pos = positions[oxygen_mask]

    if oxygen_pos.size == 0:
        raise ValueError(
            f"No oxygen atoms found in frame {frame_index}. "
            f"Check oxygen_type={oxygen_type}."
        )

    if len(np.unique(oxygen_ids)) != len(oxygen_ids):
        raise ValueError(f"Duplicate oxygen atom IDs in frame {frame_index}.")

    box = validate_box(box, f"frame {frame_index}")

    regions = classify_regions(oxygen_pos[:, 2], region_bounds)

    # Find pairs once at the largest cutoff.
    pairs, distances = self_capped_distance(
        oxygen_pos,
        max_cutoff=float(cutoffs_nm[-1]) * 10.0,
        box=box,
        return_distances=True,
    )

    print(
        f"Frame {frame_index} | timestep {timestep} | "
        f"oxygen atoms {len(oxygen_pos)} | "
        f"pairs at {cutoffs_nm[-1]:.2f} nm {len(pairs)}"
    )

    for cutoff_nm in cutoffs_nm:

        cn = coordination_numbers(pairs, distances, cutoff_nm, len(oxygen_pos))

        for region in region_bounds:
            update_stats(stats[float(cutoff_nm)][region], cn[regions == region])


def update_stats(stats, values):
    values = np.asarray(values, dtype=np.int64)
    if values.size == 0:
        return
    stats["sum"] += int(np.sum(values, dtype=np.int64))
    stats["sum_sq"] += float(np.sum(values.astype(np.float64) ** 2))
    stats["count"] += int(values.size)
    stats["frequency"].update(values.tolist())


def write_readme(
    input_description,
    frame_count,
    oxygen_type_value,
    region_config,
    input_mode,
    cutoffs_nm,
    output_dir,
):
    path = os.path.join(output_dir, "README.txt")
    with open(path, "w", encoding="utf-8") as f:
        f.write("Water oxygen coordination analysis\n")
        f.write("==================================\n\n")
        f.write(f"Input mode: {input_mode}\n")
        f.write(f"Input data: {input_description}\n")
        f.write(f"Frames processed: {frame_count}\n")
        f.write(f"Oxygen atom TYPE: {oxygen_type_value}\n")
        f.write("Coordinate and box-length units: Angstrom\n")
        f.write("Cutoff units: nm (converted to Angstrom for distance calculations)\n")
        f.write(f"Cutoffs (nm): {', '.join(f'{x:.2f}' for x in cutoffs_nm)}\n\n")
        f.write("Region definitions (z_min, z_max), Angstrom:\n")
        for region, limits in region_config.items():
            f.write(f"  {region}: {limits}\n")

        f.write(f"\nAll outputs are stored in: {os.path.abspath(output_dir)}\n")


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

    return parser.parse_args()


# MAIN


def main():

    args = parse_arguments()

    config = importlib.import_module(f"input_type.{args.simtype}")

    print(f"\nSimulation type: {args.simtype}")
    print(f"Configuration: {config.__name__}")

    os.makedirs(config.output_dir, exist_ok=True)

    if config.input_mode not in {"trajectory", "frame_files"}:
        raise ValueError("input_mode must be 'trajectory' or 'frame_files'.")

    if not config.region_bounds:
        raise ValueError("region_bounds must contain at least one region.")

    # Each cutoff has aggregate stats and a CN histogram for each region.
    stats = {
        float(cutoff): {
            region: {"sum": 0, "sum_sq": 0.0, "count": 0, "frequency": Counter()}
            for region in config.region_bounds
        }
        for cutoff in config.cutoffs_nm
    }

    frame_count = 0
    input_description = ""

    if config.input_mode == "trajectory":
        if not os.path.isfile(config.trajectory_file):
            raise FileNotFoundError(
                f"Trajectory file not found: {config.trajectory_file}"
            )

        input_description = os.path.abspath(config.trajectory_file)
        universe = mda.Universe(config.trajectory_file, format="LAMMPSDUMP")
        print("Trajectory loaded.")
        print("Number of atoms:", len(universe.atoms))
        print("Number of frames:", len(universe.trajectory))
        print(
            "Available atom types:",
            np.unique(np.asarray(universe.atoms.types, dtype=int)),
        )

        all_types = np.asarray(universe.atoms.types, dtype=int)
        if not np.any(all_types == int(config.oxygen_type)):
            raise ValueError(
                f"No atoms with type {config.oxygen_type} found. "
                f"Available types: {np.unique(all_types)}"
            )

        for ts in universe.trajectory:
            timestep = ts.data.get("step", ts.frame)
            process_frame(
                frame_index=ts.frame,
                timestep=timestep,
                atom_ids=universe.atoms.ids,
                positions=universe.atoms.positions,
                types=universe.atoms.types,
                box=ts.dimensions,
                oxygen_type=config.oxygen_type,
                region_bounds=config.region_bounds,
                cutoffs_nm=config.cutoffs_nm,
                stats=stats,
            )

            frame_count += 1

    else:
        frame_files = get_frame_files(
            pattern=config.frame_files_pattern,
            use_helpers_get_files=config.use_helpers_get_files,
            machine=config.machine,
            file_type=config.file_type,
        )
        if config.use_helpers_get_files:
            input_source = f"helpers.get_files(machine={config.machine!r}, file_type={config.file_type!r})"
        else:
            input_source = f"Glob pattern: {config.frame_files_pattern}"
        input_description = (
            f"{input_source}\n"
            f"First file: {os.path.abspath(frame_files[0])}\n"
            f"Last file: {os.path.abspath(frame_files[-1])}"
        )
        print("Number of single-frame files:", len(frame_files))

        for frame_index, path in enumerate(frame_files):
            atom_ids, positions, types, masses, domain, timestep = (
                pmmoto.io.data_read.read_lammps_atoms(
                    path,
                    unwrapped_coordinates=config.unwrapped_coordinates,
                    initial_iteration=config.initial_iteration,
                    simulation_timestep=config.simulation_timestep,
                )
            )
            process_frame(
                frame_index=frame_index,
                timestep=timestep,
                atom_ids=atom_ids,
                positions=positions,
                types=types,
                box=config.frame_files_box,
                oxygen_type=config.oxygen_type,
                region_bounds=config.region_bounds,
                cutoffs_nm=config.cutoffs_nm,
                stats=stats,
            )

            frame_count += 1

    # Summary statistics and output files
    summary_path = os.path.join(config.output_dir, "coordination_summary.dat")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write(
            "# cutoff_nm mean_coordination_number " "std_coordination_number region\n"
        )
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
                f.write(
                    f"{cutoff_nm:.2f} "
                    f"{mean:.8f} {std:.8f} {region}\n"
                    # f"{cutoff_nm:.2f} {region} {count:d} " f"{mean:.8f} {std:.8f}\n"
                )
                print(
                    f"Cutoff {cutoff_nm:.2f} nm | {region} | "
                    f"N={count} | mean CN={mean:.4f} | std={std:.4f}"
                )

    # Single file outputs for each cutoff and region: coordination number, frequency, probability
    for cutoff_nm in config.cutoffs_nm:
        cutoff_label = f"{cutoff_nm:.2f}nm"
        for region in config.region_bounds:
            frequency = stats[float(cutoff_nm)][region]["frequency"]
            region_label = safe_label(region)
            distribution_path = os.path.join(
                config.output_dir,
                f"coordination_distribution_{cutoff_label}_{region_label}.dat",
            )
            total = sum(frequency.values())
            with open(distribution_path, "w", encoding="utf-8") as f:
                # f.write("# coordination_number frequency probability\n")
                for cn_value in sorted(frequency):
                    count = frequency[cn_value]
                    probability = count / total if total else np.nan
                    f.write(f"{cn_value:d} {count:d} {probability:.12g}\n")

    write_readme(
        input_description,
        frame_count,
        config.oxygen_type,
        config.region_bounds,
        config.input_mode,
        config.cutoffs_nm,
        config.output_dir,
    )

    print("\nCoordination analysis completed.")
    print("Frames processed:", frame_count)
    print("Summary:", os.path.abspath(summary_path))
    print("README:", os.path.abspath(os.path.join(config.output_dir, "README.txt")))
    print("All output files:", os.path.abspath(config.output_dir))


if __name__ == "__main__":
    main()
