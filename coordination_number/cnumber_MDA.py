"""water_coordination_mdanalysis.py"""

import os
import csv
import numpy as np
import pmmoto
from membranes.helpers import helpers
import matplotlib.pyplot as plt

from MDAnalysis.lib.distances import self_capped_distance

# --------------------------------------------------
# SETTINGS
# --------------------------------------------------
oxygen_type = 15
cutoffs_nm = np.arange(30, 81, 2, dtype=np.int64) / 100.0

z_lower = -100.0
z_upper = 100.0

# Periodic box dimensions in Angstrom
box_lower = np.array([0.0, 0.0, -287.0])
box_upper = np.array([175.76827931, 175.76827931, 237.0])
box_lengths = box_upper - box_lower

# MDAnalysis box: [lx, ly, lz, alpha, beta, gamma]
box = np.array([box_lengths[0], box_lengths[1], box_lengths[2], 90.0, 90.0, 90.0])

output_dir = "data_out/coordination_number"
os.makedirs(output_dir, exist_ok=True)

region_names = ["feed", "membrane", "permeate"]


def coordination_numbers(pair_indices, pair_distances, cutoff_nm, oxygen_count):
    """Count neighbors for each oxygen at a cutoff given in nanometers."""
    cutoff_angstrom = np.nextafter(cutoff_nm * 10.0, np.inf)
    selected_pairs = pair_indices[pair_distances <= cutoff_angstrom]
    coordination = np.zeros(oxygen_count, dtype=np.int64)

    if len(selected_pairs) > 0:
        np.add.at(coordination, selected_pairs[:, 0], 1)
        np.add.at(coordination, selected_pairs[:, 1], 1)

    return coordination


# --------------------------------------------------
# GET INPUT FILES USING PMMOTO WORKFLOW
# --------------------------------------------------
water_files = helpers.get_files(machine="mac", file_type="gzip")

water_files = sorted(water_files, key=lambda x: int(x.split(".")[1]))

# --------------------------------------------------
# PROCESS EVERY FRAME
# --------------------------------------------------
cutoff_stats = {
    cutoff_nm: {
        region: {"sum": 0.0, "sum_sq": 0.0, "count": 0} for region in region_names
    }
    for cutoff_nm in cutoffs_nm
}

file_handles = []
writers = {}
summary_writers = {}

try:
    for cutoff_nm in cutoffs_nm:
        cutoff_label = f"{cutoff_nm:.2f}nm"
        per_water_file = os.path.join(
            output_dir, f"water_coordination_per_frame_cutoff_{cutoff_label}.csv"
        )
        summary_file = os.path.join(
            output_dir, f"water_coordination_summary_cutoff_{cutoff_label}.csv"
        )
        per_water_handle = open(per_water_file, "w", newline="")
        file_handles.append(per_water_handle)
        per_water_writer = csv.writer(per_water_handle)
        per_water_writer.writerow(
            [
                "frame",
                "timestep",
                "oxygen_atom_id",
                "x_A",
                "y_A",
                "z_A",
                "region",
                "coordination_number",
            ]
        )
        writers[cutoff_nm] = per_water_writer
        summary_handle = open(summary_file, "w", newline="")
        file_handles.append(summary_handle)
        summary_writers[cutoff_nm] = csv.writer(summary_handle)

    for frame, file in enumerate(water_files):

        # Read the PMMoTo-compatible input file
        atom_ids, positions, types, masses, domain, timestep = (
            pmmoto.io.data_read.read_lammps_atoms(
                file,
                unwrapped_coordinates=False,
                initial_iteration=71045000,
                simulation_timestep=2e-6,
            )
        )

        atom_ids = np.asarray(atom_ids)
        positions = np.asarray(positions)
        types = np.asarray(types)

        # ----------------------------------------------
        # SELECT WATER OXYGEN ATOMS
        # ----------------------------------------------
        oxygen_mask = types == oxygen_type

        oxygen_atom_ids = atom_ids[oxygen_mask]
        oxygen_pos = np.asarray(positions[oxygen_mask], dtype=np.float64)

        if len(oxygen_pos) == 0:
            raise ValueError("No oxygen atoms found. Check oxygen_type = 15.")

        if len(np.unique(oxygen_atom_ids)) != len(oxygen_atom_ids):
            raise ValueError("Duplicate oxygen atom IDs found.")

        # ----------------------------------------------
        # CLASSIFY REGIONS USING CURRENT POSITIONS
        # ----------------------------------------------
        z = oxygen_pos[:, 2]

        regions = np.full(len(z), "membrane", dtype=object)

        regions[z < z_lower] = "feed"
        regions[z > z_upper] = "permeate"

        # ----------------------------------------------
        # FIND NEIGHBORING OXYGEN PAIRS
        # ----------------------------------------------
        # Calculate pair distances once at the largest cutoff (0.80 nm = 8 A).
        pairs, distances = self_capped_distance(
            oxygen_pos,
            max_cutoff=cutoffs_nm[-1] * 10.0,
            box=box,
            return_distances=True,
        )

        print(
            "Frame:",
            frame,
            "Timestep:",
            timestep,
            "Oxygen atoms:",
            len(oxygen_pos),
            "Pairs within 0.80 nm:",
            len(pairs),
            "Feed:",
            np.sum(regions == "feed"),
            "Membrane:",
            np.sum(regions == "membrane"),
            "Permeate:",
            np.sum(regions == "permeate"),
        )

        for cutoff_nm in cutoffs_nm:
            cn = coordination_numbers(pairs, distances, cutoff_nm, len(oxygen_pos))
            writers[cutoff_nm].writerows(
                (
                    frame,
                    timestep,
                    oxygen_atom_ids[i],
                    oxygen_pos[i, 0],
                    oxygen_pos[i, 1],
                    oxygen_pos[i, 2],
                    regions[i],
                    cn[i],
                )
                for i in range(len(cn))
            )

            for region in region_names:
                values = cn[regions == region]
                stats = cutoff_stats[cutoff_nm][region]
                stats["sum"] += np.sum(values)
                stats["sum_sq"] += np.sum(values**2)
                stats["count"] += len(values)

    for cutoff_nm in cutoffs_nm:
        summary_writer = summary_writers[cutoff_nm]
        summary_writer.writerow(
            [
                "region",
                "water_frame_observations",
                "mean_coordination_number",
                "std_coordination_number",
            ]
        )

        for region in region_names:
            stats = cutoff_stats[cutoff_nm][region]
            count = stats["count"]
            if count == 0:
                mean_cn = np.nan
                std_cn = np.nan
            else:
                mean_cn = stats["sum"] / count
                variance = stats["sum_sq"] / count - mean_cn**2
                std_cn = np.sqrt(max(variance, 0.0))

            summary_writer.writerow([region, count, mean_cn, std_cn])
            print(
                f"Cutoff {cutoff_nm:.2f} nm, {region}: Mean CN {mean_cn}, "
                f"Std CN {std_cn}, Observations {count}"
            )

        print(
            "Saved per-water results:",
            os.path.join(
                output_dir,
                f"water_coordination_per_frame_cutoff_{cutoff_nm:.2f}nm.csv",
            ),
        )
        print(
            "Saved regional summary:",
            os.path.join(
                output_dir,
                f"water_coordination_summary_cutoff_{cutoff_nm:.2f}nm.csv",
            ),
        )
finally:
    for file_handle in file_handles:
        file_handle.close()
