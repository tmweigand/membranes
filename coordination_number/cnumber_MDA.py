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
cutoff = 3.6  # Angstrom = 0.36 nm

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

per_water_file = os.path.join(output_dir, "water_coordination_per_frame.csv")
summary_file = os.path.join(output_dir, "water_coordination_summary.csv")

region_names = ["feed", "membrane", "permeate"]

region_sum = {r: 0.0 for r in region_names}
region_sum_sq = {r: 0.0 for r in region_names}
region_count = {r: 0 for r in region_names}

# Coordination-number frequency for each region
region_cn_frequency = {region: {} for region in region_names}

# Mean coordination number for each region in each frame
frame_region_means = {region: [] for region in region_names}


# --------------------------------------------------
# GET INPUT FILES USING PMMOTO WORKFLOW
# --------------------------------------------------
water_files = helpers.get_files(machine="mac", file_type="gzip")

water_files = sorted(water_files, key=lambda x: int(x.split(".")[1]))

# --------------------------------------------------
# PROCESS EVERY FRAME
# --------------------------------------------------
with open(per_water_file, "w", newline="") as fout:

    writer = csv.writer(fout)

    writer.writerow(
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
        # Returns unique pairs within 3.6 Angstrom.
        # Applies periodic minimum-image distances.
        pairs = self_capped_distance(
            oxygen_pos, max_cutoff=cutoff, box=box, return_distances=False
        )

        # Each pair contributes one neighbor to each oxygen
        cn = np.zeros(len(oxygen_pos), dtype=np.int64)

        if len(pairs) > 0:
            np.add.at(cn, pairs[:, 0], 1)
            np.add.at(cn, pairs[:, 1], 1)

        # ----------------------------------------------
        # WRITE PER-OXYGEN RESULTS
        # ----------------------------------------------
        writer.writerows(
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

        # ----------------------------------------------
        # ACCUMULATE REGIONAL STATISTICS
        # ----------------------------------------------
        for region in region_names:

            values = cn[regions == region]

            region_sum[region] += np.sum(values)
            region_sum_sq[region] += np.sum(values**2)
            region_count[region] += len(values)

        print(
            "Frame:",
            frame,
            "Timestep:",
            timestep,
            "Oxygen atoms:",
            len(oxygen_pos),
            "Neighbor pairs:",
            len(pairs),
            "Feed:",
            np.sum(regions == "feed"),
            "Membrane:",
            np.sum(regions == "membrane"),
            "Permeate:",
            np.sum(regions == "permeate"),
        )

        for region in region_names:
            values = cn[regions == region]

            unique_cn, counts = np.unique(values, return_counts=True)

            for coordination, count in zip(unique_cn, counts):
                coordination = int(coordination)
                frequency = region_cn_frequency[region]
                frequency[coordination] = frequency.get(coordination, 0) + int(count)

            if len(values) > 0:
                frame_region_means[region].append(np.mean(values))
            else:
                frame_region_means[region].append(np.nan)

# --------------------------------------------------
# WRITE REGIONAL SUMMARY
# --------------------------------------------------
with open(summary_file, "w", newline="") as fout:

    writer = csv.writer(fout)

    writer.writerow(
        [
            "region",
            "water_frame_observations",
            "mean_coordination_number",
            "std_coordination_number",
        ]
    )

    for region in region_names:

        n = region_count[region]

        if n == 0:
            mean_cn = np.nan
            std_cn = np.nan
        else:
            mean_cn = region_sum[region] / n

            variance = region_sum_sq[region] / n - mean_cn**2

            std_cn = np.sqrt(max(variance, 0.0))

        writer.writerow([region, n, mean_cn, std_cn])

        print(region, "Mean CN:", mean_cn, "Std CN:", std_cn, "Observations:", n)

print("Saved per-water results:", per_water_file)
print("Saved regional summary:", summary_file)


# # --------------------------------------------------
# # Plot 1: Frequency vs coordination number
# # Plot 2: Probability vs coordination number
# # --------------------------------------------------

# for region in region_names:
#     frequency = region_cn_frequency[region]

#     if not frequency:
#         continue

#     coordination_numbers = np.array(sorted(frequency.keys()))
#     frequencies = np.array(
#         [frequency[cn] for cn in coordination_numbers],
#         dtype=np.int64
#     )

#     probabilities = frequencies / frequencies.sum()

#     # Frequency plot
#     plt.figure()
#     plt.bar(coordination_numbers, frequencies, width=0.8)
#     plt.xlabel("Coordination number")
#     plt.ylabel("Frequency")
#     plt.title("Coordination number frequency: " + region)
#     plt.tight_layout()
#     plt.savefig(
#         os.path.join(
#             output_dir,
#             "coordination_frequency_" + region + ".png"
#         ),
#         dpi=300
#     )
#     plt.close()

#     # Probability plot
#     plt.figure()
#     plt.bar(coordination_numbers, probabilities, width=0.8)
#     plt.xlabel("Coordination number")
#     plt.ylabel("Probability")
#     plt.title("Coordination number probability: " + region)
#     plt.tight_layout()
#     plt.savefig(
#         os.path.join(
#             output_dir,
#             "coordination_probability_" + region + ".png"
#         ),
#         dpi=300
#     )
#     plt.close()

# # --------------------------------------------------
# # Plot 3: Mean coordination number vs frame
# # --------------------------------------------------

# plt.figure()

# for region in region_names:
#     plt.plot(
#         range(len(frame_region_means[region])),
#         frame_region_means[region],
#         label=region
#     )

# plt.xlabel("Frame")
# plt.ylabel("Mean coordination number")
# plt.title("Mean coordination number vs frame")
# plt.legend()
# plt.tight_layout()
# plt.savefig(
#     os.path.join(output_dir, "mean_coordination_vs_frame.png"),
#     dpi=300
# )
# plt.close()

# print("All coordination plots saved in:", output_dir)
