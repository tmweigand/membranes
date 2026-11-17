"""coordination_number.py"""

import os
import argparse
from collections import Counter
import numpy as np
from MDAnalysis.lib.distances import self_capped_distance
import membranes
from membranes.input_config import InputConfig

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


class CoordinationStats:
    """Accumulates coordination-number observations per (cutoff, region).

    Cutoffs are identified by their integer index into `cutoffs_angstrom`,
    which avoids floating-point dictionary keys.
    """

    def __init__(self, cutoffs_angstrom, regions):
        self.cutoffs_angstrom = np.asarray(cutoffs_angstrom, dtype=np.float64)
        self.regions = list(regions)
        self._data = {
            (i, region): {"sum": 0, "sum_sq": 0.0, "count": 0, "frequency": Counter()}
            for i in range(len(self.cutoffs_angstrom))
            for region in self.regions
        }

    @property
    def cutoffs_nm(self):
        return self.cutoffs_angstrom / ANGSTROM_PER_NM

    def update(self, cutoff_index, region, values):
        values = np.asarray(values, dtype=np.int64)
        if values.size == 0:
            return
        item = self._data[(cutoff_index, region)]
        item["sum"] += int(values.sum())
        item["sum_sq"] += float(np.sum(values.astype(np.float64) ** 2))
        item["count"] += int(values.size)
        item["frequency"].update(values.tolist())

    def summary(self, cutoff_index, region):
        """Return count, mean, and (population) std for one cutoff/region."""
        item = self._data[(cutoff_index, region)]
        count = item["count"]
        if count == 0:
            return {"count": 0, "mean": np.nan, "std": np.nan}
        mean = item["sum"] / count
        variance = max(item["sum_sq"] / count - mean**2, 0.0)
        return {"count": count, "mean": mean, "std": float(np.sqrt(variance))}

    def distribution(self, cutoff_index, region):
        """Return {cn: {"count": int, "probability": float}} sorted by cn."""
        frequency = self._data[(cutoff_index, region)]["frequency"]
        total = sum(frequency.values())
        return {
            cn: {
                "count": count,
                "probability": count / total if total else np.nan,
            }
            for cn, count in sorted(frequency.items())
        }


class CoordinationNumber:
    """Computes coordination numbers frame by frame and feeds the stats."""

    def __init__(self, config: InputConfig):
        self.config = config
        sim = config.sim
        self.oxygen_type = sim.oxygen_type
        self.region_bounds = sim.region_bounds
        self.frames_processed = 0
        self.box = membranes.Box(domain=config.box)
        self.regions = membranes.RegionClassifier(self.region_bounds, self.box)

        start, stop = sim.cutoffs
        step = sim.cutoff_step
        if start <= 0 or stop < start:
            raise ValueError("Invalid coordination cutoff range.")
        if step <= 0:
            raise ValueError("cutoff_step must be positive.")

        cutoffs = np.arange(start, stop + step / 2, step)
        self.stats = CoordinationStats(cutoffs, self.region_bounds)

    @property
    def max_cutoff_angstrom(self):
        return float(self.stats.cutoffs_angstrom.max())

    def find_pairs(self, positions):
        """Neighbour search only: all pairs within the largest cutoff."""
        # Small pad so pairs exactly at the max cutoff are not lost to round-off.
        return self_capped_distance(
            positions,
            max_cutoff=self.max_cutoff_angstrom + 1e-6,
            box=self.box.MD_box_from_domain(),
            return_distances=True,
        )

    def accumulate_frame(self, pairs, distances, regions, atom_count):
        """Compute CN at every cutoff and add it to the stats, per region."""
        for i, cutoff in enumerate(self.stats.cutoffs_angstrom):
            cn = coordination_numbers(pairs, distances, cutoff, atom_count)
            for region in self.region_bounds:
                self.stats.update(i, region, cn[regions == region])

    def process_frame(self, atom_data):
        """Full pipeline for one frame. Returns the number of pairs found."""
        positions = atom_data["coords_unwrapped"]

        labels = self.regions.classify_regions(positions[:, 2])

        pairs, distances = self.find_pairs(positions)
        self.accumulate_frame(pairs, distances, labels, len(positions))
        self.frames_processed += 1
        return len(pairs)


class CoordinationWriter:
    """Writes summary, distributions, and README for a finished analysis."""

    def __init__(
        self,
        config: InputConfig,
        stats: CoordinationStats,
        region_bounds,
        frames_processed,
    ):
        self.config = config
        self.stats = stats
        self.region_bounds = region_bounds
        self.frames_processed = frames_processed

        os.makedirs(self.config.out_dir, exist_ok=True)

    def write_all(self):
        self.write_summary()
        self.write_distributions()
        self.write_readme()

    def write_summary(self):
        path = os.path.join(self.config.out_dir, "coordination_summary.dat")
        with open(path, "w", encoding="utf-8") as f:
            f.write(
                "# cutoff_nm mean_coordination_number "
                "std_coordination_number observations region\n"
            )
            for i, cutoff_nm in enumerate(self.stats.cutoffs_nm):
                for region in self.region_bounds:
                    s = self.stats.summary(i, region)
                    f.write(
                        f"{cutoff_nm:.2f} {s['mean']:.8f} {s['std']:.8f} "
                        f"{s['count']:d} {region}\n"
                    )

    def write_distributions(self):
        for i, cutoff_nm in enumerate(self.stats.cutoffs_nm):
            cutoff_label = f"{cutoff_nm:.2f}nm"
            for region in self.region_bounds:
                region_label = membranes.io.safe_label(region)
                path = os.path.join(
                    self.config.out_dir,
                    f"coordination_distribution_{cutoff_label}_{region_label}.dat",
                )
                distribution = self.stats.distribution(i, region)
                with open(path, "w", encoding="utf-8") as f:
                    f.write("# coordination_number frequency probability\n")
                    for cn, item in distribution.items():
                        f.write(
                            f"{cn:d} {item['count']:d} {item['probability']:.12g}\n"
                        )

    def write_readme(self):
        path = os.path.join(self.config.out_dir, "README.txt")
        cutoffs = ", ".join(f"{x:.2f}" for x in self.stats.cutoffs_nm)
        with open(path, "w", encoding="utf-8") as f:
            f.write("Water oxygen coordination analysis\n")
            f.write("==================================\n\n")
            f.write(f"Input data: {self.config.in_dir}\n")
            f.write(f"Frames processed: {self.frames_processed}\n")
            f.write("Coordinate and box-length units: Angstrom\n")
            f.write("Cutoff units: nm (converted to Angstrom for distances)\n")
            f.write(f"Cutoffs (nm): {cutoffs}\n\n")
            f.write("Region definitions (z_min, z_max), Angstrom:\n")
            for region, limits in self.region_bounds.items():
                f.write(f"  {region}: {limits}\n")
            f.write("\nOutputs\n-------\n")
            f.write(
                "coordination_summary.dat: observations, mean and standard "
                "deviation for each cutoff and region.\n"
            )
            f.write(
                "coordination_distribution_<cutoff>_<region>.dat: "
                "coordination number, frequency and probability.\n"
            )
            f.write(
                "\nNotes\n-----\n"
                "Standard deviations are population standard deviations over "
                "all oxygen observations pooled across frames. Consecutive "
                "frames are correlated, so this is not an estimate from "
                "independent samples.\n"
            )
            f.write(
                f"\nAll outputs are stored in: {os.path.abspath(self.config.out_dir)}\n"
            )

    def write_outputs(self):
        self.write_summary()
        self.write_distributions()
        self.write_readme()


def coordination_numbers(pair_indices, pair_distances, cutoff_nm, atom_count):
    cutoff_angstrom = np.nextafter(float(cutoff_nm) * 10.0, np.inf)
    selected_pairs = pair_indices[pair_distances <= cutoff_angstrom]
    coordination = np.zeros(atom_count, dtype=np.int64)

    if len(selected_pairs):
        np.add.at(coordination, selected_pairs[:, 0], 1)
        np.add.at(coordination, selected_pairs[:, 1], 1)

    return coordination


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

    analysis = CoordinationNumber(config)

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
    CoordinationWriter(
        config,
        analysis.stats,
        analysis.region_bounds,
        analysis.frames_processed,
    ).write_all()


if __name__ == "__main__":
    main()
