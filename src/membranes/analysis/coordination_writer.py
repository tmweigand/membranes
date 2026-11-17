from .coordination_stats import CoordinationStats
from ..input_config import InputConfig
from ..io import safe_label

import os


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
                region_label = safe_label(region)
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
