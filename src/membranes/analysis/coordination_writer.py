from .coordination_stats import CoordinationStats
import os


class CoordinationWriter:
    """Writes summary, distributions, and README for a finished analysis."""

    def __init__(
        self,
        in_dir,
        out_dir,
        frames_processed,
        stats: CoordinationStats,
    ):
        self.stats = stats
        self.frames_processed = frames_processed
        self.in_dir = in_dir
        self.out_dir = out_dir

        os.makedirs(self.out_dir, exist_ok=True)

    def _write_summary(self):
        path = os.path.join(self.out_dir, "coordination_summary.dat")
        summary = self.stats.summary()

        with open(path, "w", encoding="utf-8") as f:
            f.write(
                "# cutoff region mean_coordination_number "
                "std_coordination_number observations\n"
            )
            for (cutoff, region), s in summary.items():
                f.write(
                    f"{cutoff:.2f} {region} "
                    f"{s['mean']:.8f} {s['std']:.8f} {s['count']:d}\n"
                )

    def _write_distributions(self):
        distributions = self.stats.distribution()
        for (cutoff, region), distribution in distributions.items():
            label = f"{cutoff:.2f}_{region}"

            path = os.path.join(
                self.out_dir,
                f"coordination_distribution_{label}.dat",
            )

            with open(path, "w", encoding="utf-8") as f:
                f.write("# coordination_number probability frequency\n")

                for cn, item in distribution.items():
                    f.write(
                        f"{cn:d} {item['probability']:.12g} " f"{item['count']:d}\n"
                    )

    def _write_readme(self):
        path = os.path.join(self.out_dir, "README.txt")
        cutoffs = (
            f"{self.stats.cutoffs:.2f}"
            if isinstance(self.stats.cutoffs, (float, int))
            else ", ".join(f"{x:.2f}" for x in self.stats.cutoffs)
        )
        with open(path, "w", encoding="utf-8") as f:
            f.write("Water oxygen coordination analysis\n")
            f.write("==================================\n\n")
            f.write(f"Input data: {self.in_dir}\n")
            f.write(f"Frames processed: {self.frames_processed}\n")
            f.write("Coordinate and box-length units: Angstrom\n")
            f.write(f"Cutoffs (A): {cutoffs}\n\n")
            if self.stats.regions is not None:
                f.write("Region definitions (z_min, z_max), Angstrom:\n")
                for region, limits in self.stats.regions.bounds.items():
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
            f.write(f"\nAll outputs are stored in: {os.path.abspath(self.out_dir)}\n")

    def write_outputs(self):
        self._write_summary()
        self._write_distributions()
        self._write_readme()
