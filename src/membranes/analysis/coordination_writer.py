import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

from .coordination_stats import CoordinationStats

RC_PARAMS = {
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "legend.fontsize": 9,
    "legend.frameon": False,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.6,
    "figure.constrained_layout.use": True,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
}


class CoordinationWriter:
    """Writes summary, distributions, README, and plots for a finished analysis."""

    def __init__(
        self,
        in_dir,
        out_dir,
        stats: CoordinationStats,
        per_file_distributions=True,
    ):
        self.stats = stats
        self.frames_processed = stats.counter
        self.in_dir = in_dir
        self.out_dir = out_dir
        self.per_file_distributions = per_file_distributions
        self.plots_dir = os.path.join(out_dir, "plots")

        os.makedirs(self.plots_dir, exist_ok=True)

    # ------------------------------------------------------------------ helpers

    def _region_names(self):
        if self.stats.regions is None:
            return ["all"]
        return list(self.stats.regions.bounds)

    def _cutoff_list(self):
        return np.atleast_1d(self.stats.cutoffs).tolist()

    def _colors(self):
        """One fixed colour per region so every figure uses the same mapping."""
        cmap = plt.get_cmap("tab10")
        return {r: cmap(i % 10) for i, r in enumerate(self._region_names())}

    @staticmethod
    def _label(cutoff, region):
        return f"{cutoff:.2f}_{region}".replace("/", "_").replace(" ", "_")

    @staticmethod
    def _save(fig, path):
        fig.savefig(path)
        plt.close(fig)

    def _plot_summary(self):
        """Mean coordination number vs cutoff, one line per region.

        The shaded band is +-1 SD of the pooled distribution (the spread of
        individual atoms), not an uncertainty on the mean.
        """
        summary = self.stats.summary()
        colors = self._colors()
        cutoffs = self._cutoff_list()

        with plt.rc_context(RC_PARAMS):
            fig, ax = plt.subplots(figsize=(6, 4.2))

            for region in self._region_names():
                pts = [
                    (c, summary[(c, region)])
                    for c in cutoffs
                    if summary[(c, region)]["count"]
                ]
                if not pts:
                    continue
                x = np.array([c for c, _ in pts])
                mean = np.array([s["mean"] for _, s in pts])
                std = np.array([s["std"] for _, s in pts])

                ax.plot(
                    x,
                    mean,
                    marker="o",
                    ms=5,
                    lw=1.8,
                    color=colors[region],
                    label=region,
                )
                ax.fill_between(
                    x, mean - std, mean + std, color=colors[region], alpha=0.15, lw=0
                )

            ax.set_xlabel("Cutoff (Å)")
            ax.set_ylabel("Mean coordination number")
            ax.set_title("Coordination number vs cutoff")
            ax.set_ylim(bottom=0)
            ax.legend(title="Region")
            ax.text(
                0.99,
                0.01,
                "shaded: ±1 SD of distribution",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=8,
                color="0.4",
            )
            self._save(fig, os.path.join(self.plots_dir, "coordination_summary.png"))

    def _plot_region(self):
        """One figure per cutoff with all regions overlaid on the same axes."""
        plot_dir = self.plots_dir + "/distributions"
        os.makedirs(plot_dir, exist_ok=True)

        distributions = self.stats.distribution()
        summary = self.stats.summary()
        colors = self._colors()

        with plt.rc_context(RC_PARAMS):
            for cutoff in self._cutoff_list():
                fig, ax = plt.subplots(figsize=(6.2, 4.2))
                drawn = False

                for region in self._region_names():
                    dist = distributions[(cutoff, region)]
                    if not dist:
                        continue
                    s = summary[(cutoff, region)]
                    cn = np.array(list(dist))
                    p = np.array([item["probability"] for item in dist.values()])

                    ax.fill_between(cn, p, color=colors[region], alpha=0.15, lw=0)
                    ax.plot(
                        cn,
                        p,
                        marker="o",
                        ms=4,
                        lw=1.8,
                        color=colors[region],
                        label=(
                            f"{region}  (mean={s['mean']:.2f}, "
                            f"SD={s['std']:.2f}, N={s['count']:,})"
                        ),
                    )
                    ax.axvline(s["mean"], color=colors[region], ls="--", lw=1.1)
                    drawn = True

                if not drawn:
                    plt.close(fig)
                    continue

                ax.set_title(f"Coordination number by region, cutoff {cutoff:.2f} Å")
                ax.set_xlabel("Coordination number")
                ax.set_ylabel("Probability")
                ax.set_ylim(bottom=0)
                ax.xaxis.set_major_locator(MaxNLocator(integer=True))
                ax.legend(title="Region (dashed: mean)")

                self._save(
                    fig,
                    os.path.join(
                        plot_dir,
                        f"coordination_regions_{cutoff:.2f}.png",
                    ),
                )

    def _plot_distributions(self):
        """Individual (cutoff, region) histograms with mean/SD annotation."""
        plot_dir = self.plots_dir + "/individual_distributions"
        os.makedirs(plot_dir, exist_ok=True)

        distributions = self.stats.distribution()
        summary = self.stats.summary()
        colors = self._colors()

        with plt.rc_context(RC_PARAMS):
            for (cutoff, region), dist in distributions.items():
                if region == "unassigned" or not dist:
                    continue

                s = summary[(cutoff, region)]
                cn = np.array(list(dist))
                p = np.array([item["probability"] for item in dist.values()])

                fig, ax = plt.subplots(figsize=(5, 3.6))
                ax.bar(cn, p, width=0.8, color=colors[region], alpha=0.85)
                ax.axvline(s["mean"], color="k", ls="--", lw=1.2)
                ax.axvspan(
                    s["mean"] - s["std"], s["mean"] + s["std"], color="k", alpha=0.07
                )

                ax.set_title(f"Cutoff {cutoff:.2f} Å, region {region}")
                ax.set_xlabel("Coordination number")
                ax.set_ylabel("Probability")
                ax.xaxis.set_major_locator(MaxNLocator(integer=True))
                ax.text(
                    0.97,
                    0.95,
                    f"mean = {s['mean']:.2f}\nSD = {s['std']:.2f}\nN = {s['count']:,}",
                    transform=ax.transAxes,
                    ha="right",
                    va="top",
                    fontsize=9,
                    bbox=dict(boxstyle="round", fc="white", ec="0.8"),
                )
                self._save(
                    fig,
                    os.path.join(
                        plot_dir,
                        f"coordination_distribution_{self._label(cutoff, region)}.png",
                    ),
                )

    def plots(self):
        self._plot_summary()
        self._plot_region()
        if self.per_file_distributions:
            self._plot_distributions()
