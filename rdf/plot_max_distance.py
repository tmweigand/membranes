import numpy as np
import glob
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def read_results(filename):
    results = []
    connected = []
    count = 0
    with open(filename, "r") as f:
        for line in f:
            if "All:" in line:
                parts = line.split("All:")
                if len(parts) > 1:
                    try:
                        result = float(
                            parts[1].strip().split()[0]
                        )  # Get the number after "All:"
                        results.append(result)
                    except ValueError:
                        pass  # Ignore lines where parsing fails
            if "Connected:" in line:
                parts = line.split("Connected:")
                if len(parts) > 1:
                    try:
                        result2 = float(
                            parts[1].strip().split()[0]
                        )  # Get the number after "All:"
                        connected.append(result2)
                    except ValueError:
                        pass  # Ignore lines where parsing fails
            if "All:" in line and "Connected:" in line:
                if np.abs(result - result2) > 1.0e-3:
                    count += 1

    print(count, len(results), count / len(results) * 100)

    return results, connected


if __name__ == "__main__":
    file_in = [
        # "rdf/bridges_results/max_distance_1.4.out",
        "rdf/bridges_results/max_distance_3.61.out",
        # "rdf/bridges_results/max_distance_17.315.out",
    ]

    fig, ax = plt.subplots()

    G = [1.40, 3.61, 17.32]  # Order according to sort
    labels = [f"G = {g:.2f} kJ/mol" for g in G]
    colors = ["tab:orange", "tab:blue", "tab:green"]

    # First pass: collect all times
    all_times = []
    for file in file_in:
        pmf_results, connected_results = read_results(file)
        print(
            np.average(pmf_results),
            np.median(pmf_results),
            np.max(pmf_results),
            np.min(pmf_results),
        )
        all_times.append(pmf_results)

    plot = False
    if plot:

        # Flatten all_times and define common bins
        all_data = np.concatenate(all_times)
        min_time = np.min(all_data)
        max_time = np.max(all_data)
        bins = np.linspace(min_time, max_time, 71)  # 50 bins
        bin_width = bins[1] - bins[0]
        bin_centers = bins[:-1] + bin_width / 2

        # Offset settings
        n = len(G)
        offset_width = bin_width / n

        # Create subplots
        fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True)

        colors = ["tab:orange", "tab:blue", "tab:green"]

        for ax, label, times, color in zip(axes, labels, all_times, colors):
            ax.hist(
                times,
                bins=bins,
                density=True,
                color=color,
                edgecolor="black",
                label=label,
            )
            ax.tick_params(axis="y", labelsize=18)

            # Custom legend: text only, no color/symbol
            text_only_handle = Patch(facecolor="none", edgecolor="none", label=label)
            ax.legend(handles=[text_only_handle], prop={"size": 18}, frameon=False)

        axes[-1].set_xlabel("Maximum Pore Size (Å)", fontsize=18)
        axes[-1].tick_params(axis="x", labelsize=18)

        # Add single y-axis label
        fig.text(
            0.04,
            0.5,
            "Probability Density",
            va="center",
            rotation="vertical",
            fontsize=16,
        )

        plt.tight_layout(rect=[0.06, 0, 1, 1])  # leave space for y-label

        # plt.title("Minimum G(r) for Connected Path (kj/mol)")
        plt.savefig(
            "data_out/max_distance_in_dense_region.pdf",
            dpi=300,
            bbox_inches="tight",
        )
