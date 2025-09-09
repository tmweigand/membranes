import numpy as np
import glob
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def read_results(filename):
    results = []
    with open(filename, "r") as f:
        for line in f:
            if "Result:" in line:
                parts = line.split("Result:")
                if len(parts) > 1:
                    try:
                        result = int(parts[1].strip())
                        results.append(result)
                    except ValueError:
                        pass  # Ignore lines where result isn't an integer
    return results


if __name__ == "__main__":
    file_in = [
        "rdf/bridges_results/persistent_paths_random_1.4.out",
        "rdf/bridges_results/persistent_paths_random_3.61.out",
        "rdf/bridges_results/persistent_paths_random_17.315.out",
    ]

    fig, ax = plt.subplots()

    G = [1.40, 3.61, 17.32]  # Order according to sort
    labels = [f"G = {g:.2f} kJ/mol" for g in G]
    colors = ["tab:orange", "tab:blue", "tab:green"]

    # First pass: collect all times
    all_times = []
    for pmf_results in file_in:
        pmf_results = read_results(pmf_results)
        times = np.asarray(pmf_results) * 0.01
        print(np.average(times), np.median(times), np.max(times), np.min(times))
        all_times.append(times)

    # Flatten all_times and define common bins
    all_data = np.concatenate(all_times)
    min_time = 0
    max_time = np.max(all_data)
    bins = np.linspace(min_time, max_time, 71)  # 50 bins
    bin_width = bins[1] - bins[0]
    bin_centers = bins[:-1] + bin_width / 2

    # Offset settings
    n = len(G)
    offset_width = bin_width / n

    # # Plot each histogram with offset
    # for i, (label, times, color) in enumerate(zip(labels, all_times, colors)):
    #     counts, _ = np.histogram(times, bins=bins, density=False)
    #     offset = (i - n / 2) * offset_width + offset_width / 2
    #     ax.bar(
    #         bin_centers + offset,
    #         np.log(counts),
    #         width=offset_width,
    #         label=label,
    #         color=color,
    #         edgecolor="black",
    #     )

    # for i, (label, times, color) in enumerate(zip(labels, all_times, colors)):

    #     ax.hist(
    #         times,
    #         bins=bins,
    #         density=True,
    #         alpha=0.6,
    #         facecolor=color,
    #         edgecolor="black",
    #         label=label,
    #     )

    # ax.set_xlim(left=0)
    # ax.legend(prop={"size": 12})
    # ax.set_xlabel("Persistence Time (ns)", fontsize=16)
    # ax.set_ylabel("Probability Density", fontsize=16)
    # ax.tick_params(axis="x", labelsize=16)
    # ax.tick_params(axis="y", labelsize=16)

    # Create subplots
    fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True)

    colors = ["tab:orange", "tab:blue", "tab:green"]

    for ax, label, times, color in zip(axes, labels, all_times, colors):
        ax.hist(
            times, bins=bins, density=True, color=color, edgecolor="black", label=label
        )
        ax.tick_params(axis="y", labelsize=18)

        # Custom legend: text only, no color/symbol
        text_only_handle = Patch(facecolor="none", edgecolor="none", label=label)
        ax.legend(handles=[text_only_handle], prop={"size": 18}, frameon=False)

    axes[-1].set_xlabel("Persistence Time (ns)", fontsize=18)
    axes[-1].tick_params(axis="x", labelsize=18)

    # Add single y-axis label
    fig.text(
        0.04, 0.5, "Probability Density", va="center", rotation="vertical", fontsize=16
    )

    plt.tight_layout(rect=[0.06, 0, 1, 1])  # leave space for y-label

    # plt.title("Minimum G(r) for Connected Path (kj/mol)")
    plt.savefig(
        "data_out/persistent_path_length.pdf",
        dpi=300,
        bbox_inches="tight",
    )
