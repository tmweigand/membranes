import os
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np

cutoffs_nm = np.arange(30, 81, 2, dtype=np.int64) / 100.0
output_dir = "data_out/coordination_number"
chunksize = 100_000

regions = ["feed", "membrane", "permeate"]

for cutoff_nm in cutoffs_nm:
    cutoff_label = f"{cutoff_nm:.2f}nm"
    filename = os.path.join(
        output_dir, f"water_coordination_per_frame_cutoff_{cutoff_label}.csv"
    )
    cn_frequency = {region: Counter() for region in regions}

    # Read only the columns required for plotting
    for chunk in pd.read_csv(
        filename, usecols=["region", "coordination_number"], chunksize=chunksize
    ):
        for region in regions:
            values = chunk.loc[chunk["region"] == region, "coordination_number"]

            cn_frequency[region].update(values.dropna().astype(int).tolist())

    # Plot frequency and probability
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    for region in regions:
        counts = cn_frequency[region]
        if not counts:
            continue

        cn_values = sorted(counts)
        frequencies = [counts[cn] for cn in cn_values]
        total = sum(frequencies)
        probabilities = [count / total for count in frequencies]

        frequency_file = os.path.join(
            output_dir,
            f"coordination_number_frequency_{region}_cutoff_{cutoff_label}.dat",
        )
        probability_file = os.path.join(
            output_dir,
            f"coordination_number_probability_{region}_cutoff_{cutoff_label}.dat",
        )
        np.savetxt(frequency_file, np.column_stack((cn_values, frequencies)))
        print("Saved coordination number frequency data to", frequency_file)
        np.savetxt(probability_file, np.column_stack((cn_values, probabilities)))
        print("Saved coordination number probability data to", probability_file)

        axes[0].plot(cn_values, frequencies, marker="o", label=region)
        axes[1].plot(cn_values, probabilities, marker="o", label=region)

    axes[0].set_xlabel("Coordination number")
    axes[0].set_ylabel("Frequency")
    axes[0].set_title(f"Frequency vs. coordination number ({cutoff_label})")

    axes[1].set_xlabel("Coordination number")
    axes[1].set_ylabel("Probability")
    axes[1].set_title(f"Probability vs. coordination number ({cutoff_label})")

    for ax in axes:
        ax.legend()
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plot_file = os.path.join(
        output_dir, f"coordination_number_plots_cutoff_{cutoff_label}.png"
    )
    plt.savefig(plot_file, dpi=300)
    plt.close(fig)
    print("Saved coordination number plots to", plot_file)
