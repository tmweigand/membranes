import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np

filename = "data_out/coordination_number/water_coordination_per_frame.csv"
chunksize = 100_000

regions = ["feed", "membrane", "permeate"]
cn_frequency = {region: Counter() for region in regions}

# Read only the columns required for plotting
for chunk in pd.read_csv(
    filename,
    usecols=["region", "coordination_number"],
    chunksize=chunksize
):
    for region in regions:
        values = chunk.loc[
            chunk["region"] == region,
            "coordination_number"
        ]

        cn_frequency[region].update(
            values.dropna().astype(int).tolist()
        )

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

    np.savetxt("data_out/coordination_number/coordination_number_frequency_"+region+".dat", np.column_stack((cn_values, frequencies)))
    print ("Saved coordination number frequency data to data_out/coordination_number/coordination_number_frequency_"+region+".dat" )
    np.savetxt("data_out/coordination_number/coordination_number_probability_"+region+".dat", np.column_stack((cn_values, probabilities)))
    print ("Saved coordination number probability data to data_out/coordination_number/coordination_number_probability_"+region+".dat" )

    axes[0].plot(
        cn_values, frequencies, marker="o", label=region
    )
    axes[1].plot(
        cn_values, probabilities, marker="o", label=region
    )

axes[0].set_xlabel("Coordination number")
axes[0].set_ylabel("Frequency")
axes[0].set_title("Frequency vs. coordination number")

axes[1].set_xlabel("Coordination number")
axes[1].set_ylabel("Probability")
axes[1].set_title("Probability vs. coordination number")

for ax in axes:
    ax.legend()
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("data_out/coordination_number/coordination_number_plots.png", dpi=300)
plt.show()

