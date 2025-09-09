"""plot_densities.py"""

import glob
import numpy as np
import pmmoto
import matplotlib.pyplot as plt
from typing import Literal, Dict


def generate_plots():
    """
    Loop through generated rdf data for each atom and plot
    """
    box = [
        [0.0, 176],
        [0.0, 176],
        [-288, 238],  # Ignore water reservoirs
    ]

    n_files = 13947  # From files processed - see /ocean/projects/cts200024p/ajotcham/membranes/mass_output.out
    membrane_file = glob.glob(
        "density_calculations/bridges_results/membrane_bins/bin_membrane.txt"
    )
    water_file = glob.glob(
        "density_calculations/bridges_results/water_bins/bin_water.txt"
    )

    membrane_bins, membrane_density = np.genfromtxt(membrane_file[0], unpack=True)
    water_bins, water_density = np.genfromtxt(water_file[0], unpack=True)

    avogadro = 6.022e23  # Avogadro's number
    angstrom_to_cm = 1e-8
    bin_width = membrane_bins[1] - membrane_bins[0]
    bin_area = (box[1][1] - box[1][0]) * (box[0][1] - box[0][0])
    bin_volume = (bin_width * bin_area) * (angstrom_to_cm**3)

    # print(bin_width, bin_area, bin_volume)

    membrane_density = membrane_density / avogadro / n_files
    water_density = water_density / avogadro / n_files

    fig, ax = plt.subplots()
    fig.set_figwidth(9)
    fig.set_figheight(5)

    ax.plot(
        membrane_bins,
        membrane_density / bin_volume + water_density / bin_volume,
        "k",
        label="Total",
    )
    ax.plot(
        membrane_bins, membrane_density / bin_volume, color="0.5", label="Polyamide"
    )
    ax.plot(water_bins, water_density / bin_volume, "blue", label="Water")

    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    ax.set_xlabel("z-coordinate (Å)", fontsize=18)
    ax.set_ylabel("Density (g/cm$^3$)", fontsize=18)

    ax.set_xlim(left=-200, right=150)
    ax.set_ylim(bottom=0.0)

    # Grid, legend
    ax.grid(True, linestyle=":", linewidth=1)
    ax.legend(loc="lower left", fontsize=18)

    # Add vertical lines
    ax.axvline(x=-100, color="orange", linewidth=1.5)
    ax.axvline(x=100, color="orange", linewidth=1.5)

    # ax.axvline(x=-35, color="salmon", linestyle="dashed", linewidth=2)
    # ax.axvline(x=65, color="salmon", linestyle="dashed", linewidth=2)

    # Save the figure
    fig.tight_layout()
    plt.savefig(
        "density_calculations/membrane_density.png",
        dpi=300,
        bbox_inches="tight",
    )


if __name__ == "__main__":
    generate_plots()
