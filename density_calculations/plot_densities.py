"""plot_densities.py"""

import glob
import numpy as np
import pmmoto
import matplotlib.pyplot as plt


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

    print(bin_width, bin_area, bin_volume)

    membrane_density = membrane_density / avogadro / n_files
    water_density = water_density / avogadro / n_files

    fig, ax = plt.subplots()
    fig.set_size_inches(10, 4)

    ax.plot(
        membrane_bins,
        membrane_density / bin_volume + water_density / bin_volume,
        "k",
        label="Total",
    )
    ax.plot(membrane_bins, membrane_density / bin_volume, color="0.5", label="Membrane")
    ax.plot(water_bins, water_density / bin_volume, "blue", label="Water")

    ax.set_xlabel("z-coordinate (Å)")
    ax.set_ylabel("Density (g/cm$^3$)")

    ax.set_xlim(left=-200, right=150)
    # ax.set_ylim(top=1.3)

    # Grid, legend
    ax.grid(True, linestyle=":", linewidth=1)
    ax.legend()

    # Add vertical lines
    ax.axvline(x=-100, color="gold", linestyle="dashed", linewidth=2)
    ax.axvline(x=100, color="gold", linestyle="dashed", linewidth=2)

    ax.axvline(x=-35, color="salmon", linestyle="dashed", linewidth=2)
    ax.axvline(x=65, color="salmon", linestyle="dashed", linewidth=2)

    # Save the figure
    plt.savefig("density_calculations/membrane_density.png")


if __name__ == "__main__":
    generate_plots()
