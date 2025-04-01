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

    n_files = 4

    membrane_file = glob.glob("data_out/membrane_bins/bin_membrane.txt")
    water_file = glob.glob("data_out/water_bins/bin_water.txt")

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

    plt.plot(membrane_bins, membrane_density / bin_volume)
    plt.plot(water_bins, water_density / bin_volume)
    plt.plot(membrane_bins, membrane_density / bin_volume + water_density / bin_volume)
    plt.show()


if __name__ == "__main__":
    generate_plots()
