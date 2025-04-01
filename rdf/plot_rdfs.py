"""plot_rdfs.py"""

import glob
import numpy as np
import pmmoto
import matplotlib.pyplot as plt


density = {
    1: 0.00043010966306420854,
    2: 0.003319592339478703,
    3: 0.004498271614748887,
    4: 0.003755463286713287,
    5: 0.003015336935791481,
    6: 0.0037503973299427844,
    7: 0.003749702002542912,
    8: 0.0037620192307692307,
    9: 0.002356663223140496,
    10: 0.005113636363636364,
    11: 0.003735994119516847,
    12: 0.0004311029879211697,
    13: 0.0033264462809917354,
    14: 0.00416828909726637,
    15: 0.0004295136681500318,
    16: 0.00042802368086458995,
}


def generate_bin_plot(bin, folder):

    pmmoto.io.io_utils.check_file_path(folder)

    plt.figure(figsize=(8, 5))
    plt.plot(bin.centers, bin.values, linewidth=2, color="navy")

    plt.xlabel("Distance (Å)", fontsize=12)
    plt.ylabel("Observations", fontsize=12)
    plt.title(
        f"Distance Observation for Atom Type {bin.name}",
        fontsize=14,
    )
    plt.grid(True, alpha=0.3)
    plt.xlim(left=0)
    plt.ylim(bottom=0)
    plt.tight_layout()

    out_file = folder + f"bin_count_{bin.name}.pdf"
    plt.savefig(
        out_file,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


def generate_rdf_plot(bin, rdf, folder):

    pmmoto.io.io_utils.check_file_path(folder)

    plt.figure(figsize=(8, 5))
    plt.plot(bin.centers, rdf, linewidth=2, color="navy")

    plt.xlabel("Distance (Å)", fontsize=12)
    plt.ylabel("g(r)", fontsize=12)
    plt.title(f"RDF for Atom Type {bin.name}", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.xlim(left=0)
    plt.ylim(bottom=0)

    plt.tight_layout()
    out_file = folder + f"bin_count_{bin.name}.pdf"
    plt.savefig(
        out_file,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


def generate_plots():
    """
    Loop through generated rdf data for each atom and plot
    """

    bin_files = glob.glob("data_out/bins_extended/*")
    # bin_files = glob.glob("rdf/bridges_results/bins/*")

    for bin_file in bin_files:

        # Read the first two lines separately
        with open(bin_file, "r") as f:
            atom_type = f.readline().split(":")[-1]
            atom_label = int(f.readline().split(":")[-1])

        bins, binned_distances = np.genfromtxt(bin_file, skip_header=2, unpack=True)

        bin = pmmoto.analysis.bins.Bin(
            bins[0], bins[-1], len(bins), atom_type, binned_distances
        )

        # Convert to rdf
        rdf = bin.generate_rdf()

        new_rdf = rdf / density[atom_label]

        generate_bin_plot(bin, "data_out/bin_count_plots/")
        generate_rdf_plot(bin, rdf, "data_out/rdf_plots/")
        generate_rdf_plot(bin, new_rdf, "data_out/new_rdf_plots/")


if __name__ == "__main__":
    generate_plots()
