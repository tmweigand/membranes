"""plot_rdfs.py"""

import glob
import numpy as np
import pmmoto
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


def generate_bin_plot(bin, folder):

    pmmoto.io.io_utils.check_file_path(folder)

    plt.figure(figsize=(8, 5))
    plt.plot(bin.centers, bin.values, linewidth=2, color="navy")

    plt.xlabel("Distance (Å)", fontsize=12)
    plt.ylabel("ObservationCounts", fontsize=12)
    plt.title(
        f"Distance Observations between Water and Atom Type {bin.name}",
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


def generate_bin_and_rdf_plot(bin, rdf, folder):
    pmmoto.io.io_utils.check_file_path(folder)

    fig, ax1 = plt.subplots(figsize=(8, 5))

    # Plot bin counts
    color1 = "navy"
    ax1.plot(
        bin.centers,
        bin.values / 1e6,
        color=color1,
        linewidth=2,
        label="Observation Counts",
    )
    ax1.set_xlabel("Radial Distance (Å)", fontsize=16)
    ax1.set_ylabel("Observation Counts (10$^6$)", fontsize=16, color=color1)
    ax1.tick_params(axis="y", labelcolor=color1)
    ax1.set_xlim(left=0)
    ax1.set_ylim(bottom=0)
    ax1.grid(True, alpha=0.3)

    # Create second y-axis
    ax2 = ax1.twinx()
    color2 = "darkred"
    ax2.plot(rdf.radii, rdf.rdf, color=color2, linestyle="--", linewidth=2, label="RDF")
    ax2.set_ylabel("Radial Distribution Function", fontsize=16, color=color2)
    ax2.tick_params(axis="y", labelcolor=color2)
    ax2.set_ylim(bottom=0)

    # Title and layout
    plt.title(
        f"Water and {bin.name}",
        fontsize=16,
    )
    fig.tight_layout()

    ax1.tick_params(axis="x", labelsize=16)
    ax1.tick_params(axis="y", labelsize=16)
    ax2.tick_params(axis="y", labelsize=16)

    # Save
    out_file = folder + f"bin_count_and_rdf_{bin.name}.pdf"
    plt.savefig(
        out_file,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


def generate_rdf_plot(rdf, folder, x_line=None, y_line=None):
    """
    Generate and save an RDF plot.

    Parameters:
    - bin: Object containing bin information, including bin centers and name.
    - rdf: Array-like, radial distribution function values.
    - folder: String, directory where the plot will be saved.
    - x_line: Float (optional), location on the x-axis to draw a vertical line.
    """

    pmmoto.io.io_utils.check_file_path(folder)

    plt.figure(figsize=(8, 5))
    plt.plot(rdf.radii, rdf.rdf, linewidth=2, color="navy", label="RDF")

    # Draw vertical line if x_line is provided
    if x_line is not None:
        plt.axvline(
            x=x_line,
            color="red",
            linestyle="--",
            linewidth=1.5,
            label=f"r = {x_line:.2f}",
        )

    if y_line is not None:
        plt.axhline(
            y=y_line,
            color="blue",
            linestyle="--",
            linewidth=1.5,
            label=f"g = {y_line:.2f}",
        )

    plt.xlabel("Distance (Å)", fontsize=12)
    plt.ylabel("g(r)", fontsize=12)
    plt.title(f"RDF for Atom Type {rdf.name}", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.xlim(left=0)
    plt.ylim(bottom=0)

    plt.legend()
    plt.tight_layout()

    out_file = folder + f"bin_count_{rdf.name}.pdf"
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close()


def generate_free_energy_plot(rdf, folder, x_line=None, y_line=None):
    """
    Generate and save an RDF plot.

    Parameters:
    - bin: Object containing bin information, including bin centers and name.
    - rdf: Array-like, radial distribution function values.
    - folder: String, directory where the plot will be saved.
    - x_line: Float (optional), location on the x-axis to draw a vertical line.
    """

    pmmoto.io.io_utils.check_file_path(folder)

    plt.figure(figsize=(8, 5))
    plt.plot(
        rdf.radii,
        rdf.potential_mean_force(),
        linewidth=2,
        color="navy",
        label="Permeation Results",
    )

    # Draw vertical line if x_line is provided
    if x_line is not None:
        plt.axvline(
            x=x_line,
            color="red",
            linestyle="--",
            linewidth=1.5,
            label=f"Equilibrium Distance = {x_line:.2f} Å",
        )

    if y_line is not None:
        radius = rdf.interpolate_radius_from_pmf(y_line)
        plt.plot(radius, y_line, "ko", label=f"G({y_line:.2f}) = {radius:.2f} Å")
        # plt.axhline(
        #     y=y_line,
        #     color="blue",
        #     linestyle="--",
        #     linewidth=1.5,
        #     label=f"G({y_line:.2f}) = {radius:.2f} Å",
        # )

    plt.xlabel("Radial Distance (Å)", fontsize=16)
    plt.ylabel("G(r) (kJ/mol)", fontsize=16)
    plt.title(f"Potential of Mean Force for Water and {rdf.name}", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.xlim(left=0)
    # plt.ylim(top=20)

    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    plt.legend()
    plt.tight_layout()

    out_file = folder + f"G_{rdf.name}.pdf"
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close()


def generate_plots():
    """
    Loop through generated rdf data for each atom and plot
    """

    # bin_files = glob.glob("data_out/bins_extended/*")
    bin_files = glob.glob("rdf/bridges_results/bins_extended/*.rdf")

    atom_folder = "rdf/bridges_results/bins_extended/"
    atom_map, _ = pmmoto.io.data_read.read_rdf(atom_folder)

    for bin_file in bin_files:

        atom_type = bin_file.split("/")[-1].split(".")[0]
        print(atom_type)
        for key, values in atom_map.items():
            if values["name"] == atom_type:
                atom_label = key

        # Read the first two lines separately
        # with open(bin_file, "r") as f:
        #     atom_type = f.readline().split(":")[-1].strip()
        #     atom_label = int(f.readline().split(":")[-1])

        # bins, binned_distances = np.genfromtxt(bin_file, skip_header=2, unpack=True)
        bins, binned_distances = np.genfromtxt(bin_file, skip_header=0, unpack=True)

        bin = pmmoto.analysis.bins.Bin(
            bins[0], bins[-1], len(bins), atom_type, binned_distances
        )

        # Convert to rdf
        _rdf = bin.generate_rdf()

        rdf = pmmoto.domain_generation.rdf.RDF(
            name=bin.name,
            atom_id=atom_label,
            radii=bin.centers,
            rdf=_rdf,
        )

        bounded_rdf = pmmoto.domain_generation.rdf.Bounded_RDF(
            name=bin.name,
            atom_id=atom_label,
            radii=bin.centers,
            rdf=_rdf,
            eps=1.0e-3,
        )

        element = atom_map[atom_label]["element"]
        element_number = pmmoto.particles.convert_atoms_elements_to_ids(
            atom_map[atom_label]["element"]
        )

        radii = pmmoto.particles.uff_radius(atom_names=element)
        equil_radius = radii[element_number[0]] + 1.4

        generate_bin_plot(bin, "data_out/bin_count_plots/")
        generate_rdf_plot(rdf, "data_out/rdf_plots/", equil_radius)
        generate_free_energy_plot(
            bounded_rdf, "data_out/bounded_free_energy_plots/", equil_radius, 5
        )
        generate_bin_and_rdf_plot(bin, rdf, "data_out/bin_count_and_rdf_plots/")


if __name__ == "__main__":
    generate_plots()
