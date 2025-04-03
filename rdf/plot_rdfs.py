"""plot_rdfs.py"""

import glob
import numpy as np
import pmmoto
import matplotlib.pyplot as plt


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


def generate_rdf_plot(bin, rdf, folder, x_line=None, y_line=None):
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
    plt.plot(bin.centers, rdf, linewidth=2, color="navy", label="RDF")

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
    plt.title(f"RDF for Atom Type {bin.name}", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.xlim(left=0)
    plt.ylim(bottom=0)

    plt.legend()
    plt.tight_layout()

    out_file = folder + f"bin_count_{bin.name}.pdf"
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close()


def generate_free_energy_plot(bin, rdf, folder, x_line=None, y_line=None):
    """
    Generate and save an RDF plot.

    Parameters:
    - bin: Object containing bin information, including bin centers and name.
    - rdf: Array-like, radial distribution function values.
    - folder: String, directory where the plot will be saved.
    - x_line: Float (optional), location on the x-axis to draw a vertical line.
    """

    k_b = 0.0083144621
    T = 300

    with np.errstate(divide="ignore"):
        G = -k_b * T * np.log(rdf)

    # G = G / np.sum(rdf)

    # Normalize to set minimum free energy to zero
    # G -= np.nanmin(G)

    pmmoto.io.io_utils.check_file_path(folder)

    plt.figure(figsize=(8, 5))
    plt.plot(bin.centers, G, linewidth=2, color="navy", label="Simulation")

    # Draw vertical line if x_line is provided
    if x_line is not None:
        plt.axvline(
            x=x_line,
            color="red",
            linestyle="--",
            linewidth=1.5,
            label=f"VdW = {x_line:.2f}",
        )

    if y_line is not None:
        plt.axhline(
            y=y_line,
            color="blue",
            linestyle="--",
            linewidth=1.5,
            label=f"G = {y_line:.2f}",
        )

    plt.xlabel("Distance (Å)", fontsize=12)
    plt.ylabel("G(r) kJ/mol", fontsize=12)
    plt.title(f"Free Energy for Atom Type {bin.name}", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.xlim(left=0)
    # plt.ylim(top=20)

    plt.legend()
    plt.tight_layout()

    out_file = folder + f"G_{bin.name}.pdf"
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
            if values["label"] == atom_type:
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
        rdf = bin.generate_rdf()

        element = atom_map[atom_label]["element"]
        element_number = pmmoto.particles.convert_atoms_elements_to_ids(
            atom_map[atom_label]["element"]
        )

        radii = pmmoto.particles.uff_radius(atom_names=element)
        equil_radius = radii[element_number[0]] + 1.4

        # SAve rdf files
        data = np.column_stack((bins, rdf))

        # Save with header
        folder = "data_out/rdf_data/"
        pmmoto.io.io_utils.check_file_path(folder)
        out_file = folder + f"{atom_type}.rdf"

        # header = f"Atom Type: {bin.name} \nAtom Label: {label}"
        np.savetxt(out_file, data, delimiter="\t")

        generate_bin_plot(bin, "data_out/bin_count_plots/")
        generate_rdf_plot(bin, rdf, "data_out/rdf_plots/", equil_radius)
        generate_free_energy_plot(
            bin, rdf, "data_out/free_energy_plots/", equil_radius, 5
        )


if __name__ == "__main__":
    generate_plots()
