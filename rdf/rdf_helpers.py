"""rdf_helpers.py"""

import os
import glob
import numpy as np
import re


def check_files(files_1, files_2):
    """
    Ensure the membrane and water files match in length and timestamps.
    """
    if len(files_1) != len(files_2):
        raise ValueError("File lengths do not match.")

    for f1, f2 in zip(files_1, files_2):
        value = -1
        if f1.split(".")[-1] == "gz":
            value = -2
        time_1 = int(f1.split(".")[value])
        time_2 = int(f2.split(".")[value])

        if time_1 != time_2:
            raise ValueError(f"These times dont match! {time_1} {time_2}")


def save_rdf(sd, bins, rdf_bins, atom_labels_to_name, time):
    """
    Save the RDF data to files.

    Parameters
    ----------
    sd : SubDomain
        The subdomain object containing rank information
    bins : Bins
        The bins object containing bin centers
    rdf_bins : dict
        Dictionary containing RDF values for each atom type
    atom_labels_to_name : dict
        Dictionary mapping atom labels to atom names
    """
    if sd.rank == 0:
        # Create output directory if it doesn't exist
        os.makedirs("data_out/rdf_data", exist_ok=True)

        # Save data for each atom type
        for label, rdf in rdf_bins.rdf_bins.items():
            atom_name = atom_labels_to_name[label]["label"]
            centers = bins.bin_centers[label]

            # Stack the data into columns
            data = np.column_stack((centers, rdf))

            # Save with header
            np.savetxt(
                f"data_out/rdf_data/rdf_{atom_name}.txt",
                data,
                header=f"Distance(Å) RDF at {time}",
                delimiter="\t",
                comments="#",
            )


def get_bridges_files():
    """
    Helper to grab the bridges files
    """
    membrane_files = glob.glob(
        "/ocean/projects/cts200024p/rvickers/RV-P3/64x/03_90/perm_v2/membrane/*"
    )
    water_files = glob.glob(
        "/ocean/projects/cts200024p/rvickers/RV-P3/64x/03_90/perm_v2/pressure/*"
    )

    # Remove specific unwanted files
    unwanted_files = [
        "/ocean/projects/cts200024p/rvickers/RV-P3/64x/03_90/perm_v2/pressure/pressuredata.100000000.gz",
        "/ocean/projects/cts200024p/rvickers/RV-P3/64x/03_90/perm_v2/pressure/pressuredata.110005000.gz",
        "/ocean/projects/cts200024p/rvickers/RV-P3/64x/03_90/perm_v2/pressure/pressuredata.99995000.gz",
    ]
    for file in unwanted_files:
        if file in water_files:
            water_files.remove(file)

    # Checked and working
    membrane_files.sort(key=lambda x: int(re.search(r"(\d+)\.gz$", x).group(1)))
    water_files.sort(key=lambda x: int(re.search(r"(\d+)\.gz$", x).group(1)))

    return membrane_files, water_files
