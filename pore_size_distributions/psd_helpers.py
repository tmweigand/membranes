"""Helper functions for PSD"""

import glob


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

    return membrane_files, water_files
