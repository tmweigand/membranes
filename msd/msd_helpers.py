import pmmoto
import glob
import re


def get_files(machine):
    """
    Helper to grab the bridges files
    """

    if machine == "bridges":
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

        water_files.sort(key=lambda x: int(re.search(r"(\d+)\.gz$", x).group(1)))

    elif machine == "mac":
        water_files = sorted(glob.glob("data/water_data/*"))

    return water_files


def initialize_domain(voxels):
    """
    Initialize the membrane domain
    """

    if proc_size == 8:
        subdomains = (2, 2, 2)
    elif proc_size == 27:
        subdomains = (3, 3, 3)
    elif proc_size == 64:
        subdomains = (4, 4, 4)
    elif proc_size == 125:
        subdomains = (5, 5, 5)
    elif proc_size == 216:
        subdomains = (6, 6, 6)
    else:
        print("oops")
        subdomains = (0, 0, 0)

    # Full domain with reservoirs
    box = [
        [0.0, 176],
        [0.0, 176],
        [-100, 100],  # Ignore water reservoirs
    ]

    sd = pmmoto.initialize(
        voxels=voxels,
        box=box,
        rank=rank,
        subdomains=subdomains,
        boundary_types=((2, 2), (2, 2), (0, 0)),
        verlet_domains=(20, 20, 20),
        inlet=((0, 0), (0, 0), (1, 0)),
        outlet=((0, 0), (0, 0), (0, 1)),
    )

    return sd
