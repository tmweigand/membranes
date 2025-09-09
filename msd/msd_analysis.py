"""msd_analysis"""

import pmmoto
import glob
from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
proc_size = comm.Get_size()


def get_bridges_files():
    """
    Helper to grab the bridges files
    """

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


if __name__ == "__main__":

    bridges = False

    voxels_in = (800, 800, 800)
    if bridges:
        water_files = get_bridges_files()
    else:
        water_files = sorted(glob.glob("data/water_data/*"))

    msd_values = []
    times = []

    initial_time = True
    for file in water_files:
        ids, positions, types, domain, time = pmmoto.io.data_read.read_lammps_atoms(
            file
        )

        order = np.argsort(ids)
        ids, positions = ids[order], positions[order]

        if initial_time:
            num_atoms = len(ids)
            initial_positions = positions
            initial_ids = ids
            initial_time = False
        else:
            assert np.array_equal(ids, initial_ids)

            disp = positions - initial_positions
            sq_disp = np.sum(disp**2, axis=1)
            msd = np.mean(sq_disp)

            msd_values.append(msd)
            times.append(time)

            print(f"Time {time}: MSD = {msd}")


# Convert to arrays
msd_values = np.array(msd_values)
times = np.array(times)


# Save results to CSV
output_file = "msd/msd_results.csv"
np.savetxt(
    output_file,
    np.column_stack([times, msd_values]),
    delimiter=",",
    header="time,msd",
    comments="",
)
print(f"Saved MSD results to {output_file}")

# Plot MSD vs time
plt.figure(figsize=(6, 4))
plt.plot(times, msd_values, "--")
plt.xlabel("Time")
plt.ylabel("Mean Squared Displacement (MSD)")
plt.title("MSD vs Time")
plt.grid(True)
plt.tight_layout()
plt.savefig("msd/msd_plot.png", dpi=300)
