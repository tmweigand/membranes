"""msd_analysis"""

import pmmoto
from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt
import msd_helpers

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
proc_size = comm.Get_size()


def main():

    # water_files = msd_helpers.get_files(machine="mac")
    water_files = msd_helpers.get_files(machine="mac", file_type="gzip")

    msd_values = []
    times = []

    water_files = water_files * 10

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

            print(f"Time {time}: MSD = {msd} file {file}")

    # Convert to arrays
    msd_values_arr = np.array(msd_values)
    times_arr = np.array(times) - np.min(times)

    # Save results to CSV
    output_file = "msd/msd_results.csv"
    np.savetxt(
        output_file,
        np.column_stack([times_arr, msd_values_arr]),
        delimiter=",",
        header="time,msd",
        comments="",
    )
    print(f"Saved MSD results to {output_file}")

    # Plot MSD vs time
    plt.figure(figsize=(6, 4))
    plt.plot(times_arr, msd_values_arr, "--")
    plt.xlabel("Time")
    plt.ylabel("Mean Squared Displacement (MSD)")
    plt.title("MSD vs Time")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("msd/msd_plot.png", dpi=300)


if __name__ == "__main__":
    import cProfile
    import pstats

    profiler = cProfile.Profile()
    profiler.enable()
    main()
    profiler.disable()
    stats = pstats.Stats(profiler).sort_stats("cumtime")
    stats.print_stats(20)  # Print top 20 slowest functions
