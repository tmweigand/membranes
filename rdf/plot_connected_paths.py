import numpy as np
import matplotlib.pyplot as plt


def read_connected_paths(file_in):
    # Read and parse the file
    results = []

    with open(file_in, "r") as f:
        # Skip the first line ("Connection Results")
        next(f)
        for line in f:
            parts = line.strip().split(", ")
            if len(parts) < 4:
                continue  # skip incomplete lines
            index = int(parts[0])
            file_path = parts[1].split(": ", 1)[1]
            status = parts[2].split(": ", 1)[1]
            result_str = parts[3].split(": ", 1)[1]
            if result_str == "None":
                result = None
            else:
                result = float(result_str)

            results.append(
                {"index": index, "file": file_path, "status": status, "result": result}
            )

    return results


if __name__ == "__main__":
    file_in = "rdf/bridges_results/connected_paths.out"
    results = read_connected_paths(file_in)
    radii = [entry["result"] for entry in results]

    print(len(radii))

    print(np.average(radii), np.median(radii), np.max(radii), np.min(radii))

    # # Plot it
    # fig, ax = plt.subplots()
    # ax.hist(radii, density=True, bins=50, edgecolor="black")
    # # plt.xlim(left=0)
    # ax.set_xlabel("Minimum G(r) for connected path", fontsize=16)
    # ax.set_ylabel("Probability Density", fontsize=16)
    # ax.tick_params(axis="x", labelsize=16)
    # ax.tick_params(axis="y", labelsize=16)
    # # plt.title("Minimum G(r) for Connected Path (kj/mol)")
    # plt.savefig(
    #     "data_out/connected_path_frequency.pdf",
    #     dpi=300,
    #     bbox_inches="tight",
    # )
