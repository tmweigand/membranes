import numpy as np
import glob
import matplotlib.pyplot as plt


def get_time(file):
    try:
        return int(file.split(".")[-2])
    except ValueError:
        print("Check File name for the time. ")


def read_lammps_data(filename):

    with open(filename, "r") as f:
        lines = f.readlines()

    # skip to atoms
    for i, line in enumerate(lines):
        if line.startswith("ITEM: ATOMS"):
            start_idx = i + 1
            break

    # Read the ATOMS section into a NumPy array
    data = np.loadtxt(lines[start_idx:], dtype=int)  # Convert to integer array

    return data


files = glob.glob("./bridges_runs/xlink_95/out/tim_mol.*")
files = sorted(files, key=get_time)

time = 0
for file in files:
    new_time = int(file.split(".")[-2])
    assert new_time > time
    time = new_time

num_files = len(files)
molecule_info = np.zeros([num_files, 3])
for n, file in enumerate(files):
    data = read_lammps_data(file)
    num_atoms = data.shape[0]
    molecule_info[n, 0] = int(file.split(".")[-2])
    molecule_info[n, 1] = np.max(data[:, 1])
    molecule_info[n, 2] = np.max(data[:, 2])

    bins = np.linspace(0, 300, 30)
    mol_sizes = np.unique(data[:, 2], axis=0)
    plt.figure(figsize=(8, 6))
    plt.hist(mol_sizes, bins, edgecolor="black", alpha=0.7)
    plt.xlabel("Size of Molecules")
    plt.ylabel("Frequency")
    plt.savefig(f"data_out/small_mols/mol_sizes_{n}.pdf")
    plt.close()

    bins = np.linspace(0, 17500, 500)
    plt.figure(figsize=(8, 6))
    plt.hist(mol_sizes, bins, edgecolor="black", alpha=0.7)
    plt.xlabel("Size of Molecules")
    plt.ylabel("Frequency")
    plt.savefig(f"data_out/all_mols/mol_sizes_{n}.pdf")
    plt.close()

plt.figure(figsize=(8, 6))
plt.semilogx(molecule_info[:, 0], molecule_info[:, 1])
plt.xlabel("Simulation Time")
plt.ylabel("Number of Molecules")
plt.savefig(f"data_out/number_molecules.pdf")
plt.close()

plt.figure(figsize=(8, 6))
plt.semilogx(molecule_info[:, 0], molecule_info[:, 2])
plt.xlabel("Simulation Time")
plt.ylabel("Max Atoms per Molecule")
plt.savefig(f"data_out/max_atoms_in_molecule.pdf")
plt.close()
