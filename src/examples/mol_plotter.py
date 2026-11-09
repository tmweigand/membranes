import sys

sys.path.insert(0, "/Users/tim/Desktop/membranes/code/membranes/src")

from membranes.monomers.molecule_plotter import MoleculePlotter


# Parse type labels from the ff.lammps file.
# Supports two formats:
#   1. "Atom Type Labels" section:  1  ca
#   2. Masses section with comments: 1 12.011  # ca
def read_type_labels(ff_path):
    labels = {}

    # Pass 1: look for "Atom Type Labels" section
    in_section = False
    with open(ff_path) as f:
        for line in f:
            line = line.strip()
            if line == "Atom Type Labels":
                in_section = True
                continue
            if in_section:
                if not line:
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        labels[int(parts[0])] = parts[1]
                    except ValueError:
                        break
                else:
                    break

    if labels:
        return labels

    # Pass 2: fall back to Masses section comments:  1  12.011  # ca
    in_masses = False
    with open(ff_path) as f:
        for line in f:
            stripped = line.strip()
            if stripped == "Masses":
                in_masses = True
                continue
            if in_masses:
                if not stripped:
                    continue
                # stop at next section header (non-numeric first token)
                parts = stripped.split()
                try:
                    tid = int(parts[0])
                except ValueError:
                    break
                if "#" in stripped:
                    name = stripped.split("#", 1)[1].strip().split()[0]
                    labels[tid] = name

    return labels


type_labels = read_type_labels("data_out/monomeric_states/ff.lammps")
print(type_labels)
MoleculePlotter.plot_mol_file(
    "data_out/monomeric_states/pre_reaction.mol",
    type_labels=type_labels,
    out_png="data_out/monomeric_states/pre_reaction.png",
)

MoleculePlotter.plot_mol_file(
    "data_out/monomeric_states/post_reaction.mol",
    type_labels=type_labels,
    out_png="data_out/monomeric_states/post_reaction.png",
)
