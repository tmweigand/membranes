"""Generate MPD/TMC pre- and post-reaction LAMMPS molecule templates."""

import os
from dataclasses import replace as dc_replace

from membranes.monomers import generate_molecule_data
from membranes.monomers.monomer_data import Bond, Angle, Dihedral
from membranes.monomers.molecule import Molecule

data_out = "data_out/monomeric_states"
os.makedirs(data_out, exist_ok=True)
os.makedirs(f"{data_out}/run", exist_ok=True)


def find_atom(structure, atom_type, bonded_to=None):
    """Return 1-based ID of first atom matching type; optionally bonded to bonded_to type."""
    for a in structure.atoms:
        if a.type != atom_type:
            continue
        if bonded_to is None or any(n.type == bonded_to for n in a.bond_partners):
            return a.idx + 1
    raise RuntimeError(
        f"No '{atom_type}'" + (f" bonded to '{bonded_to}'" if bonded_to else "")
    )


def neighborhood(structure, center_id, max_bonds=3):
    """Return sorted 1-based atom IDs within max_bonds graph hops of center_id."""
    center = structure.atoms[center_id - 1]
    seen = {center}
    frontier = [center]
    for _ in range(max_bonds):
        nxt = [n for a in frontier for n in a.bond_partners if n not in seen]
        seen.update(nxt)
        frontier = nxt
    return sorted(a.idx + 1 for a in seen)


def slice_mol(mol, keep_ids, name, mol_id):
    """Extract atoms in keep_ids, reindex from 1. Returns (Molecule, old_id -> new_id map)."""
    keep = set(keep_ids)
    ordered = [a for a in mol.atoms if a.id in keep]
    rm = {a.id: i + 1 for i, a in enumerate(ordered)}

    def remap(term, i):
        return dc_replace(term, id=i + 1, atom_ids=tuple(rm[x] for x in term.atom_ids))

    atoms = [dc_replace(a, id=rm[a.id], mol=mol_id) for a in ordered]
    bonds = [
        remap(t, i)
        for i, t in enumerate(
            t for t in mol.bonds if all(x in keep for x in t.atom_ids)
        )
    ]
    angles = [
        remap(t, i)
        for i, t in enumerate(
            t for t in mol.angles if all(x in keep for x in t.atom_ids)
        )
    ]
    dihedrals = [
        remap(t, i)
        for i, t in enumerate(
            t for t in mol.dihedrals if all(x in keep for x in t.atom_ids)
        )
    ]
    impropers = [
        remap(t, i)
        for i, t in enumerate(
            t for t in mol.impropers if all(x in keep for x in t.atom_ids)
        )
    ]
    return Molecule(name, atoms, bonds, angles, dihedrals, impropers), rm


def merge(name, frags):
    """Concatenate fragment Molecules into one, offsetting all atom IDs."""
    all_atoms, all_bonds, all_angles, all_dihedrals, all_impropers = [], [], [], [], []
    offset = 0
    for frag_idx, frag in enumerate(frags, 1):
        all_atoms.extend(
            dc_replace(a, id=a.id + offset, mol=frag_idx) for a in frag.atoms
        )
        for lst, src in [
            (all_bonds, frag.bonds),
            (all_angles, frag.angles),
            (all_dihedrals, frag.dihedrals),
            (all_impropers, frag.impropers),
        ]:
            for t in src:
                lst.append(
                    dc_replace(
                        t,
                        id=len(lst) + 1,
                        atom_ids=tuple(x + offset for x in t.atom_ids),
                    )
                )
        offset += frag.n_atoms
    return Molecule(
        name, all_atoms, all_bonds, all_angles, all_dihedrals, all_impropers
    )


def bonded_ids(mol, atom_id):
    """Return IDs of atoms directly bonded to atom_id."""
    result = []
    for b in mol.bonds:
        a, c = b.atom_ids
        if a == atom_id:
            result.append(c)
        elif c == atom_id:
            result.append(a)
    return result


def build_adj(mol):
    """Adjacency dict: atom_id -> [neighbor_ids]."""
    adj = {a.id: [] for a in mol.atoms}
    for b in mol.bonds:
        adj[b.atom_ids[0]].append(b.atom_ids[1])
        adj[b.atom_ids[1]].append(b.atom_ids[0])
    return adj


def bfs_order(adj, types, root, max_bonds=3):
    """BFS from root, neighbours sorted by atom type at each level. Returns ordered ID list."""
    seen, frontier, order = {root}, [root], [root]
    for _ in range(max_bonds):
        nxt = []
        for aid in frontier:
            for nei in sorted(adj.get(aid, []), key=lambda x: types.get(x, "")):
                if nei not in seen:
                    seen.add(nei)
                    nxt.append(nei)
                    order.append(nei)
        frontier = nxt
    return order


# ── generate monomers ──────────────────────────────────────────────────────

oligomers = {
    "mpd0": "Nc1cccc(N)c1",
    "tmc0": "ClC(=O)c1cc(C(Cl)=O)cc(C(Cl)=O)c1",
    "mpd1_tmc1": "Nc1cccc(NC(=O)c2cc(C(Cl)=O)cc(C(Cl)=O)c2)c1",
    # "mpd12_tmc2": "Nc1cccc(NC(=O)c2cc(C(Cl)=O)cc(C(=O)Nc3cccc(N)c3)c2)c1",
    # "mpd13_tmc3": "Nc1cccc(NC(=O)c2cc(C(=O)Nc3cccc(N)c3)cc(C(=O)Nc3cccc(N)c3)c2)c1",
    # "mpd2_tmc12": "O=C(Nc1cccc(NC(=O)c2cc(C(Cl)=O)cc(C(Cl)=O)c2)c1)c1cc(C(Cl)=O)cc(C(Cl)=O)c1",
    # "mpd2_tmc22": (
    #     "O=C(Nc1cccc(NC(=O)c2cc(C(Cl)=O)cc(C(=O)Nc3cccc(N)c3)c2)c1)"
    #     "c1cc(C(Cl)=O)cc(C(=O)Nc2cccc(N)c2)c1"
    # ),
    # "mpd2_tmc23": (
    #     "O=C(Nc1cccc(NC(=O)c2cc(C(=O)Nc3cccc(N)c3)cc(C(=O)Nc3cccc(N)c3)c2)c1)"
    #     "c1cc(C(=O)Nc2cccc(N)c2)cc(C(=O)Nc2cccc(N)c2)c1"
    # ),
}
generated = {}

for mol_name, smiles in oligomers.items():
    g = generate_molecule_data.GenMolecule(
        mol_name,
        smiles,
        charge=0,
        forcefield="gaff2",
        charge_model="bcc",
        outdir=f"{data_out}/run",
    )
    g.generate_forcefield()
    g.create_pmd()
    generated[mol_name] = g
    g.plot_2d_parameters(kind="charge", out_png=f"{data_out}/{mol_name}_charge.png")

mpd, tmc = generated["mpd0"], generated["tmc0"]
mpd_mol = mpd.convert_to_molecule(mol_id=1)
tmc_mol = tmc.convert_to_molecule(mol_id=2)

# ── reactive sites & local neighborhoods ──────────────────────────────────

nv_id = find_atom(mpd.pmd_structure, "nv")
c_id = find_atom(tmc.pmd_structure, "c", bonded_to="cl")

mpd_ids = neighborhood(mpd.pmd_structure, nv_id)
tmc_ids = neighborhood(tmc.pmd_structure, c_id)
print(f"MPD local ({len(mpd_ids)} atoms): {mpd_ids}")
print(f"TMC local ({len(tmc_ids)} atoms): {tmc_ids}")

mpd_frag, mpd_rm = slice_mol(mpd_mol, mpd_ids, "mpd_local", mol_id=1)
tmc_frag, tmc_rm = slice_mol(tmc_mol, tmc_ids, "tmc_local", mol_id=2)

# ── pre-reaction ───────────────────────────────────────────────────────────

pre = merge("mpd_tmc_pre_reaction", [mpd_frag, tmc_frag])
pre.write_mol(f"{data_out}/pre_reaction.mol")
print(f"Wrote pre_reaction.mol  ({pre.n_atoms} atoms, {pre.n_bonds} bonds)")

# ── post-reaction: nv-c bond formed, one hn + cl deleted ──────────────────

# IDs of nv and c inside the merged mol
nv_pre = mpd_rm[nv_id]
c_pre = mpd_frag.n_atoms + tmc_rm[c_id]

atom_type = {a.id: a.type for a in pre.atoms}

hn_del = next(x for x in bonded_ids(pre, nv_pre) if atom_type[x] == "hn")
cl_del = next(x for x in bonded_ids(pre, c_pre) if atom_type[x] == "cl")
drop = {hn_del, cl_del}

keep_post = [a.id for a in pre.atoms if a.id not in drop]
keep_set = set(keep_post)
rm_post = {old: new for new, old in enumerate(keep_post, 1)}

new_nv = rm_post[nv_pre]
new_c = rm_post[c_pre]
atype_post = {rm_post[old]: atom_type[old] for old in keep_post}


def remap_post(term, i):
    return dc_replace(term, id=i + 1, atom_ids=tuple(rm_post[x] for x in term.atom_ids))


atoms_p = [
    dc_replace(a, id=rm_post[a.id], mol=1) for a in pre.atoms if a.id in keep_set
]
bonds_p = [
    remap_post(t, i)
    for i, t in enumerate(
        t for t in pre.bonds if all(x in keep_set for x in t.atom_ids)
    )
]
angles_p = [
    remap_post(t, i)
    for i, t in enumerate(
        t for t in pre.angles if all(x in keep_set for x in t.atom_ids)
    )
]
dihedrals_p = [
    remap_post(t, i)
    for i, t in enumerate(
        t for t in pre.dihedrals if all(x in keep_set for x in t.atom_ids)
    )
]
impropers_p = [
    remap_post(t, i)
    for i, t in enumerate(
        t for t in pre.impropers if all(x in keep_set for x in t.atom_ids)
    )
]

# new nv-c bond
bonds_p.append(
    Bond(
        id=len(bonds_p) + 1,
        type=("nv", "c"),
        atom_ids=(new_nv, new_c),
        atom_names=("nv", "c"),
        k=0.0,
        r_eq=0.0,
    )
)

# neighbors of nv/c in post (original bonds, minus dropped atoms)
nv_nbrs = [rm_post[x] for x in bonded_ids(pre, nv_pre) if x in keep_set]
c_nbrs = [rm_post[x] for x in bonded_ids(pre, c_pre) if x in keep_set]

# angles through new bond: X-nv-c  and  nv-c-Y
for x in nv_nbrs:
    angles_p.append(
        Angle(
            id=len(angles_p) + 1,
            type=(atype_post[x], "nv", "c"),
            atom_ids=(x, new_nv, new_c),
            atom_names=(atype_post[x], "nv", "c"),
            k=0.0,
            theta_eq=0.0,
        )
    )
for y in c_nbrs:
    angles_p.append(
        Angle(
            id=len(angles_p) + 1,
            type=("nv", "c", atype_post[y]),
            atom_ids=(new_nv, new_c, y),
            atom_names=("nv", "c", atype_post[y]),
            k=0.0,
            theta_eq=0.0,
        )
    )

# dihedrals through new bond: X-nv-c-Y
for x in nv_nbrs:
    for y in c_nbrs:
        dihedrals_p.append(
            Dihedral(
                id=len(dihedrals_p) + 1,
                type=(atype_post[x], "nv", "c", atype_post[y]),
                atom_ids=(x, new_nv, new_c, y),
                atom_names=(atype_post[x], "nv", "c", atype_post[y]),
                k=0.0,
                n=2,
                d=0,
                weight=1.0,
            )
        )

post = Molecule(
    "mpd_tmc_post_reaction", atoms_p, bonds_p, angles_p, dihedrals_p, impropers_p
)

# ── update post-reaction charges from mpd1_tmc1 ───────────────────────────

ref_gen = generated["mpd1_tmc1"]
ref_mol_full = ref_gen.convert_to_molecule(mol_id=1)

# find amide N in mpd1_tmc1: any N bonded to carbonyl c (not the unreacted nv)
amide_n_ref = next(
    a.idx + 1
    for a in ref_gen.pmd_structure.atoms
    if a.type != "nv" and any(n.type == "c" for n in a.bond_partners)
)

ref_types_map = {a.id: a.type for a in ref_mol_full.atoms}
ref_charge_map = {a.id: a.charge for a in ref_mol_full.atoms}
ref_bfs = bfs_order(build_adj(ref_mol_full), ref_types_map, amide_n_ref, max_bonds=3)

post_types_map = {a.id: a.type for a in post.atoms}
post_bfs = bfs_order(build_adj(post), post_types_map, new_nv, max_bonds=3)

n_match = min(len(post_bfs), len(ref_bfs))
charge_update = {post_bfs[i]: ref_charge_map[ref_bfs[i]] for i in range(n_match)}
print(f"Updated {n_match} atom charges from mpd1_tmc1:")
for pid, rid in zip(post_bfs[:n_match], ref_bfs[:n_match]):
    print(f"  atom {pid} ({post_types_map[pid]}) {ref_charge_map[rid]:+.4f}")

updated_atoms = [
    dc_replace(a, charge=charge_update.get(a.id, a.charge)) for a in post.atoms
]
post = Molecule(
    post.name, updated_atoms, post.bonds, post.angles, post.dihedrals, post.impropers
)

post.write_mol(f"{data_out}/post_reaction.mol")
print(f"Wrote post_reaction.mol ({post.n_atoms} atoms, {post.n_bonds} bonds)")
