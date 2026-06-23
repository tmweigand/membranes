"""polymerization.py

Defines the PolymerizationSimulation class for TMC-MPD interfacial polymerization.
"""

import os
import numpy as np
from lammps import LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR

from membranes.domain_generation.lammps_init import LAMMPSInitialize
from membranes.domain_generation.style_properties import StyleProperties
from membranes.domain_generation.system_properties import SystemProperties


def random_int() -> int:
    return np.random.randint(1, 32768)


class PolymerizationSimulation:
    """
    Encapsulates a TMC-MPD interfacial polymerization simulation in LAMMPS.

    Stages
    ------
    1. pack_molecules  — place MPD and TMC into the simulation box
    2. minimize        — steepest-descent then conjugate-gradient minimization
    3. polymerize_pa   — amide (PA) bond formation via bond/react
    4. add_hydroxide   — insert OH molecules to cap unreacted chloride sites
    5. polymerize_oh   — TMC–OH bond formation via bond/react
    6. cleanup         — remove excess atoms and write final data file
    """

    # Atom counts per molecule
    MPD_ATOMS_PER_MOL = 16
    TMC_ATOMS_PER_MOL = 18
    MPD_TMC_RATIO = 3.0 / 2.0

    # Atom type indices
    AMIDE_TYPE = 7
    CL_TYPE = 18
    EXCESS_TYPES = (11, 13, 18)

    def __init__(
        self,
        in_dir: str = "rv",
        multiple: float = 0.1,
        xlink: float = 0.80,
        temperature: float = 300,
        bond_frequency: int = 50,
        bond_distance: tuple[float, float] = (0.0, 5.0),
        stabilization: float = 0.03,
        max_cycles: int = 375,
    ):
        self.in_dir = in_dir
        self.multiple = multiple
        self.xlink = xlink
        self.temperature = temperature
        self.bond_frequency = bond_frequency
        self.bond_distance = bond_distance
        self.stabilization = stabilization
        self.max_cycles = max_cycles

        self.out_dir = f"data_out/{in_dir}"
        os.makedirs(self.out_dir, exist_ok=True)

        self._init_lammps()

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------

    def _init_lammps(self):
        wrapper = LAMMPSInitialize(
            log_file="log.polymerize", units="real", atom_style="full"
        )
        self.lmp = wrapper.lmp

        if self.in_dir == "rv":
            StyleProperties(
                lmp=self.lmp,
                dihedral_style="charmm",
                improper_style="harmonic",
                special_bonds="charmm",
            )
        else:
            StyleProperties(lmp=self.lmp)

        self.system = SystemProperties(
            self.lmp, dimension=3, n_atom_types=30, boundary=("p", "p", "p")
        )

    def _data_path(self, filename: str) -> str:
        return f"data_in/{self.in_dir}/{filename}"

    # ------------------------------------------------------------------
    # Stage 1 — molecule packing
    # ------------------------------------------------------------------

    def pack_molecules(self):
        """Read FF data, place MPD and TMC molecules into the box."""
        box_size = 90 * self.multiple ** (1 / 3)
        num_mpd = int(300 * self.multiple)
        num_tmc = int(num_mpd / self.MPD_TMC_RATIO)
        target_atoms = (
            self.TMC_ATOMS_PER_MOL * num_tmc + self.MPD_ATOMS_PER_MOL * num_mpd
        )

        # Load force field
        if self.in_dir == "rv":
            self.lmp.cmd.read_data(
                self._data_path("ff.lammps"), "extra/improper/per/atom", 2
            )
        else:
            self.lmp.cmd.read_data(self._data_path("ff.lammps"))

        self.lmp.cmd.molecule("MPD", self._data_path("lammps_mpd.mol"))
        self.lmp.cmd.molecule("TMC", self._data_path("lammps_tmc.mol"))

        # Resize box
        self.lmp.cmd.change_box(
            "all",
            "x",
            "final",
            0,
            box_size,
            "y",
            "final",
            0,
            box_size,
            "z",
            "final",
            0,
            box_size,
            "units",
            "box",
        )
        self.lmp.cmd.delete_atoms("group", "all")

        # Place MPD
        self.lmp.cmd.lattice("fcc", 11.3)
        self.lmp.cmd.create_atoms(
            0, "box", "subset", num_mpd, random_int(), "mol", "MPD", random_int()
        )
        self.lmp.cmd.group("MPDs", "union", "all")

        # Iteratively place TMC, removing overlaps each round
        for p in range(1, 1001):
            current_atoms = self.lmp.get_natoms()
            num_tmc_needed = int(
                (target_atoms - current_atoms) / self.TMC_ATOMS_PER_MOL
            )
            if num_tmc_needed <= 0:
                print("All TMC molecules placed")
                break

            self.lmp.cmd.create_atoms(
                0,
                "box",
                "subset",
                num_tmc_needed,
                random_int(),
                "mol",
                "TMC",
                random_int(),
            )
            self.lmp.cmd.group("TMCs", "subtract", "all", "MPDs")
            self.lmp.cmd.delete_atoms("overlap", 1.0, "TMCs", "all", "mol", "yes")

            self.lmp.cmd.fix("1", "all", "nve")
            self.lmp.cmd.run(0)
            self.lmp.cmd.unfix("1")

            current_atoms = self.lmp.get_natoms()
            print(f"  Loop {p}: {current_atoms}/{target_atoms} atoms placed")
            if current_atoms >= target_atoms:
                break

        self.lmp.cmd.reset_atoms("id")

        # Dump packing structure
        self._dump("2", "logs/packing_structure.lammpstrj", freq=1000)
        self.lmp.cmd.run(0)
        self.lmp.cmd.undump("2")

    # ------------------------------------------------------------------
    # Stage 2 — minimization
    # ------------------------------------------------------------------

    def minimize(self):
        """Two-stage minimization: steepest descent then conjugate gradient."""
        self.lmp.cmd.dielectric(1.0)
        self.lmp.cmd.neighbor(2.0, "bin")
        self.lmp.cmd.neigh_modify("delay", 0, "every", 1, "check", "yes")
        self.lmp.cmd.timestep(10.0)
        self.lmp.cmd.run_style("verlet")

        self.lmp.cmd.min_style("sd")
        self.lmp.cmd.minimize(1.0e-5, 1.0e-5, 10000, 100000)

        self.lmp.cmd.min_style("cg")
        self.lmp.cmd.min_modify("line", "quadratic")
        self.lmp.cmd.minimize(1.0e-8, 1.0e-8, 10000, 100000)

        self.lmp.cmd.write_data(f"{self.out_dir}/data.lmps")
        self.lmp.cmd.timestep(1.0)

    # ------------------------------------------------------------------
    # Stage 3 — PA polymerization
    # ------------------------------------------------------------------

    def polymerize_pa(self) -> float:
        """
        First cross-linking stage: TMC–MPD amide bond formation.

        Returns
        -------
        max_bonds : float
            Total number of available amide bonding sites (used to size the
            hydroxide addition in the next stage).
        """
        self.lmp.cmd.molecule(
            "pre_reaction_TMC_MPD", self._data_path("mpd_tmc_pre_reaction.mol")
        )
        self.lmp.cmd.molecule(
            "post_reaction_TMC_MPD", self._data_path("mpd_tmc_post_reaction.mol")
        )

        self.lmp.cmd.group("amide", "type", self.AMIDE_TYPE)
        self.lmp.cmd.variable("n_amides", "equal", "count(amide)")
        max_bonds = self.lmp.extract_variable("n_amides")
        target_bonds = max_bonds * self.xlink
        print(f"MAX_BONDS = {max_bonds}  |  TARGET = {target_bonds}")

        self._setup_bond_react(
            "PArxn",
            "pre_reaction_TMC_MPD",
            "post_reaction_TMC_MPD",
            "mpd_tmx_rxnmap.in",
        )
        self.lmp.cmd.thermo(100)
        self.lmp.cmd.thermo_style(
            "custom", "step", "temp", "press", "density", "f_PArxn[1]"
        )

        bonds_made = self._polymerize_loop("PArxn", target_bonds, label="NVT PA")
        print(f"Total PA bonds formed = {bonds_made}")
        self.lmp.cmd.unfix("PArxn")

        # Equilibrate
        self.lmp.cmd.velocity("all", "scale", self.temperature)
        self.lmp.cmd.fix(
            1, "all", "nvt", "temp", self.temperature, self.temperature, 100.0
        )
        self._dump(
            3, "logs/polymerization_structure.lammpstrj", freq=self.bond_frequency
        )

        print("PA EQUILIBRATION RUN")
        self.lmp.cmd.thermo_style("custom", "step", "temp", "press", "density")
        self.lmp.cmd.run(25000, "upto")
        self.lmp.cmd.unfix(1)
        self.lmp.cmd.undump(3)

        self.lmp.cmd.reset_atoms("id")
        self.lmp.cmd.write_data("logs/polym_clean.lmps")

        return max_bonds

    # ------------------------------------------------------------------
    # Stage 4 — add hydroxide
    # ------------------------------------------------------------------

    def add_hydroxide(self, max_bonds: float):
        """Insert OH molecules to cap unreacted chloride sites."""
        num_oh = int(max_bonds * 0.8)
        print(f"NUMBER of OH: {num_oh}")

        self.lmp.cmd.region(
            "new_box", "block", "INF", "INF", "INF", "INF", "INF", "INF"
        )
        self.lmp.cmd.molecule("OH", self._data_path("hydroxide.mol"))
        self.lmp.cmd.create_atoms(
            0, "random", num_oh, random_int(), "new_box", "mol", "OH", random_int()
        )
        self.lmp.cmd.minimize(1.0e-4, 1.0e-4, 1000, 100000)
        self.lmp.cmd.write_data("logs/term_data.lmps")

    # ------------------------------------------------------------------
    # Stage 5 — OH polymerization
    # ------------------------------------------------------------------

    def polymerize_oh(self):
        """Second cross-linking stage: TMC–OH bond formation."""
        self.lmp.cmd.molecule("mol3", self._data_path("TMC_OH_prerxn.mol"))
        self.lmp.cmd.molecule("mol4", self._data_path("TMC_OH_postrxn.mol"))

        self._setup_bond_react("OHrxn", "mol3", "mol4", "TMC_OH_rxnmap.in")

        self.lmp.cmd.group("Cl", "type", self.CL_TYPE)
        self.lmp.cmd.variable("n_Cl", "equal", "count(Cl)")
        open_cl = self.lmp.extract_variable("n_Cl")
        print(f"OPEN CHLORIDES: {open_cl}")

        self.lmp.cmd.thermo_style(
            "custom", "step", "temp", "press", "density", "f_OHrxn[1]"
        )

        bonds_made = self._polymerize_loop("OHrxn", open_cl, label="OH")
        print(f"Total OH bonds formed = {bonds_made}")
        self.lmp.cmd.unfix("OHrxn")

        # Final equilibration
        self.lmp.cmd.fix(
            1, "all", "nvt", "temp", self.temperature, self.temperature, 100.0
        )
        self._dump(4, "logs/OH_structure.lammpstrj", freq=10000)

        print("FINAL EQUILIBRATION RUN")
        self.lmp.cmd.thermo_style("custom", "step temp press density")
        self.lmp.cmd.run(50000, "upto")
        self.lmp.cmd.undump(4)
        self.lmp.cmd.unfix(1)

    # ------------------------------------------------------------------
    # Stage 6 — cleanup
    # ------------------------------------------------------------------

    def cleanup(self):
        """Remove excess atoms and write the final data file."""
        excess_types = " ".join(str(t) for t in self.EXCESS_TYPES)
        self.lmp.cmd.group("EXCESSOHCL", "type", excess_types)
        self.lmp.cmd.delete_atoms("group", "EXCESSOHCL", "bond", "yes")
        self.lmp.cmd.reset_atoms("id")
        self.lmp.cmd.write_data("logs/term_final.lmps")

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _setup_bond_react(
        self, fix_name: str, pre_mol: str, post_mol: str, rxnmap: str
    ):
        """Register a bond/react fix."""
        self.lmp.cmd.fix(
            fix_name,
            "all",
            "bond/react",
            "stabilization",
            "yes",
            "statted_grp",
            self.stabilization,
            "react",
            "rxn1",
            "all",
            self.bond_frequency,
            self.bond_distance[0],
            self.bond_distance[1],
            pre_mol,
            post_mol,
            self._data_path(rxnmap),
        )

    def _run_nvt_npt_cycle(self):
        """One NVT + NPT cycle used inside each polymerization loop."""
        T = self.temperature

        self.lmp.cmd.velocity("all", "create", T, random_int())

        self.lmp.cmd.fix(1, "statted_grp_REACT", "nvt", "temp", T, T, 100.0)
        self.lmp.cmd.fix(4, "bond_react_MASTER_group", "temp/rescale", 50, T, T, 10, 1)
        self.lmp.cmd.run(1000)
        self.lmp.cmd.unfix(1)
        self.lmp.cmd.unfix(4)

        self.lmp.cmd.fix(
            1,
            "statted_grp_REACT",
            "npt",
            "temp",
            T,
            T,
            100.0,
            "iso",
            0.5,
            0.5,
            100.0,
        )
        self.lmp.cmd.fix(4, "bond_react_MASTER_group", "temp/rescale", 50, T, T, 10, 1)
        self.lmp.cmd.run(1000)
        self.lmp.cmd.unfix(1)
        self.lmp.cmd.unfix(4)

        self.lmp.cmd.reset_atoms("id")

    def _polymerize_loop(self, fix_name: str, target: float, label: str) -> float:
        """Generic bond-formation loop. Returns the total bonds made."""
        bonds_made = 0
        for iteration in range(self.max_cycles):
            print("==================================")
            print(
                f"[{label}] iteration {iteration}  |  bonds = {bonds_made}  |  target = {target}"
            )
            print("==================================")

            self._run_nvt_npt_cycle()
            bonds_made = self.lmp.extract_fix(
                fix_name, LMP_STYLE_GLOBAL, LMP_TYPE_VECTOR
            )
            print(f"Bonds made = {bonds_made}")

            if bonds_made >= target:
                print("Target reached → stopping")
                break

        return bonds_made

    def _dump(self, dump_id, path: str, freq: int = 1000):
        self.lmp.cmd.dump(
            dump_id,
            "all",
            "custom",
            freq,
            path,
            "id",
            "type",
            "x",
            "y",
            "z",
            "vx",
            "vy",
            "vz",
        )
