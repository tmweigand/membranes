"""hydration.py

Hydrates an equilibrated PA membrane by adding water to feed/permeate
regions and running a pressurized filtration simulation.
"""

import os
from membranes.domain_generation.initialize import initialize_simulation
from membranes.utils import random_int


class HydrationSimulation:
    """
    Hydrates a PA membrane and runs a pressurized filtration experiment.

    Stages
    ------
    1. load          — read equilibrated structure, shift box to origin, minimize
    2. unwrap        — cut periodic membrane into a slab, extend box in z
    3. add_water     — place TIP3P water in feed and permeate regions
    4. add_pistons   — add graphene piston planes and polysulfone backing layer
    5. run_hydration — apply pressure forces and run filtration MD
    """

    # Unit conversion constants
    FINMD_TO_FINN = 6.94768717535204e-11
    ANG_TO_M = 1e-10
    ATM_TO_PA = 101325

    # Graphene lattice constant (Å)
    GRAPHENE_LC = 2.46

    # Atom types
    WATER_TYPES = (15, 16)
    PISTON_TYPE = 17
    WATER_BOND = 17
    WATER_ANGLE = 22
    WATER_ATOMS = (15, 16)
    HYDROGEN_ATOMS_RV = (5, 8, 14)

    def __init__(
        self,
        in_dir: str,
        input_data: str,
        multiple: float = 0.1,
        temperature: float = 300,
        feed_pressure_atm: float = 0.5,
        perm_pressure_atm: float = 0.5,
        hydration_steps: int = 3_000_000,
        production_steps: int = 25_000_000,
        seed: int | None = None,
    ):
        self.in_dir = in_dir
        self.multiple = multiple
        self.input_data = input_data
        self.temperature = temperature
        self.feed_pressure_atm = feed_pressure_atm
        self.perm_pressure_atm = perm_pressure_atm
        self.hydration_steps = hydration_steps
        self.production_steps = production_steps
        self.seed = seed

        self.mult1d = multiple ** (1 / 3)
        self.mult2d = multiple ** (2 / 3)
        self.z_delta = 200 * self.mult1d
        self.num_h2o = int(3600 * self.mult2d)

        self.out_dir = f"data_out/{in_dir}/hydration"
        os.makedirs(self.out_dir, exist_ok=True)
        os.makedirs(f"{self.out_dir}/restarts", exist_ok=True)

        self.wrapper = initialize_simulation(self.in_dir, log_file="logs/log.hydration")
        self.lmp = self.wrapper.lmp

    def _data_path(self, filename: str) -> str:
        return f"data_in/{self.in_dir}/{filename}"

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _box(self) -> dict:
        """Return current box bounds as plain Python floats."""
        return {
            "xlo": self.lmp.get_thermo("xlo"),
            "xhi": self.lmp.get_thermo("xhi"),
            "ylo": self.lmp.get_thermo("ylo"),
            "yhi": self.lmp.get_thermo("yhi"),
            "zlo": self.lmp.get_thermo("zlo"),
            "zhi": self.lmp.get_thermo("zhi"),
        }

    def _minimize(self, pre_timestep: float = 10.0, post_timestep: float = 2.0):
        self.lmp.cmd.timestep(pre_timestep)
        self.wrapper.minimize(
            etol=1.0e-8,
            ftol=1.0e-8,
            maxiter=1000,
            maxeval=100000,
        )
        self.lmp.cmd.timestep(post_timestep)

    def _nvt(self, nsteps: int, temp: float):
        self.wrapper.fix_nvt(
            fix_id=1, group_id="all", temperature_start=temp, temperature_end=temp
        )
        self.lmp.cmd.run(nsteps)
        self.lmp.cmd.unfix(1)

    def _make_group(self, name: str, region: str, atom_type: int):
        """Intersect a region and atom type to cleanly define a plane group."""
        self.lmp.cmd.group(f"{name}1", "type", atom_type)
        self.lmp.cmd.group(f"{name}2", "region", region)
        self.lmp.cmd.group(name, "intersect", f"{name}1", f"{name}2")
        self.lmp.cmd.group(f"{name}1", "delete")
        self.lmp.cmd.group(f"{name}2", "delete")

    # ------------------------------------------------------------------
    # Stage 1 — load and prepare
    # ------------------------------------------------------------------

    def load(self):
        """Read equilibrated structure, shift box origin to zero, and minimize."""
        self.lmp.cmd.read_data(self.input_data)

        self.lmp.cmd.dielectric(1.0)
        self.lmp.cmd.neighbor(4.0, "bin")
        self.lmp.cmd.neigh_modify("delay", 0, "every", 1, "check", "yes")
        self.lmp.cmd.run_style("verlet")
        self.lmp.cmd.lattice("fcc", 4.3)

        # Shift box so all coordinates start at 0. commands_string is required
        # here so LAMMPS evaluates the $(xlo) inline expressions itself.
        self.lmp.commands_string(
            "change_box all x delta -$(xlo) -$(xlo) y delta -$(ylo) -$(ylo) "
            "z delta -$(zlo) -$(zlo) units box\n"
            "displace_atoms all move -$(xlo) -$(ylo) -$(zlo) units box"
        )

        self.lmp.cmd.velocity("all", "create", 300, random_int(self.seed))

        self._minimize(pre_timestep=10.0, post_timestep=2.0)
        self._nvt(5000, self.temperature)
        self._minimize(pre_timestep=10.0, post_timestep=2.0)
        self._nvt(5000, self.temperature)

    # ------------------------------------------------------------------
    # Stage 2 — unwrap membrane into a slab
    # ------------------------------------------------------------------

    def unwrap(self):
        """Cut the periodic membrane into a finite slab and extend the box in z."""
        box = self._box()

        # Trim atoms outside the membrane slab
        self.lmp.cmd.region(
            "membranebox",
            "block",
            box["xlo"],
            box["xhi"],
            box["ylo"],
            box["yhi"],
            box["zlo"],
            box["zhi"] - 2.0 * self.mult1d,
            "side",
            "out",
            "units",
            "box",
        )
        self.lmp.cmd.delete_atoms("region", "membranebox", "bond", "yes")
        self.lmp.cmd.reset_atoms("id", "sort", "yes")
        self.lmp.cmd.region("membranebox", "delete")

        # Tag remaining atoms as the membrane group
        self.lmp.cmd.region(
            "membranebox", "block", "INF", "INF", "INF", "INF", "INF", "INF"
        )
        self.lmp.cmd.group("membrane", "region", "membranebox")
        self.lmp.cmd.region("membranebox", "delete")

        # Switch to non-periodic z and extend the box for water compartments
        self.lmp.cmd.kspace_modify("slab", 3.0)
        self.lmp.cmd.change_box(
            "all",
            "boundary",
            "p",
            "p",
            "f",
            "z",
            "delta",
            -self.z_delta,
            self.z_delta,
            "units",
            "box",
        )

        # Re-read bounds after change_box
        box = self._box()

        # Temporary LJ walls to contain atoms during minimization of cut edges
        self.lmp.cmd.fix(
            "zwalls",
            "all",
            "wall/lj126",
            "zlo",
            box["zlo"] + self.z_delta - 9 * self.mult1d,
            0.06844,
            3.40700,
            9,
            "zhi",
            box["zhi"] - self.z_delta + 9 * self.mult1d,
            0.06844,
            3.40700,
            9,
        )
        self.lmp.cmd.fix_modify("zwalls", "energy", "yes")

        self._minimize(pre_timestep=10.0, post_timestep=1.0)
        self._nvt(100, self.temperature)
        self._minimize(pre_timestep=10.0, post_timestep=2.0)
        self.lmp.cmd.unfix("zwalls")

        # Remove floating fragments outside the water compartment bounds
        box = self._box()
        self.lmp.cmd.region(
            "membranebox",
            "block",
            box["xlo"],
            box["xhi"],
            box["ylo"],
            box["yhi"],
            box["zlo"] + self.z_delta,
            box["zhi"] - self.z_delta,
            "side",
            "out",
            "units",
            "box",
        )
        self.lmp.cmd.delete_atoms("region", "membranebox", "bond", "yes")
        self.lmp.cmd.region("membranebox", "delete")
        self.lmp.cmd.reset_atoms("id")

    # ------------------------------------------------------------------
    # Stage 3 — add water
    # ------------------------------------------------------------------

    def add_water(self):
        """Place TIP3P water molecules in the feed and permeate compartments."""
        b = self._box()

        self.lmp.cmd.region(
            "FEED",
            "block",
            "INF",
            "INF",
            "INF",
            "INF",
            b["zlo"] + 10 * self.mult1d,
            b["zlo"] + self.z_delta - 10 * self.mult1d,
            "units",
            "box",
        )
        self.lmp.cmd.region(
            "PERM",
            "block",
            "INF",
            "INF",
            "INF",
            "INF",
            b["zhi"] - self.z_delta + 10 * self.mult1d,
            b["zhi"] - 10 * self.mult1d,
            "units",
            "box",
        )

        self.lmp.cmd.molecule("mol1", self._data_path("tip3p.mol"))
        self.lmp.cmd.create_atoms(
            0,
            "region",
            "FEED",
            "subset",
            self.num_h2o,
            random_int(self.seed),
            "mol",
            "mol1",
            random_int(self.seed),
        )
        self.lmp.cmd.create_atoms(
            0,
            "region",
            "PERM",
            "subset",
            self.num_h2o,
            random_int(self.seed),
            "mol",
            "mol1",
            random_int(self.seed),
        )

        self.lmp.cmd.group("FEEDWATER", "region", "FEED")
        self.lmp.cmd.group("PERMWATER", "region", "PERM")
        self.lmp.cmd.group("WATER", "type", *self.WATER_TYPES)

        self._minimize(pre_timestep=10.0, post_timestep=2.0)

        self.lmp.cmd.write_data(f"{self.out_dir}/water_added.lmps")

    # ------------------------------------------------------------------
    # Stage 4 — add pistons and backing layer
    # ------------------------------------------------------------------

    def add_pistons(self):
        """Add graphene-like piston planes and a polysulfone backing layer."""

        self.lmp.cmd.lattice("hcp", self.GRAPHENE_LC)
        box = self._box()

        # --- Backing layer (polysulfone pin) on the permeate side ---
        pin_lo = box["zhi"] - self.z_delta + 5 * self.mult1d
        pin_hi = pin_lo + self.GRAPHENE_LC
        self.lmp.cmd.region(
            "PINlayer",
            "block",
            "INF",
            "INF",
            "INF",
            "INF",
            pin_lo,
            pin_hi,
            "units",
            "box",
        )
        self.lmp.cmd.create_atoms(self.PISTON_TYPE, "region", "PINlayer")
        self._make_group("PINplane", "PINlayer", self.PISTON_TYPE)
        self.lmp.cmd.fix("PINFORCE", "PINplane", "setforce", 0.0, 0.0, 0.0)
        self.lmp.cmd.neigh_modify("exclude", "group", "WATER", "PINplane")

        # --- Feed-side (low-z) piston ---
        self.lmp.cmd.region(
            "zwallLO",
            "block",
            "INF",
            "INF",
            "INF",
            "INF",
            box["zlo"],
            box["zlo"] + self.GRAPHENE_LC * 0.75,
            "units",
            "box",
        )
        self.lmp.cmd.create_atoms(self.PISTON_TYPE, "region", "zwallLO")
        self._make_group("LOzwall", "zwallLO", self.PISTON_TYPE)

        # --- Permeate-side (high-z) piston ---
        self.lmp.cmd.region(
            "zwallHI",
            "block",
            "INF",
            "INF",
            "INF",
            "INF",
            box["zhi"] - self.GRAPHENE_LC * 0.75,
            box["zhi"],
            "units",
            "box",
        )
        self.lmp.cmd.create_atoms(self.PISTON_TYPE, "region", "zwallHI")
        self._make_group("HIzwall", "zwallHI", self.PISTON_TYPE)

        self.lmp.cmd.group("zwalls", "union", "LOzwall", "HIzwall")
        self.lmp.cmd.group("mobile", "subtract", "all", "zwalls")

        # Expand box and add hard bounding LJ walls for the pistons
        self.lmp.cmd.change_box(
            "all", "z", "delta", -200 * self.mult1d, 200 * self.mult1d, "units", "box"
        )
        box = self._box()
        self.lmp.cmd.fix(
            "zwalls1",
            "all",
            "wall/lj126",
            "zlo",
            box["zlo"],
            0.06844,
            3.40700,
            9,
            "zhi",
            box["zhi"],
            0.06844,
            3.40700,
            9,
        )
        self.lmp.cmd.fix_modify("zwalls1", "energy", "yes")

        # --- Pressure forces on pistons ---
        # Compute piston area once — box xy dimensions are fixed from here on
        box = self._box()
        piston_area = (box["xhi"] - box["xlo"]) * (box["yhi"] - box["ylo"])

        feed_p = self.feed_pressure_atm * self.ATM_TO_PA
        perm_p = self.perm_pressure_atm * self.ATM_TO_PA

        pos_force = piston_area * self.ANG_TO_M**2 * feed_p / self.FINMD_TO_FINN
        neg_force = -piston_area * self.ANG_TO_M**2 * perm_p / self.FINMD_TO_FINN

        self.lmp.cmd.fix("FEEDFORCE", "LOzwall", "addforce", 0.0, 0.0, pos_force)
        self.lmp.cmd.fix("PERMFORCE", "HIzwall", "addforce", 0.0, 0.0, neg_force)

        self.lmp.cmd.fix(
            "LOPISTON",
            "LOzwall",
            "rigid",
            "single",
            "torque",
            "*",
            "off",
            "off",
            "off",
            "force",
            "*",
            "off",
            "off",
            "on",
        )
        self.lmp.cmd.fix(
            "HIPISTON",
            "HIzwall",
            "rigid",
            "single",
            "torque",
            "*",
            "off",
            "off",
            "off",
            "force",
            "*",
            "off",
            "off",
            "on",
        )

        self.lmp.cmd.write_data(f"{self.out_dir}/pistons_added.lmps")

    # ------------------------------------------------------------------
    # Stage 5 — run hydration experiment
    # ------------------------------------------------------------------

    def run_hydration(self):
        """Set up computes, dumps, SHAKE, and run the full filtration MD."""
        self._setup_computes()
        self._setup_dumps()

        self.lmp.cmd.reset_atoms("id")
        self.lmp.cmd.write_data("test.lmps")
        self.lmp.cmd.restart(5000, f"{self.out_dir}/restarts/hydration_restart.*")

        # Exclude piston self-interactions
        self.lmp.cmd.neigh_modify("exclude", "type", self.PISTON_TYPE, self.PISTON_TYPE)

        # SHAKE on WATER only — piston atoms (type 17) have no bonds
        self.wrapper.fix_shake(
            fix_id="FXSHAKE",
            group_id="all",
            tol=0.0001,
            iterations=20,
            n_stats=0,
            bond_type=self.WATER_BOND,
            angle_type=self.WATER_ANGLE,
            atom_type=self.HYDROGEN_ATOMS_RV,
        )

        # NVT on mobile atoms only — pistons are driven by aveforce
        self.wrapper.fix_nvt(
            fix_id=1,
            group_id="mobile",
            temperature_start=self.temperature,
            temperature_end=self.temperature,
            temperature_damp=200,
        )
        self.lmp.cmd.thermo(1000)
        self.lmp.cmd.thermo_style(
            "custom", "step", "temp", "press", "etotal", "ke", "pe"
        )

        self.lmp.cmd.timestep(1.0)

        # Dump atom types at step 0 to help diagnose any future SHAKE errors
        self.lmp.cmd.dump(
            "debugdump",
            "all",
            "custom",
            100,
            f"{self.out_dir}/pre_run_atoms.txt",
            "id",
            "type",
            "mol",
            "x",
            "y",
            "z",
        )

        self.lmp.cmd.write_restart(f"{self.out_dir}/restarts/pre_hydration.restart")
        self.lmp.cmd.run(self.hydration_steps)
        self.lmp.cmd.write_restart(f"{self.out_dir}/restarts/initial_hydration.restart")

        self.lmp.cmd.timestep(2.0)
        self.lmp.cmd.run(self.production_steps)

        self.lmp.cmd.unfix(1)
        self.lmp.cmd.unfix("FXSHAKE")
        self.lmp.cmd.write_data(f"{self.out_dir}/hydrated_data.lmps")

    def _setup_computes(self):
        """Register per-atom stress and Voronoi computes for water molecules."""
        for name, keywords in [
            ("peratomstressvol", ["NULL"]),
            ("peratomstressvolke", ["NULL", "ke"]),
            ("peratomstressvolkspace", ["NULL", "kspace"]),
            ("peratomstressvolfix", ["NULL", "fix"]),
        ]:
            self.lmp.cmd.compute(name, "WATER", "stress/atom", *keywords)
        self.lmp.cmd.compute("peratomvol", "WATER", "voronoi/atom")

    def _setup_dumps(self):
        """Register trajectory dumps for water, membrane, and full system."""
        stress_fields = [
            f"c_{name}[{i}]"
            for name in (
                "peratomstressvol",
                "peratomstressvolke",
                "peratomstressvolkspace",
                "peratomstressvolfix",
            )
            for i in (1, 2, 3)
        ]
        self.lmp.cmd.dump(
            "waterdata",
            "WATER",
            "custom/gz",
            100000,
            f"{self.out_dir}/pressuredata.*.gz",
            "id",
            "mol",
            "type",
            "mass",
            "q",
            "x",
            "y",
            "z",
            "vx",
            "vy",
            "vz",
            *stress_fields,
            "c_peratomvol[1]",
        )
        self.lmp.cmd.dump_modify("waterdata", "sort", 2)

        self.lmp.cmd.dump(
            "membranedata",
            "membrane",
            "custom/gz",
            100000,
            f"{self.out_dir}/membranedata.*.gz",
            "id",
            "mol",
            "type",
            "mass",
            "q",
            "x",
            "y",
            "z",
        )
        self.lmp.cmd.dump_modify("membranedata", "sort", 1)

        self.lmp.cmd.dump(
            "systemdata",
            "all",
            "custom/gz",
            100000,
            f"{self.out_dir}/systemdata.*.gz",
            "id",
            "mol",
            "type",
            "mass",
            "q",
            "x",
            "y",
            "z",
        )
        self.lmp.cmd.dump_modify("systemdata", "sort", 1)
