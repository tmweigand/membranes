"""hydration.py

Hydrates an equilibrated PA membrane by adding water to feed/permeate
regions and running a pressurized filtration simulation.
"""

from membranes.domain_generation.lammps_init import LAMMPSInitialize
from membranes.domain_generation.style_properties import StyleProperties
from membranes.domain_generation.system_properties import SystemProperties


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

    def __init__(
        self,
        in_dir: str = "rv",
        mult: float = 0.1,
        rand: int = 1,
        input_data: str = "logs/equil_polymer.lmps",
        feed_pressure_atm: float = 0.5,
        perm_pressure_atm: float = 0.5,
        hydration_steps: int = 3_000_000,
        production_steps: int = 25_000_000,
    ):
        self.in_dir = in_dir
        self.mult = mult
        self.rand = rand
        self.input_data = input_data
        self.feed_pressure_atm = feed_pressure_atm
        self.perm_pressure_atm = perm_pressure_atm
        self.hydration_steps = hydration_steps
        self.production_steps = production_steps

        self.mult1d = mult ** (1 / 3)
        self.mult2d = mult ** (2 / 3)
        self.z_delta = 200 * self.mult1d
        self.num_h2o = int(3600 * self.mult2d)

        self._init_lammps()

    def _data_path(self, filename: str) -> str:
        return f"data_in/{self.in_dir}/{filename}"

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------

    def _init_lammps(self):
        wrapper = LAMMPSInitialize(
            log_file="logs/log.hydration", units="real", atom_style="full"
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

        SystemProperties(
            self.lmp, dimension=3, n_atom_types=30, boundary=("p", "p", "p")
        )

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
        self.lmp.cmd.minimize(1.0e-8, 1.0e-8, 1000, 100000)
        self.lmp.cmd.timestep(post_timestep)

    def _nvt(self, nsteps: int, temp: float = 300, damping: float = 100):
        self.lmp.cmd.fix(1, "all", "nvt", "temp", temp, temp, damping)
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

        self.lmp.cmd.velocity("all", "create", 300, 123 * self.rand)

        self._minimize(pre_timestep=10.0, post_timestep=2.0)
        self._nvt(5000)
        self._minimize(pre_timestep=10.0, post_timestep=2.0)
        self._nvt(5000)

    # ------------------------------------------------------------------
    # Stage 2 — unwrap membrane into a slab
    # ------------------------------------------------------------------

    def unwrap(self):
        """Cut the periodic membrane into a finite slab and extend the box in z."""
        m1 = self.mult1d
        z_delta = self.z_delta
        b = self._box()

        # Trim atoms outside the membrane slab
        self.lmp.cmd.region(
            "membranebox",
            "block",
            b["xlo"],
            b["xhi"],
            b["ylo"],
            b["yhi"],
            b["zlo"],
            b["zhi"] - 2.0 * m1,
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
            -z_delta,
            z_delta,
            "units",
            "box",
        )

        # Re-read bounds after change_box
        b = self._box()

        # Temporary LJ walls to contain atoms during minimization of cut edges
        self.lmp.cmd.fix(
            "zwalls",
            "all",
            "wall/lj126",
            "zlo",
            b["zlo"] + z_delta - 9 * m1,
            0.06844,
            3.40700,
            9,
            "zhi",
            b["zhi"] - z_delta + 9 * m1,
            0.06844,
            3.40700,
            9,
        )
        self.lmp.cmd.fix_modify("zwalls", "energy", "yes")

        self._minimize(pre_timestep=10.0, post_timestep=1.0)
        self._nvt(100)
        self._minimize(pre_timestep=10.0, post_timestep=2.0)
        self.lmp.cmd.unfix("zwalls")

        # Remove floating fragments outside the water compartment bounds
        b = self._box()
        self.lmp.cmd.region(
            "membranebox",
            "block",
            b["xlo"],
            b["xhi"],
            b["ylo"],
            b["yhi"],
            b["zlo"] + z_delta,
            b["zhi"] - z_delta,
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
        m1 = self.mult1d
        z_delta = self.z_delta
        b = self._box()

        self.lmp.cmd.region(
            "FEED",
            "block",
            "INF",
            "INF",
            "INF",
            "INF",
            b["zlo"] + 10 * m1,
            b["zlo"] + z_delta - 10 * m1,
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
            b["zhi"] - z_delta + 10 * m1,
            b["zhi"] - 10 * m1,
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
            521 + self.rand,
            "mol",
            "mol1",
            322 + self.rand,
        )
        self.lmp.cmd.create_atoms(
            0,
            "region",
            "PERM",
            "subset",
            self.num_h2o,
            512 + self.rand,
            "mol",
            "mol1",
            632 + self.rand,
        )

        self.lmp.cmd.group("FEEDWATER", "region", "FEED")
        self.lmp.cmd.group("PERMWATER", "region", "PERM")
        self.lmp.cmd.group("WATER", "type", *self.WATER_TYPES)

        self._minimize(pre_timestep=10.0, post_timestep=2.0)

    # ------------------------------------------------------------------
    # Stage 4 — add pistons and backing layer
    # ------------------------------------------------------------------

    def add_pistons(self):
        """Add graphene-like piston planes and a polysulfone backing layer."""
        m1 = self.mult1d
        z_delta = self.z_delta
        glc = self.GRAPHENE_LC

        self.lmp.cmd.lattice("hcp", glc)
        b = self._box()

        # --- Backing layer (polysulfone pin) on the permeate side ---
        pin_lo = b["zhi"] - z_delta + 5 * m1
        pin_hi = pin_lo + glc
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
            b["zlo"],
            b["zlo"] + glc * 0.75,
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
            b["zhi"] - glc * 0.75,
            b["zhi"],
            "units",
            "box",
        )
        self.lmp.cmd.create_atoms(self.PISTON_TYPE, "region", "zwallHI")
        self._make_group("HIzwall", "zwallHI", self.PISTON_TYPE)

        self.lmp.cmd.group("zwalls", "union", "LOzwall", "HIzwall")
        self.lmp.cmd.group("mobile", "subtract", "all", "zwalls")

        # Expand box and add hard bounding LJ walls for the pistons
        self.lmp.cmd.change_box(
            "all", "z", "delta", -200 * m1, 200 * m1, "units", "box"
        )
        b = self._box()
        self.lmp.cmd.fix(
            "zwalls1",
            "all",
            "wall/lj126",
            "zlo",
            b["zlo"],
            0.06844,
            3.40700,
            9,
            "zhi",
            b["zhi"],
            0.06844,
            3.40700,
            9,
        )
        self.lmp.cmd.fix_modify("zwalls1", "energy", "yes")

        # --- Pressure forces on pistons ---
        # Compute piston area once — box xy dimensions are fixed from here on
        b = self._box()
        piston_area = (b["xhi"] - b["xlo"]) * (b["yhi"] - b["ylo"])

        a2m = self.ANG_TO_M
        f2n = self.FINMD_TO_FINN
        feed_p = self.feed_pressure_atm * self.ATM_TO_PA
        perm_p = self.perm_pressure_atm * self.ATM_TO_PA

        pos_force = piston_area * a2m**2 * feed_p / f2n
        neg_force = -piston_area * a2m**2 * perm_p / f2n

        # Zero x/y forces and set z forces on each piston plane so they move
        # only axially. Using setforce + addforce avoids fix/rigid, which
        # conflicts with fix/shake when atom groups overlap across procs.
        self.lmp.cmd.fix("FEEDFORCE", "LOzwall", "setforce", 0.0, 0.0, "NULL")
        self.lmp.cmd.fix("PERMFORCE", "HIzwall", "setforce", 0.0, 0.0, "NULL")
        self.lmp.cmd.fix("FEEDPRESS", "LOzwall", "aveforce", 0.0, 0.0, pos_force)
        self.lmp.cmd.fix("PERMPRESS", "HIzwall", "aveforce", 0.0, 0.0, neg_force)

    # ------------------------------------------------------------------
    # Stage 5 — run hydration experiment
    # ------------------------------------------------------------------

    def run_hydration(self):
        """Set up computes, dumps, SHAKE, and run the full filtration MD."""
        self._setup_computes()
        self._setup_dumps()

        self.lmp.cmd.reset_atoms("id")
        self.lmp.cmd.write_data("test.lmps")
        self.lmp.cmd.restart(5000, "restarts_hydr/restart.*")

        # Exclude piston self-interactions
        self.lmp.cmd.neigh_modify("exclude", "type", self.PISTON_TYPE, self.PISTON_TYPE)

        # SHAKE on WATER only — piston atoms (type 17) have no bonds
        self.lmp.cmd.fix(
            "FXSHAKE",
            "WATER",
            "shake",
            0.0001,
            20,
            0,
            "b",
            17,
            "a",
            22,
            "t",
            5,
            8,
            14,
        )

        # NVT on mobile atoms only — pistons are driven by aveforce
        self.lmp.cmd.fix(1, "mobile", "nvt", "temp", 300, 300, 200)
        self.lmp.cmd.thermo(100)
        self.lmp.cmd.thermo_style(
            "custom", "step", "temp", "press", "etotal", "ke", "pe"
        )

        # Ramp timestep from 0.5 → 1.0 → 2.0 fs.
        # Water molecules placed near pistons can have large initial forces;
        # a sudden 1 fs step causes H atoms to fly off (bond explosion) which
        # is what SHAKE reports as "missing atoms".
        self.lmp.cmd.timestep(0.5)

        # Dump atom types at step 0 to help diagnose any future SHAKE errors
        self.lmp.cmd.dump(
            "debugdump",
            "all",
            "custom",
            1,
            "logs/pre_run_atoms.txt",
            "id",
            "type",
            "mol",
            "x",
            "y",
            "z",
        )
        self.lmp.cmd.run(0)
        self.lmp.cmd.undump("debugdump")

        self.lmp.cmd.run(200)
        self.lmp.cmd.timestep(1.0)
        self.lmp.cmd.run(500)

        self.lmp.cmd.write_restart("restarts_hydr/hydr_pre.restart")
        self.lmp.cmd.run(self.hydration_steps)
        self.lmp.cmd.write_restart("restarts_hydr/hydr_init.restart")

        self.lmp.cmd.timestep(2.0)
        self.lmp.cmd.run(self.production_steps)

        self.lmp.cmd.unfix(1)
        self.lmp.cmd.unfix("FXSHAKE")
        self.lmp.cmd.write_data("hydrated_data.lmps")

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
            "logs/pressuredata.*.gz",
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
            "logs/membranedata.*.gz",
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
            "logs/systemdata.*.gz",
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

    # ------------------------------------------------------------------
    # Top-level runner
    # ------------------------------------------------------------------

    def run(self):
        self.load()
        self.unwrap()
        self.add_water()
        self.add_pistons()
        self.run_hydration()
