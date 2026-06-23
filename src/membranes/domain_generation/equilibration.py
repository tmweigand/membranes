"""equilibration.py

Equilibrates a polymerized PA membrane to experimental density.
"""

import numpy as np
from membranes.domain_generation.lammps_init import LAMMPSInitialize
from membranes.domain_generation.style_properties import StyleProperties
from membranes.domain_generation.system_properties import SystemProperties

# Pressure ramp schedule (atm)
PRESSURE_RAMP = [
    1,
    1,
    100,
    1000,
    2000,
    5000,
    10000,
    15000,
    20000,
    25000,
    30000,
    35000,
    40000,
    45000,
    50000,
]

# Atom types to remove before equilibration
EXCESS_TYPES = (11, 13, 18)


class EquilibrationSimulation:
    """
    Post-polymerization equilibration of a PA membrane.

    Cycles through NVT (1000 K) → NVT (300 K) → NPT (300 K, ramped pressure)
    until the system reaches `rho_target` g/cm³, then runs a final 500 000-step
    NPT at 1 bar to confirm stability.
    """

    def __init__(
        self,
        in_dir: str = "rv",
        input_data: str = "logs/term_final.lmps",
        rho_target: float = 1.24,
        temperature: float = 300,
        nsteps: int = 50_000,
        initial_seed: int = 58447419,
    ):
        self.in_dir = in_dir
        self.input_data = input_data
        self.rho_target = rho_target
        self.temperature = temperature
        self.nsteps = nsteps
        self.initial_seed = initial_seed

        self._init_lammps()

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------

    def _init_lammps(self):
        wrapper = LAMMPSInitialize(
            log_file="logs/log.EQ", units="real", atom_style="full"
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

    def _density(self) -> float:
        self.lmp.cmd.variable("rho", "equal", "density")
        return self.lmp.extract_variable("rho")

    def _minimize(self):
        self.lmp.cmd.min_style("cg")
        self.lmp.cmd.min_modify("line", "quadratic")
        self.lmp.cmd.minimize(1.0e-4, 1.0e-4, 10000, 100000)

    def _apply_shake(self):
        self.lmp.cmd.fix(
            "FXSHAKE", "all", "shake", 0.0001, 10, 0, "b", 17, "a", 22, "t", 5, 8, 14
        )

    def _nvt(self, temp: int, nsteps: int, seed: int | None = None):
        self.lmp.cmd.fix(1, "all", "nvt", "temp", temp, temp, 100)
        if seed is not None:
            self.lmp.cmd.velocity("all", "create", temp, seed)
        else:
            self.lmp.cmd.velocity("all", "scale", temp)
        self.lmp.cmd.run(nsteps)
        self.lmp.cmd.unfix(1)

    def _npt(self, p_start: float, p_end: float, nsteps: int):
        self.lmp.cmd.unfix("FXSHAKE")
        self._minimize()
        self._apply_shake()
        self.lmp.cmd.fix(
            1,
            "all",
            "npt",
            "temp",
            self.temperature,
            self.temperature,
            100,
            "iso",
            p_start,
            p_end,
            100,
        )
        self.lmp.cmd.velocity("all", "scale", self.temperature)
        self.lmp.cmd.run(nsteps)
        self.lmp.cmd.unfix(1)

    # ------------------------------------------------------------------
    # Stages
    # ------------------------------------------------------------------

    def load(self):
        """Read structure, strip excess atoms, minimize, and apply SHAKE."""
        self.lmp.cmd.read_data(self.input_data)

        self.lmp.cmd.dielectric(1.0)
        self.lmp.cmd.neighbor(2.0, "bin")
        self.lmp.cmd.neigh_modify("delay", 0, "every", 1, "check", "yes")
        self.lmp.cmd.timestep(1.0)
        self.lmp.cmd.run_style("verlet")
        self.lmp.cmd.thermo(5000)

        self._minimize()
        self._apply_shake()

    def equilibrate(self):
        """Cycle NVT/NPT stages until experimental density is reached.

        Pressure advances one step per cycle, mirroring the original script's
        `next Pvar` iterator: each NPT run goes from PRESSURE_RAMP[i] to
        PRESSURE_RAMP[i+1].  When the schedule is exhausted the index wraps
        back to the start so the ramp repeats if needed.
        """
        pressure_idx = 0
        n_up, n_down, n_max = 1, 0, 0
        first_run = True

        while True:
            p_start = PRESSURE_RAMP[pressure_idx]
            p_end = PRESSURE_RAMP[pressure_idx + 1]
            rho = self._density()

            print("==================================")
            print(f"  Ramp-up / down / max : {n_up} / {n_down} / {n_max}")
            print(f"  Density : {rho:.4f} g/cm3  (target >= {self.rho_target})")
            print(f"  Pressure: {p_start} -> {p_end} atm")
            print("==================================")

            # NVT at 1000 K
            self._nvt(1000, self.nsteps, seed=self.initial_seed if first_run else None)
            first_run = False

            # NVT at T (2x steps, matching original script)
            self._nvt(self.temperature, self.nsteps * 2)

            # NPT at current pressure step
            self._npt(p_start, p_end, self.nsteps)

            rho = self._density()

            # Advance pressure index; wrap and increment n_max when exhausted
            if pressure_idx < len(PRESSURE_RAMP) - 2:
                pressure_idx += 1
                n_up += 1
            else:
                n_max += 1
                pressure_idx = 0  # restart ramp

            # Once dense enough at meaningful pressure, test stability at 1 bar
            if rho >= self.rho_target and p_start > 2000:
                self.lmp.cmd.fix(
                    1,
                    "all",
                    "npt",
                    "temp",
                    self.temperature,
                    self.temperature,
                    100,
                    "iso",
                    1,
                    1,
                    100,
                )
                self.lmp.cmd.velocity("all", "scale", self.temperature)
                self.lmp.cmd.run(5000)
                self.lmp.cmd.unfix(1)
                self.lmp.cmd.write_data("init_equil_polym.lmps")

                rho = self._density()
                if rho >= self.rho_target:
                    return  # Stable — proceed to finalize
                n_down += 1  # Density dropped — keep looping

    def finalize(self):
        """Long 1-bar NPT run and write final structure."""
        print("FINAL EQUILIBRATION RUN")
        self.lmp.cmd.fix(
            1,
            "all",
            "npt",
            "temp",
            self.temperature,
            self.temperature,
            100,
            "iso",
            1,
            1,
            100,
        )
        self.lmp.cmd.velocity("all", "scale", self.temperature)
        self.lmp.cmd.run(self.nsteps)
        self.lmp.cmd.unfix(1)
        rho = self._density()
        print(f"Final density = {rho:.4f} g/cm³")
        self.lmp.cmd.write_data("logs/equil_polymer.lmps")

    def run(self):
        self.load()
        self.equilibrate()
        self.finalize()
