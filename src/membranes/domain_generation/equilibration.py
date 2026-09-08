"""equilibration.py

Equilibrates a polymerized PA membrane to experimental density.
"""

import os
from membranes.domain_generation.initialize import initialize_simulation
from membranes.utils import random_int

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
WATER_BOND = 17
WATER_ANGLE = 22
WATER_ATOMS = (15, 16)
HYDROGEN_ATOMS_RV = (5, 8, 14)


class EquilibrationSimulation:
    """
    Post-polymerization equilibration of a PA membrane.

    Cycles through NVT (1000 K) → NVT (300 K) → NPT (300 K, ramped pressure)
    until the system reaches `rho_target` g/cm³, then runs a final 500 000-step
    NPT at 1 bar to confirm stability.
    """

    def __init__(
        self,
        in_dir: str,
        input_data: str,
        rho_target: float = 1.24,
        temperature: float = 300,
        nsteps: int = 50_000,
        seed: int | None = None,
    ):
        self.in_dir = in_dir
        self.input_data = input_data
        self.rho_target = rho_target
        self.temperature = temperature
        self.nsteps = nsteps
        self.seed = seed

        self.out_dir = f"data_out/{in_dir}/equilibration"
        os.makedirs(self.out_dir, exist_ok=True)

        self.wrapper = initialize_simulation(
            self.in_dir, log_file="logs/log.equilibration"
        )
        self.lmp = self.wrapper.lmp

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _density(self) -> float:
        self.lmp.cmd.variable("rho", "equal", "density")
        return self.lmp.extract_variable("rho")

    def _minimize(self):
        self.lmp.cmd.min_style("cg")
        self.lmp.cmd.min_modify("line", "quadratic")
        self.wrapper.minimize(
            etol=1.0e-4,
            ftol=1.0e-4,
            maxiter=10000,
            maxeval=100000,
        )

    def _apply_shake(self):
        self.wrapper.fix_shake(
            fix_id="FXSHAKE",
            group_id="all",
            tol=0.0001,
            iterations=10,
            n_stats=0,
            atom_type=HYDROGEN_ATOMS_RV,
        )

    def _nvt(self, temp: int, nsteps: int):
        self.wrapper.fix_nvt(
            fix_id=1,
            group_id="all",
            temperature_start=temp,
            temperature_end=temp,
        )
        self.lmp.cmd.velocity("all", "create", temp, random_int(self.seed))
        self.lmp.cmd.run(nsteps)
        self.lmp.cmd.unfix(1)

    def _npt(self, p_start: float, p_end: float, nsteps: int):
        self.lmp.cmd.unfix("FXSHAKE")
        self._minimize()
        self._apply_shake()
        self.wrapper.fix_npt(
            fix_id=1,
            group_id="all",
            temperature_start=self.temperature,
            temperature_end=self.temperature,
            pressure_start=p_start,
            pressure_end=p_end,
        )
        self.lmp.cmd.velocity("all", "scale", self.temperature)
        self.lmp.cmd.run(nsteps)
        self.lmp.cmd.unfix(1)

    # ------------------------------------------------------------------
    # Stages
    # ------------------------------------------------------------------

    def load(self):
        """Read structure, strip excess atoms"""
        self.lmp.cmd.read_data(self.input_data)

        self.lmp.cmd.dielectric(1.0)
        self.lmp.cmd.neighbor(2.0, "bin")
        self.lmp.cmd.neigh_modify("delay", 0, "every", 1, "check", "yes")
        self.lmp.cmd.timestep(1.0)
        self.lmp.cmd.run_style("verlet")
        self.lmp.cmd.thermo(5000)

    def equilibrate(self):
        """Cycle NVT/NPT stages until experimental density is reached.

        Pressure advances one step per cycle, mirroring the original script's
        `next Pvar` iterator: each NPT run goes from PRESSURE_RAMP[i] to
        PRESSURE_RAMP[i+1].  When the schedule is exhausted the index wraps
        back to the start so the ramp repeats if needed.
        """
        self._minimize()
        self._apply_shake()

        pressure_idx = 0
        n_up, n_down, n_max = 1, 0, 0
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
            self._nvt(temp=1000, nsteps=self.nsteps)

            # NVT at T (2x steps, matching original script)
            self._nvt(temp=self.temperature, nsteps=self.nsteps * 2)

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
                self.wrapper.fix_npt(
                    fix_id=1,
                    group_id="all",
                    temperature_start=self.temperature,
                    temperature_end=self.temperature,
                    pressure_start=1,
                    pressure_end=1,
                )

                self.lmp.cmd.velocity("all", "scale", self.temperature)
                self.lmp.cmd.run(5000)
                self.lmp.cmd.unfix(1)

                rho = self._density()
                if rho >= self.rho_target:
                    return  # Stable — proceed to finalize
                n_down += 1  # Density dropped — keep looping

    def finalize(self):
        """Long 1-bar NPT run and write final structure."""
        print("FINAL EQUILIBRATION RUN")
        self.wrapper.fix_npt(
            fix_id=1,
            group_id="all",
            temperature_start=self.temperature,
            temperature_end=self.temperature,
            pressure_start=1,
            pressure_end=1,
        )
        self.lmp.cmd.velocity("all", "scale", self.temperature)
        self.lmp.cmd.run(self.nsteps)
        self.lmp.cmd.unfix(1)
        rho = self._density()
        print(f"Final density = {rho:.4f} g/cm³")
        self.lmp.cmd.write_data(f"{self.out_dir}/equilibrated_polymer.lmps")
