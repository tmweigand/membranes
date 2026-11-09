"""lammps_init.py"""

from lammps import lammps


class LAMMPSInitialize:
    """Base class providing LAMMPS instance to all child classes"""

    def __init__(
        self,
        log_file: str = "lammps.log",
        units: str = "real",
        atom_style: str = "full",
        newton: str = "off",
        echo: str = "none",
    ):

        self.log_file = log_file
        self.units = units
        self.atom_style = atom_style
        self.newton = newton
        self.echo = echo

        self.lmp = lammps(
            cmdargs=[
                "-screen",
                "none",
                "-nocite",
                "-log",
                log_file,
            ]
        )
        self.lmp.cmd.units(units)
        self.lmp.cmd.atom_style(atom_style)
        self.lmp.cmd.newton(newton)
        self.lmp.cmd.echo(echo)

    def __del__(self):
        """Cleanup LAMMPS instance"""
        if hasattr(self, "lammps"):
            self.lmp.close()

    # ── Readable wrappers for common LAMMPS commands ───────────────────

    def minimize(
        self,
        etol: float,
        ftol: float,
        maxiter: int,
        maxeval: int,
    ):
        """Wrapper to allow keyword argument to lammps minimize

        Parameters:
            etol: stopping tolerance for energy (unitless)
            ftol: stopping tolerance for force (force units)
            maxiter: max iterations of minimizer
            maxeval: max number of force/energy evaluations
        """
        self.lmp.cmd.minimize(etol, ftol, maxiter, maxeval)

    def fix_nvt(
        self,
        fix_id: int | str,
        group_id: str,
        temperature_start: float,
        temperature_end: float,
        temperature_damp: float = 100.0,
    ):
        """Wrapper to allow keyword argument to lammps fix nvt

        Parameters:
            fix_id: label for fix
            group_id: name of fix
            temperature_start: external temperature at start of run
            temperature_end: external temperature at end of run
            temperature_damp: temperature damping parameter (time units)
        """
        self.lmp.cmd.fix(
            fix_id,
            group_id,
            "nvt",
            "temp",
            temperature_start,
            temperature_end,
            temperature_damp,
        )

    def fix_npt(
        self,
        fix_id: int | str,
        group_id: str,
        temperature_start: float,
        temperature_end: float,
        pressure_start: float,
        pressure_end: float,
        temperature_damp: float = 100.0,
        pressure_damp: float = 100.0,
    ):
        """Wrapper to allow keyword argument to lammps fix npt

        Parameters:
            fix_id: label for fix
            group_id: name of fix
            temperature_start: external temperature at start of run
            temperature_end: external temperature at end of run
            temperature_damp: temperature damping parameter (time units)
            pressure_start: external pressure at start of run
            pressure_end: external pressure at end of run
            pressure_damp: pressure damping parameter (time units)
        """
        self.lmp.cmd.fix(
            fix_id,
            group_id,
            "npt",
            "temp",
            temperature_start,
            temperature_end,
            temperature_damp,
            "iso",
            pressure_start,
            pressure_end,
            pressure_damp,
        )

    def fix_rescale_temp(
        self,
        fix_id: int | str,
        group_id: str,
        n_steps: int,
        temperature_start: float,
        temperature_end: float,
        window: float,
        fraction: float,
    ):
        """Wrapper to allow keyword argument to lammps fix temp/rescale

        Rescaling is only performed if the difference between the current and desired temperatures is greater than the window value. The amount of rescaling that is applied is a fraction (from 0.0 to 1.0) of the difference between the actual and desired temperature. E.g. if fraction = 1.0, the temperature is reset to exactly the desired value.

        Parameters:
            fix_id: label for fix
            group_id: name of fix
            n_steps: Perform rescaling every N steps
            temperature_start: external temperature at start of run
            temperature_end: external temperature at end of run
            window: only rescale if temperature is outside this window (temperature units)
            fraction: rescale to target temperature by this fraction
        """
        self.lmp.cmd.fix(
            fix_id,
            group_id,
            "temp/rescale",
            n_steps,
            temperature_start,
            temperature_end,
            window,
            fraction,
        )

    def fix_shake(
        self,
        fix_id: int | str,
        group_id: str,
        tol: float,
        iterations: int,
        n_stats: int,
        bond_type: int | list[int] | None = None,
        angle_type: int | list[int] | None = None,
        atom_type: int | list[int] | None = None,
    ):
        """Wrapper to allow keyword argument to lammps fix shake

        Apply bond and angle constraints to specified bonds and angles in the simulation by either the SHAKE algorithm. This typically enables a longer timestep. The SHAKE constraint algorithm, however, can only be applied during molecular dynamics runs.

        Parameters:
            fix_id: label for fix
            group_id: name of fix
            tol: accuracy tolerance of SHAKE solution
            iterations: max iterations in each SHAKE solution
            n_stats: print SHAKE statistics every n timesteps
            bond_type:
            angle_type:
            atom_type:
        """
        args = [
            fix_id,
            group_id,
            "shake",
            tol,
            iterations,
            n_stats,
        ]

        if bond_type is not None:
            bond_type = _expand(bond_type)
            args += ["b", *bond_type]

        if angle_type is not None:
            angle_type = _expand(angle_type)
            args += ["a", *angle_type]

        if atom_type is not None:
            atom_type = _expand(atom_type)
            args += ["t", *atom_type]

        self.lmp.cmd.fix(*args)


def _expand(arg):
    if isinstance(arg, (list, tuple)):
        return list(arg)
    return [arg]
