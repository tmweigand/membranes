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

        self.lmp = lammps(cmdargs=["-screen", "none", "-nocite"])
        self.lmp.cmd.log(log_file)
        self.lmp.cmd.units(units)
        self.lmp.cmd.atom_style(atom_style)
        self.lmp.cmd.newton(newton)
        self.lmp.cmd.echo(echo)

    def __del__(self):
        """Cleanup LAMMPS instance"""
        if hasattr(self, "lammps"):
            self.lmp.close()
