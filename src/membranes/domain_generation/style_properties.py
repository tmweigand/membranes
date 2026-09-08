"""style_properties.py"""


class StyleProperties:
    """Style properties for LAMMPS simulation"""

    def __init__(
        self,
        lmp: "lammps",
        pair_style: str = "lj/charmm/coul/long",
        pair_cutoff_inner: float = 7.0,
        pair_cutoff_outer: float = 9.0,
        bond_style: str = "harmonic",
        angle_style: str = "harmonic",
        dihedral_style: str = "harmonic",
        improper_style: str = "cvff",
        special_bonds: str = "amber",
        pair_modify_shift: str = "no",
        pair_modify_mix: str = "sixthpower",
        kspace_style: str = "pppm",
        kspace_accuracy: float = 1.0e-4,
        kspace_diff: str = "ad",
    ):
        self.lmp = lmp
        self.pair_style = pair_style
        self.pair_cutoff_inner = pair_cutoff_inner
        self.pair_cutoff_outer = pair_cutoff_outer
        self.bond_style = bond_style
        self.angle_style = angle_style
        self.dihedral_style = dihedral_style
        self.improper_style = improper_style

        # Set LAMMPS styles
        self.lmp.cmd.pair_style(pair_style, pair_cutoff_inner, pair_cutoff_outer)
        self.lmp.cmd.pair_modify("shift", pair_modify_shift, "mix", pair_modify_mix)
        self.lmp.cmd.kspace_style(kspace_style, kspace_accuracy)
        self.lmp.cmd.kspace_modify("diff", kspace_diff)
        self.lmp.cmd.bond_style(bond_style)
        self.lmp.cmd.angle_style(angle_style)
        self.lmp.cmd.dihedral_style(dihedral_style)
        self.lmp.cmd.improper_style(improper_style)
        self.lmp.cmd.special_bonds(special_bonds)
