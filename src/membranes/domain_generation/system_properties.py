"""system_properties.py"""


class SystemProperties:
    """System properties for LAMMPS simulation"""

    def __init__(
        self,
        lmp: "lammps",
        dimension: int = 3,
        n_atom_types: int = 1,
        boundary: tuple[str, ...] = ("p", "p", "p"),
    ):
        self.lmp = lmp
        self.dimension = dimension
        self.n_atom_types = n_atom_types
        self.boundary = boundary
        self.lmp.cmd.dimension(dimension)
        self.lmp.cmd.boundary(boundary[0], boundary[1], boundary[2])

    def set_box(self, box_length: list[list[int, int], ...], box_name: str = "box"):
        """
        Set the size of the simulation box
        """
        assert len(box_length) == self.dimension

        if len(box_length) == 3:
            self.lmp.cmd.region(
                box_name,
                "block",
                box_length[0][0],
                box_length[0][1],
                box_length[1][0],
                box_length[1][1],
                box_length[2][0],
                box_length[2][1],
            )

        self.lmp.cmd.create_box(self.n_atom_types, box_name)
