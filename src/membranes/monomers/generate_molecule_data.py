"""
generate_molecule_data.py
"""

import os
from typing import Dict, Optional
import parmed as pmd
from .run_gaff import runGAFF
from .molecule_plotter import MoleculePlotter

from .molecule import convert_to_molecule, Molecule


class GenMolecule:
    """
    Runs the full GAFF pipeline for one molecule and stores the
    local (per-molecule) topology
    """

    def __init__(
        self,
        name: str,
        smiles: str,
        charge: int,
        forcefield: str = "gaff2",
        outdir: str = "/.",
    ):
        self.name = name
        self.smiles = smiles
        self.charge = charge
        self.forcefield = forcefield
        self.outdir = outdir
        self.pmd_structure: Optional[pmd.Structure] = None
        self.mol = None
        self.gaff: Optional[runGAFF] = None

    def generate_forcefield(self, seed: int = None):
        """
        From smiles, generate GAFF
        """
        self.gaff = runGAFF(
            name=self.name, smiles=self.smiles, charge=self.charge, outdir=self.outdir
        )
        self.mol = self.gaff.mol

    def create_pmd(self):
        """
        Load the AMBER files to create a Parmed Structure
        """

        self.pmd_structure = pmd.load_file(
            os.path.join(self.outdir, f"{self.name}.prmtop"),
            xyz=os.path.join(self.outdir, f"{self.name}.inpcrd"),
        )

    def convert_to_molecule(
        self,
        mol_id: int = 1,
        custom_charges: Optional[Dict[int, float]] = None,
    ) -> Molecule:
        if self.pmd_structure is None:
            raise RuntimeError(f"Call .run() before .to_molecule() for {self.name}")

        return convert_to_molecule(
            self.name, self.pmd_structure, mol_id, custom_charges
        )

    def plot_2d_parameters(self, kind="bond", width=900, height=700, out_png=None):
        return MoleculePlotter(
            self.mol, self.pmd_structure, self.name
        ).plot_2d_parameters(kind=kind, width=width, height=height, out_png=out_png)
