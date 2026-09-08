"""
run_gaff.py
"""

import os
import subprocess
import random
from rdkit import Chem as rdkit_chem
from rdkit.Chem import AllChem as rdkit_all_chem
from ..utils import check_dependencies

check_dependencies()


class runGAFF:
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
        charge_model: str = "bcc",  # Run antechamber -L for all options!
        outdir: str = "/.",
        seed=None,
    ):
        self.name = name
        self.charge = charge
        self.smiles = smiles
        self.forcefield = forcefield
        self.charge_model = charge_model
        self.outdir = outdir
        self.mol = None

        # Run pipeline with final files of {self.name}.prmtop and {self.name}.inpcrd
        self.generate_3d_structure(seed)
        self.convert_pdb_to_mol2()
        self.run_antechamber()
        self.run_parmchk2()
        self.run_tleap2()

    def generate_3d_structure(self, seed: int = None):
        """
        Generate 3-D structure from smiles
        """
        if seed is None:
            seed = random.randint(0, 2**31 - 1)

        # Generate molecule from smiles
        self.mol = rdkit_chem.MolFromSmiles(self.smiles)
        self.mol = rdkit_chem.AddHs(self.mol)

        # Use ETKDG to compute atomic coordinates in 3D
        p = rdkit_all_chem.ETKDGv3()
        p.randomSeed = seed
        rdkit_all_chem.EmbedMolecule(self.mol, p)
        rdkit_all_chem.MMFFOptimizeMolecule(self.mol)

        # Write output as pdb file
        rdkit_chem.MolToPDBFile(self.mol, os.path.join(self.outdir, f"{self.name}.pdb"))

    def convert_pdb_to_mol2(self):
        """
        Convert the generated PDB file to MOL2 format using Open Babel.

        This step infers bonding and assigns atom types. Open Babel may also
        assign approximate partial charges (e.g., Gasteiger), but these are
        not guaranteed and are typically replaced during parameterization
        (e.g., with antechamber/GAFF).
        """
        subprocess.run(
            ["obabel", f"{self.name}.pdb", "-O", f"{self.name}.mol2"],
            check=True,
            capture_output=True,
            cwd=self.outdir,
        )

    def run_antechamber(self):
        """
        Run AmberTools `antechamber` to assign GAFF/GAFF2 atom types and AM1-BCC charges.

        Input:
            - {name}.mol2 (from Open Babel conversion)

        Output:
            - {name}_gaff.mol2 (typed molecule with charges)

        Command options used:
            -i   input filename
            -fi  input format (`mol2`)
            -o   output filename
            -fo  output format (`mol2`)
            -c   charge method (`bcc` = AM1-BCC)
            -s   status/verbosity level (`2`)
            -at  atom type set (`gaff` or `gaff2`, from self.forcefield)
            -nc  net molecular charge (from self.charge)
            -m   spin multiplicity (`1` for singlet)
        """
        subprocess.run(
            [
                "antechamber",
                "-i",
                f"{self.name}.mol2",
                "-fi",
                "mol2",
                "-o",
                f"{self.name}_gaff.mol2",
                "-fo",
                "mol2",
                "-c",
                self.charge_model,
                "-s",
                "2",
                "-at",
                self.forcefield,
                "-nc",
                str(self.charge),
                "-m",
                "1",
            ],
            check=True,
            capture_output=True,
            cwd=self.outdir,
        )

    def run_parmchk2(self):
        """
        Run AmberTools `parmchk2` to generate missing force-field parameters.

        Purpose:
            `parmchk2` inspects the GAFF-typed MOL2 file and identifies parameters
            not found directly in the selected force field. It writes these guessed
            or supplementary parameters to an `.frcmod` file for use in `tleap`.

        Input:
            - {name}_gaff.mol2 (typically produced by `run_antechamber`)

        Output:
            - {name}.frcmod

        Command options used:
            -i  input filename
            -f  input format (`mol2`)
            -o  output frcmod filename
            -s  force-field family (`gaff`/`gaff2`, from self.forcefield)
        """
        subprocess.run(
            [
                "parmchk2",
                "-i",
                f"{self.name}_gaff.mol2",
                "-f",
                "mol2",
                "-o",
                f"{self.name}.frcmod",
                "-s",
                self.forcefield,
            ],
            check=True,
            capture_output=True,
            cwd=self.outdir,
        )

    def run_tleap2(self):
        """
        Run AmberTools `tleap` to build Amber topology and coordinate files.

        Purpose:
            `tleap` combines:
            - the GAFF/GAFF2 base force field (`leaprc.<forcefield>`)
            - molecule-specific parameters from `{name}.frcmod`
            - typed/charged structure from `{name}_gaff.mol2`
            and writes final Amber simulation inputs.

        Input:
            - {name}_gaff.mol2 (from `run_antechamber`)
            - {name}.frcmod (from `run_parmchk2`)

        Output:
            - {name}.prmtop (Amber topology/parameters)
            - {name}.inpcrd (Amber coordinates)
            - tleap_{name}.in (generated tleap input script)

        tleap script commands used:
            - source leaprc.{self.forcefield}
            - loadamberparams {name}.frcmod
            - MOL = loadmol2 {name}_gaff.mol2
            - check MOL
            - saveamberparm MOL {name}.prmtop {name}.inpcrd
            - quit
        """
        tleap_lines = [
            f"source leaprc.{self.forcefield}",
            f"loadamberparams {self.name}.frcmod",
            f"MOL = loadmol2 {self.name}_gaff.mol2",
            "check MOL",
            f"saveamberparm MOL {self.name}.prmtop {self.name}.inpcrd",
            "quit",
            "",
        ]

        tleap_input = os.path.join(self.outdir, f"tleap_{self.name}.in")
        with open(tleap_input, "w") as fh:
            fh.write("\n".join(tleap_lines))

        # Execute tleap in molecule output directory
        subprocess.run(
            ["tleap", "-f", os.path.basename(tleap_input)],
            check=True,
            capture_output=True,
            cwd=self.outdir,
        )
