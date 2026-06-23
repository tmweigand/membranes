"""test_generate_molecule_data.py"""

from membranes.monomers import generate_molecule_data
from membranes.monomers import monomer_data
from membranes.monomers import forcefield


def test_gen_mpd(tmp_path):
    mpd = generate_molecule_data.GenMolecule(
        "mpd", "Nc1cccc(N)c1", charge=0, forcefield="gaff2", outdir=str(tmp_path)
    )

    # Generate GAFF data
    mpd.generate_forcefield()

    # Create Parmed structure
    mpd.create_pmd()

    # Plot molecule with charges
    mpd.plot_2d_parameters(kind="charge", out_png=str(tmp_path / "mpd_charge.png"))

    # Switch from Parmed Structure to Molecule
    mol = mpd.convert_to_molecule()

    # Save as lammps .mol file
    mol.write_mol(tmp_path / "lammps_mpd.mol")

    # Save data as JSON
    mf = monomer_data.MoleculeFile.from_molecule(mol)
    mf.save(tmp_path / "MPD_molecule.json")

    ff = forcefield.ForceField([mol])
    ff.write(tmp_path / "ff.lammps")


def test_gen_tmc(tmp_path):
    tmc = generate_molecule_data.GenMolecule(
        "tmc",
        "ClC(=O)c1cc(C(Cl)=O)cc(C(Cl)=O)c1",
        charge=0,
        forcefield="gaff2",
        outdir=str(tmp_path),
    )

    # Generate GAFF data
    tmc.generate_forcefield()

    # Create Parmed structure
    tmc.create_pmd()

    # Plot molecule with charges
    tmc.plot_2d_parameters(kind="charge", out_png=str(tmp_path / "tmc_charge.png"))

    # Switch from Parmed Structure to Molecule
    mol = tmc.convert_to_molecule()

    # Save as lammps .mol file
    mol.write_mol(tmp_path / "lammps_tmc.mol")

    # Save data as JSON
    mf = monomer_data.MoleculeFile.from_molecule(mol)
    mf.save(tmp_path / "TMC_molecule.json")

    ff = forcefield.ForceField([mol])
    ff.write(tmp_path / "ff.lammps")


def test_gen_mpd_and_tmc(tmp_path):

    mpd = generate_molecule_data.GenMolecule(
        "mpd", "Nc1cccc(N)c1", charge=0, forcefield="gaff2", outdir=str(tmp_path)
    )

    mpd.generate_forcefield()
    mpd.create_pmd()
    mpd_mol = mpd.convert_to_molecule()

    # Save as lammps .mol file
    mpd_mol.write_mol(tmp_path / "lammps_mpd.mol")

    tmc = generate_molecule_data.GenMolecule(
        "tmc",
        "ClC(=O)c1cc(C(Cl)=O)cc(C(Cl)=O)c1",
        charge=0,
        forcefield="gaff2",
        outdir=str(tmp_path),
    )

    tmc.generate_forcefield()
    tmc.create_pmd()
    tmc_mol = tmc.convert_to_molecule()

    # Save as lammps .mol file
    tmc_mol.write_mol(tmp_path / "lammps_tmc.mol")

    # Create LAMMPS output
    ff = forcefield.ForceField([mpd_mol, tmc_mol])
    ff.write(tmp_path / "ff.lammps")


def test_gen_mpd1_tmc1(tmp_path):
    tmc = generate_molecule_data.GenMolecule(
        "tmc",
        "Nc1cccc(NC(=O)c2cc(C(Cl)=O)cc(C(Cl)=O)c2)c1",
        charge=0,
        forcefield="gaff2",
        outdir=str(tmp_path),
    )

    # Generate GAFF data
    tmc.generate_forcefield()

    # Create Parmed structure
    tmc.create_pmd()

    # Plot molecule with charges
    tmc.plot_2d_parameters(
        kind="charge", out_png=str(tmp_path / "mpd1_tmc1_charge.png")
    )

    # Switch from Parmed Structure to Molecule
    mol = tmc.convert_to_molecule()

    # Save as lammps .mol file
    mol.write_mol(tmp_path / "lammps_mpd1_tmc1.mol")

    ff = forcefield.ForceField([mol])
    ff.write(tmp_path / "mpd1_tmc1_ff.lammps")
