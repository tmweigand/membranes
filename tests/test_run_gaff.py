"""test_run_GAFF.py"""

from membranes.monomers import generate_molecule_data
from membranes.monomers import run_gaff


def test_run_GAFF(tmp_path):
    mpd = generate_molecule_data.GenMolecule(
        "mpd", "Nc1cccc(N)c1", charge=0, forcefield="gaff2", outdir=str(tmp_path)
    )

    mpd_gaff = run_gaff.runGAFF(
        name=mpd.name, smiles=mpd.smiles, charge=mpd.charge, outdir=mpd.outdir
    )

    # # Generate 3d structure
    assert mpd_gaff.mol is not None
    pdb_path = tmp_path / "mpd.pdb"
    assert pdb_path.exists(), "pdb file was not created"

    # Convert to mol2 file
    mpd_gaff.convert_pdb_to_mol2()
    mol2_path = tmp_path / "mpd.mol2"
    assert mol2_path.exists(), "mol2 file was not created"

    # Run antechamber
    mpd_gaff.run_antechamber()
    gaff_path = tmp_path / "mpd_gaff.mol2"
    assert gaff_path.exists(), "gaff file was not created"

    # Run parmchk2
    mpd_gaff.run_parmchk2()
    frcmod_path = tmp_path / "mpd.frcmod"
    assert frcmod_path.exists(), "frcmod file was not created"

    # Run tleap2
    mpd_gaff.run_tleap2()
