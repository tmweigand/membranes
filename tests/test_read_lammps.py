from membranes.io import read_lammps
import pytest
from types import SimpleNamespace
import os
import numpy as np


def test_safe_label():
    assert read_lammps.safe_label("region 1") == "region_1"
    assert read_lammps.safe_label("feed/1.0") == "feed_1.0"
    assert read_lammps.safe_label("region-1.2") == "region-1.2"


def test_get_timestep():
    assert read_lammps.get_timestep("pressuredata.71745000.gz") == 71745000
    assert read_lammps.get_timestep("pressuredata.71745000") == 71745000
    assert read_lammps.get_timestep("71745000.*") == 71745000


def test_get_timestep_no_number():
    with pytest.raises(ValueError, match="Could not identify timestep"):
        read_lammps.get_timestep("pressuredata.dump")


def test_check_input_path(tmp_path):

    file = tmp_path / "frame.dump"
    file.touch()

    read_lammps._check_input_path(str(file))


def test_check_input_path_missing(tmp_path):

    file = tmp_path / "missing.dump"

    with pytest.raises(FileNotFoundError, match="File not found"):
        read_lammps._check_input_path(str(file))


def test_get_frame_files(tmp_path):

    (tmp_path / "frame_300.dump").touch()
    (tmp_path / "frame_100.dump").touch()
    (tmp_path / "frame_200.dump").touch()

    config = SimpleNamespace(
        input_path=str(tmp_path),
        file_pattern="*",
    )

    files = read_lammps.get_frame_files(config)

    assert [os.path.basename(file) for file in files] == [
        "frame_100.dump",
        "frame_200.dump",
        "frame_300.dump",
    ]


def test_get_frame_files_no_files(tmp_path):

    config = SimpleNamespace(
        input_path=str(tmp_path),
        file_pattern="frame_*.dump",
    )

    with pytest.raises(FileNotFoundError, match="No frame files found"):
        read_lammps.get_frame_files(config)


def test_group_coordinates():

    atom_data = {
        "x": np.array([1.0, 2.0]),
        "y": np.array([3.0, 4.0]),
        "z": np.array([5.0, 6.0]),
        "type": np.array([1, 2]),
    }

    result = read_lammps._group_coordinates(atom_data)

    np.testing.assert_array_equal(
        result["positions"], [[1.0, 3.0, 5.0], [2.0, 4.0, 6.0]]
    )

    assert "positions" in result


def test_select_atoms_by_type():

    atom_data = {
        "id": np.array([1, 2, 3, 4]),
        "type": np.array([1, 2, 1, 2]),
        "coords": np.array(
            [
                [0, 0, 1],
                [0, 0, 2],
                [0, 0, 3],
                [0, 0, 4],
            ]
        ),
        "charge": np.array([-0.8, 0.4, -0.8, 0.4]),
    }

    result = read_lammps._select_atoms_by_type(atom_data, atom_type=1)

    assert np.array_equal(result["id"], [1, 3])
    assert np.array_equal(result["type"], [1, 1])
    assert np.array_equal(result["coords"], [[0, 0, 1], [0, 0, 3]])
    assert np.array_equal(result["charge"], [-0.8, -0.8])

    atom_data = {
        "ids": np.array([1, 2, 3]),
        "types": np.array([1, 2, 1]),
        "x": np.array([0, 0, 0]),
        "y": np.array([0, 0, 0]),
        "z": np.array([1, 2, 3]),
    }
    atom_data = read_lammps._group_coordinates(atom_data)
    result = read_lammps._select_atoms_by_type(atom_data, atom_type=1)

    assert np.array_equal(result["ids"], [1, 3])
    assert np.array_equal(result["types"], [1, 1])
    assert np.array_equal(result["positions"], [[0, 0, 1], [0, 0, 3]])


def test_read_lammps_file(monkeypatch):

    file = "pressuredata.750000"
    config = object()
    extract_atom = 1

    atom_data = {"positions": np.array([[1, 2, 3]])}
    domain = np.array([[0, 10], [0, 10], [0, 10]])
    timestep = 100

    def mock_read_frame_files(file_arg, config_arg, extract_atom_arg):
        assert file_arg == file
        assert config_arg is config
        assert extract_atom_arg == extract_atom

        return atom_data, domain, timestep

    monkeypatch.setattr(
        read_lammps,
        "_read_frame_files",
        mock_read_frame_files,
    )

    result, result_domain, result_timestep = read_lammps.read_lammps_file(
        file, config, extract_atom
    )
    assert result is atom_data
    assert np.array_equal(result["positions"], np.array([[1, 2, 3]]))
    assert np.array_equal(result_domain, domain)
    assert result_timestep == timestep
