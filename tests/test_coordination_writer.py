import membranes
import numpy as np


def test_init_coordination_writer(tmp_path):

    stats = membranes.analysis.CoordinationStats(4.0)

    out_dir = tmp_path / "coordination_output"

    writer = membranes.analysis.CoordinationWriter(
        in_dir="input_data",
        frames_processed=1,
        out_dir=str(out_dir),
        stats=stats,
    )

    assert writer.stats is stats
    assert writer.out_dir == str(out_dir)
    assert out_dir.exists()
    assert out_dir.is_dir()


def test_write_summary(tmp_path):
    cutoff = 4.0
    positions = np.array([[0, 0, 1], [0, 0, 2], [0, 0, 3], [0, 0, 9]])
    stats = membranes.analysis.CoordinationStats(cutoff)
    stats.update(positions)
    writer = membranes.analysis.CoordinationWriter(
        in_dir="input_data",
        frames_processed=1,
        out_dir=str(tmp_path),
        stats=stats,
    )
    writer.write_outputs()
    assert (tmp_path / "coordination_summary.dat").exists()
    with open(tmp_path / "coordination_summary.dat", "r", encoding="utf-8") as f:
        lines = f.readlines()
    assert lines[0].strip() == (
        "# cutoff region mean_coordination_number "
        "std_coordination_number observations"
    )
    assert lines[1].strip() == "4.00 all 1.50000000 0.86602540 4"

    cutoff = [4.0, 10.0]
    regions = {
        "region1": (0, 3),
        "region2": (3, 10),
    }
    positions = np.array([[0, 0, 1], [0, 0, 2], [0, 0, 3], [0, 0, 9]])
    region_classifier = membranes.RegionClassifier(
        bounds=regions, box=membranes.Box(domain=np.array([[0, 10], [0, 10], [0, 10]]))
    )
    stats = membranes.analysis.CoordinationStats(cutoff, region_classifier)
    stats.update(positions)
    writer = membranes.analysis.CoordinationWriter(
        in_dir="input_data",
        frames_processed=1,
        out_dir=str(tmp_path),
        stats=stats,
    )
    writer.write_outputs()
    with open(tmp_path / "coordination_summary.dat", "r", encoding="utf-8") as f:
        lines = f.readlines()
    assert lines[0].strip() == (
        "# cutoff region mean_coordination_number "
        "std_coordination_number observations"
    )
    assert lines[1].strip() == "4.00 region1 2.00000000 0.00000000 3"
    assert lines[2].strip() == "4.00 region2 0.00000000 0.00000000 1"
    assert lines[3].strip() == "10.00 region1 3.00000000 0.00000000 3"
    assert lines[4].strip() == "10.00 region2 3.00000000 0.00000000 1"
    assert len(lines) == 5


def test_write_distributions(tmp_path):
    cutoff = 4.0
    positions = np.array([[0, 0, 1], [0, 0, 2], [0, 0, 3], [0, 0, 9]])
    stats = membranes.analysis.CoordinationStats(cutoff)
    stats.update(positions)
    writer = membranes.analysis.CoordinationWriter(
        in_dir="input_data",
        frames_processed=1,
        out_dir=str(tmp_path),
        stats=stats,
    )
    writer.write_outputs()
    path = tmp_path / "coordination_distribution_4.00.dat"
    assert path.exists()
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    assert lines[0].strip() == ("# coordination_number probability frequency")
    assert lines[1].strip() == "0 0.25 1"
    assert lines[2].strip() == "2 0.75 3"
    assert len(lines) == 3

    cutoff = [4.0, 10.0]
    regions = {
        "region1": (0, 3),
        "region2": (3, 10),
    }
    positions = np.array([[0, 0, 1], [0, 0, 2], [0, 0, 3], [0, 0, 9]])
    region_classifier = membranes.RegionClassifier(
        bounds=regions, box=membranes.Box(domain=np.array([[0, 10], [0, 10], [0, 10]]))
    )
    stats = membranes.analysis.CoordinationStats(cutoff, region_classifier)
    stats.update(positions)
    writer = membranes.analysis.CoordinationWriter(
        in_dir="input_data",
        frames_processed=1,
        out_dir=str(tmp_path),
        stats=stats,
    )
    writer.write_outputs()
    expected = {
        "coordination_distribution_4.00_region1.dat": ["2 1 3"],
        "coordination_distribution_4.00_region2.dat": ["0 1 1"],
        "coordination_distribution_10.00_region1.dat": ["3 1 3"],
        "coordination_distribution_10.00_region2.dat": ["3 1 1"],
    }
    for filename, expected_lines in expected.items():
        path = tmp_path / filename
        assert path.exists()
        with open(path, "r", encoding="utf-8") as f:
            lines = f.readlines()
        assert lines[0].strip() == ("# coordination_number probability frequency")
        assert [line.strip() for line in lines[1:]] == expected_lines
        assert len(lines) == 2


def test_write_readme(tmp_path):
    cutoff = 4.0
    positions = np.array([[0, 0, 1], [0, 0, 2], [0, 0, 3], [0, 0, 9]])
    stats = membranes.analysis.CoordinationStats(cutoff)
    stats.update(positions)
    writer = membranes.analysis.CoordinationWriter(
        in_dir="input_data",
        frames_processed=1,
        out_dir=str(tmp_path),
        stats=stats,
    )
    writer.write_outputs()
    path = tmp_path / "README.txt"
    assert path.exists()
    with open(path, "r", encoding="utf-8") as f:
        content = f.read()

    assert "Water oxygen coordination analysis" in content
    assert "Cutoffs (A): 4.00" in content
    assert "coordination_summary.dat" in content
    assert "coordination_distribution_<cutoff>_<region>.dat" in content
    assert "frames are correlated, so this is not an estimate from" in content
    assert str(tmp_path) in content

    cutoff = [4.0, 10.0]

    regions = {
        "region1": (0, 3),
        "region2": (3, 10),
    }

    region_classifier = membranes.RegionClassifier(
        bounds=regions, box=membranes.Box(domain=np.array([[0, 10], [0, 10], [0, 10]]))
    )

    stats = membranes.analysis.CoordinationStats(cutoff, region_classifier)

    writer = membranes.analysis.CoordinationWriter(
        in_dir="input_data",
        frames_processed=1,
        out_dir=str(tmp_path),
        stats=stats,
    )

    writer.write_outputs()

    path = tmp_path / "README.txt"

    assert path.exists()

    with open(path, "r", encoding="utf-8") as f:
        content = f.read()
        # print(content)

    assert "Water oxygen coordination analysis" in content
    assert "Cutoffs (A): 4.00, 10.00" in content
    assert "region1: (0, 3)" in content
    assert "region2: (3, 10)" in content
    assert "coordination_summary.dat" in content
    assert "coordination_distribution_<cutoff>_<region>.dat" in content
    assert "Standard deviations are population standard deviations" in content
    assert str(tmp_path) in content
