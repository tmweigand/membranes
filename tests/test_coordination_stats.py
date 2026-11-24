import membranes
import numpy as np
from collections import Counter


def test_init_coordination_stats():
    cutoff = 4.0
    regions = {
        "region1": (0, 3),
        "region2": (3, 10),
    }
    region_classifier = membranes.RegionClassifier(
        bounds=regions, box=membranes.Box(domain=np.array([[0, 10], [0, 10], [0, 10]]))
    )
    stats = membranes.analysis.CoordinationStats(cutoff, region_classifier)
    assert (stats._data.keys()) == {
        (4.0, "region1"),
        (4.0, "region2"),
        (4.0, "unassigned"),
    }

    cutoffs = [4.0, 10.0]
    stats = membranes.analysis.CoordinationStats(cutoffs)
    assert (stats._data.keys()) == {(4.0, "all"), (10.0, "all")}

    cutoffs = [4.0, 10.0]
    regions = {
        "region1": (0, 3),
        "region2": (3, 10),
    }
    region_classifier = membranes.RegionClassifier(
        bounds=regions, box=membranes.Box(domain=np.array([[0, 10], [0, 10], [0, 10]]))
    )
    stats = membranes.analysis.CoordinationStats(cutoffs, region_classifier)
    assert (stats._data.keys()) == {
        (4.0, "region1"),
        (4.0, "region2"),
        (4.0, "unassigned"),
        (10.0, "region1"),
        (10.0, "region2"),
        (10.0, "unassigned"),
    }


def test_call_update():
    cutoff = 4.0
    regions = {
        "region1": (0, 3),
        "region2": (3, 10),
    }
    positions = np.array([[0, 0, 1], [0, 0, 2], [0, 0, 3], [0, 0, 9]])
    pos_ids = [1, 2, 3, 4]
    region_classifier = membranes.RegionClassifier(
        bounds=regions, box=membranes.Box(domain=np.array([[0, 10], [0, 10], [0, 10]]))
    )
    stats = membranes.analysis.CoordinationStats(cutoff, region_classifier)
    stats.update(positions, pos_ids)

    assert (stats._data) == {
        (4.0, "region1"): {
            "sum": 6,
            "sum_sq": 12.0,
            "count": 3,
            "frequency": Counter({2: 3}),
        },
        (4.0, "region2"): {
            "sum": 0,
            "sum_sq": 0.0,
            "count": 1,
            "frequency": Counter({0: 1}),
        },
        (4.0, "unassigned"): {
            "sum": 0,
            "sum_sq": 0.0,
            "count": 0,
            "frequency": Counter(),
        },
    }

    cutoff = [4.0, 10.0]
    regions = {
        "region1": (0, 3),
        "region2": (3, 10),
    }
    positions = np.array([[0, 0, 1], [0, 0, 2], [0, 0, 3], [0, 0, 9]])
    pos_ids = [1, 2, 3, 4]

    region_classifier = membranes.RegionClassifier(
        bounds=regions, box=membranes.Box(domain=np.array([[0, 10], [0, 10], [0, 10]]))
    )
    stats = membranes.analysis.CoordinationStats(cutoff, region_classifier)
    stats.update(positions, pos_ids)
    assert (stats._data) == {
        (4.0, "region1"): {
            "sum": 6,
            "sum_sq": 12.0,
            "count": 3,
            "frequency": Counter({2: 3}),
        },
        (4.0, "region2"): {
            "sum": 0,
            "sum_sq": 0.0,
            "count": 1,
            "frequency": Counter({0: 1}),
        },
        (4.0, "unassigned"): {
            "sum": 0,
            "sum_sq": 0.0,
            "count": 0,
            "frequency": Counter(),
        },
        (10.0, "region1"): {
            "sum": 9,
            "sum_sq": 27.0,
            "count": 3,
            "frequency": Counter({3: 3}),
        },
        (10.0, "region2"): {
            "sum": 3,
            "sum_sq": 9.0,
            "count": 1,
            "frequency": Counter({3: 1}),
        },
        (10.0, "unassigned"): {
            "sum": 0,
            "sum_sq": 0.0,
            "count": 0,
            "frequency": Counter(),
        },
    }


def test_summary():
    cutoff = 4.0
    positions = np.array([[0, 0, 1], [0, 0, 2], [0, 0, 3], [0, 0, 9]])
    pos_ids = [2, 3, 1, 4]
    stats = membranes.analysis.CoordinationStats(cutoff)
    stats.update(positions, pos_ids)
    assert stats.summary() == {
        (4.0, "all"): {"count": 4, "mean": 1.5, "std": np.sqrt(0.75)}
    }

    cutoff = [4.0, 10.0]
    regions = {
        "region1": (0, 3),
        "region2": (3, 10),
    }
    positions = np.array([[0, 0, 1], [0, 0, 2], [0, 0, 3], [0, 0, 9]])
    pos_ids = [2, 3, 1, 4]
    region_classifier = membranes.RegionClassifier(
        bounds=regions, box=membranes.Box(domain=np.array([[0, 10], [0, 10], [0, 10]]))
    )
    stats = membranes.analysis.CoordinationStats(cutoff, region_classifier)
    stats.update(positions, pos_ids)

    assert stats.summary() == {
        (4.0, "region1"): {"count": 3, "mean": 2.0, "std": 0.0},
        (4.0, "region2"): {"count": 1, "mean": 0.0, "std": 0.0},
        (4.0, "unassigned"): {"count": 0, "mean": np.nan, "std": np.nan},
        (10.0, "region1"): {"count": 3, "mean": 3.0, "std": 0.0},
        (10.0, "region2"): {"count": 1, "mean": 3.0, "std": 0.0},
        (10.0, "unassigned"): {"count": 0, "mean": np.nan, "std": np.nan},
    }


def test_distribution():
    cutoff = 4.0
    positions = np.array([[0, 0, 1], [0, 0, 2], [0, 0, 3], [0, 0, 9]])
    pos_ids = [1, 2, 3, 4]
    stats = membranes.analysis.CoordinationStats(cutoff)
    stats.update(positions, pos_ids)
    assert stats.distribution() == {
        (4.0, "all"): {
            0: {"count": 1, "probability": 0.25},
            2: {"count": 3, "probability": 0.75},
        }
    }

    cutoff = [4.0, 10.0]
    regions = {
        "region1": (0, 3),
        "region2": (3, 10),
    }
    positions = np.array([[0, 0, 1], [0, 0, 2], [0, 0, 3], [0, 0, 9]])
    pos_ids = [1, 2, 3, 4]
    region_classifier = membranes.RegionClassifier(
        bounds=regions, box=membranes.Box(domain=np.array([[0, 10], [0, 10], [0, 10]]))
    )
    stats = membranes.analysis.CoordinationStats(cutoff, region_classifier)
    stats.update(positions, pos_ids)

    assert stats.distribution() == {
        (4.0, "region1"): {2: {"count": 3, "probability": 1.0}},
        (4.0, "region2"): {0: {"count": 1, "probability": 1.0}},
        (4.0, "unassigned"): {},
        (10.0, "region1"): {3: {"count": 3, "probability": 1.0}},
        (10.0, "region2"): {3: {"count": 1, "probability": 1.0}},
        (10.0, "unassigned"): {},
    }
