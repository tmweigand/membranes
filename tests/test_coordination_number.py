import membranes
import numpy as np


def test_coordination_number():

    cutoff = 4.0
    positions = np.array([[0, 0, 1], [0, 0, 2], [0, 0, 3], [0, 0, 9]])

    cn_analyzer = membranes.analysis.CoordinationNumber(cutoff)

    cn_values = cn_analyzer.get_coordination_numbers(positions)

    assert np.array_equal(cn_values[4.0], np.array([2, 2, 2, 0]))

    cutoff = 10.0
    cn_analyzer = membranes.analysis.CoordinationNumber(cutoff)
    cn_values = cn_analyzer.get_coordination_numbers(positions)

    assert np.array_equal(cn_values[10.0], np.array([3, 3, 3, 3]))


def test_find_pairs():
    cutoff = 4.0
    positions = np.array([[0, 0, 1], [0, 0, 2], [0, 0, 3], [0, 0, 9]])

    cn_analyzer = membranes.analysis.CoordinationNumber(cutoff)

    pairs, distances = cn_analyzer.find_pairs(positions, cutoff)

    assert pairs.shape[1] == 2
    assert np.array_equal(pairs, np.array([[0, 1], [0, 2], [1, 2]]))

    distances_expected = np.array([1.0, 2.0, 1.0])
    assert np.allclose(distances, distances_expected)

    cutoff = 10.0
    pairs, distances = cn_analyzer.find_pairs(positions, cutoff)
    assert np.array_equal(
        pairs, np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])
    )


def test_coordination_numbers_cutoffs():
    cutoffs = [4.0, 10.0]
    positions = np.array([[0, 0, 1], [0, 0, 2], [0, 0, 3], [0, 0, 9]])

    cn_analyzer = membranes.analysis.CoordinationNumber(cutoffs)

    cn_values_dict = cn_analyzer.get_coordination_numbers(positions)

    assert set(cn_values_dict.keys()) == set(cutoffs)

    assert np.array_equal(cn_values_dict[4.0], np.array([2, 2, 2, 0]))
    assert np.array_equal(cn_values_dict[10.0], np.array([3, 3, 3, 3]))
