import membranes
import pytest
import numpy as np


def test_box():
    box = membranes.Box([[0, 1], [0, 1], [0, 1]])
    assert box.domain.shape == (3, 2)

    with pytest.raises(ValueError):
        box2 = membranes.Box([[0, 0, 0], [2, 2]])

    box_lengths = box.domain_lengths()
    assert np.allclose(box_lengths, [1, 1, 1])

    MD_box = box.MD_box_from_domain()
    assert np.allclose(MD_box, [1, 1, 1, 90, 90, 90])
