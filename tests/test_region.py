import membranes
import pytest
import numpy as np


def test_region_classifier():
    # Define a box with a domain
    box = membranes.Box(domain=np.array([[0, 10], [0, 10], [0, 10]]))

    # Define valid region bounds
    valid_bounds = {
        "region1": (1, 3),
        "region2": (3, 5),
        "region3": (5, 7),
        "region4": (7, 8),
    }

    # Create a RegionClassifier instance with valid bounds
    classifier = membranes.RegionClassifier(bounds=valid_bounds, box=box)

    invalid_bounds_overlap = {"region1": (1, 4), "region2": (3, 5)}

    # Check that ValueError is raised for overlapping bounds
    with pytest.raises(ValueError):
        classifier_overlap = membranes.RegionClassifier(
            bounds=invalid_bounds_overlap, box=box
        )

    # Define invalid region bounds (low >= high)
    invalid_bounds_low_high = {"region1": (5, 3)}

    # Create a RegionClassifier instance with invalid low/high bounds

    # Check that ValueError is raised for low >= high
    with pytest.raises(ValueError):
        classifier_low_high = membranes.RegionClassifier(
            bounds=invalid_bounds_low_high, box=box
        )

    out_of_bounds = {"region1": (-1, 3)}
    with pytest.raises(ValueError):
        classifier_out_of_bounds = membranes.RegionClassifier(
            bounds=out_of_bounds, box=box
        )
