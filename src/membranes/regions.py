import numpy as np
from .box import Box


class RegionClassifier:
    def __init__(self, bounds, box: Box):
        self.bounds = bounds
        self.box = box
        self.check_bounds()

    def check_bounds(self):
        """Check if the region bounds are valid and do not overlap."""

        for region, (low, high) in self.bounds.items():
            if self.box.domain[2, 0] > low or self.box.domain[2, 1] < high:
                raise ValueError(
                    f"Region {region} bounds ({low}, {high}) are outside the box domain ({self.box.domain[2, 0]}, {self.box.domain[2, 1]})."
                )
            elif low >= high:
                raise ValueError(
                    f"Invalid bounds for {region}: low ({low}) must be less than high ({high})."
                )
            else:
                overlap = False
                for other_region, (other_low, other_high) in self.bounds.items():
                    if region != other_region:
                        if not (high <= other_low or low >= other_high):
                            overlap = True
                            break
                if overlap:
                    raise ValueError(
                        f"Region {region} bounds ({low}, {high}) overlap with another region's bounds."
                    )
                elif low == high:
                    raise ValueError(
                        f"Region {region} bounds ({low}, {high}) are equal; low must be less than high."
                    )

    def classify_regions(self, z_values):
        labels = np.full(len(z_values), "", dtype=object)

        for region, (low, high) in self.bounds.items():
            mask = (z_values > low) & (z_values <= high)

            labels[mask] = region

        return labels
