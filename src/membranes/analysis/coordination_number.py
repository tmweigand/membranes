from ..box import Box
from ..regions import RegionClassifier
from collections import Counter

# from .coordination_stats import CoordinationStats
import numpy as np
from MDAnalysis.lib.distances import self_capped_distance


class CoordinationNumber:
    """Computes coordination numbers frame by frame and feeds the stats."""

    def __init__(
        self,
        cutoffs: float | np.ndarray,
        regions: RegionClassifier | None = None,
    ):
        self.cutoffs = cutoffs
        self.regions = regions

    def find_pairs(self, positions: np.ndarray, cutoff: float, eps: float = 1e-6):
        """Neighbor search only: all pairs within the cutoff."""
        return self_capped_distance(
            positions,
            max_cutoff=cutoff + eps,
            return_distances=True,
        )

    def _get_coordination_numbers(self, positions, pairs, distances, cutoff):
        """Calculate coordination numbers for specified cutoff and return the array."""

        # pairs, distances = self.find_pairs(positions, cutoff)

        selected_pairs = pairs[distances <= cutoff]

        coordination = np.zeros(len(positions), dtype=np.int64)

        if len(selected_pairs):
            np.add.at(coordination, selected_pairs[:, 0], 1)
            np.add.at(coordination, selected_pairs[:, 1], 1)

        return coordination

    def get_coordination_numbers(self, positions):
        if isinstance(self.cutoffs, (float, int)):
            pairs, distances = self.find_pairs(positions, self.cutoffs)
            if self.regions is None:
                return self._get_coordination_numbers(
                    positions, pairs, distances, self.cutoffs
                )

            cn_values = self._get_coordination_numbers(
                positions, pairs, distances, self.cutoffs
            )

            region_dict = {}
            for region in self.regions.bounds:
                region_mask = self.regions.classify_regions(positions[:, 2]) == region
                region_cn_values = cn_values[region_mask]
                region_dict[region] = region_cn_values

            return region_dict
        else:
            pairs, distances = self.find_pairs(positions, max(self.cutoffs))
            coordination_numbers = {}
            if self.regions is None:
                for cutoff in self.cutoffs:
                    cn_values = self._get_coordination_numbers(
                        positions, pairs, distances, cutoff
                    )
                    coordination_numbers[cutoff] = cn_values
                return coordination_numbers

            for cutoff in self.cutoffs:

                cn_values = self._get_coordination_numbers(
                    positions, pairs, distances, cutoff
                )

                for region in self.regions.bounds:
                    region_mask = (
                        self.regions.classify_regions(positions[:, 2]) == region
                    )
                    region_cn_values = cn_values[region_mask]
                    coordination_numbers[(cutoff, region)] = region_cn_values
            return coordination_numbers


# class CoordinationNumberStats(CoordinationNumberRegions):
#     def __init__(self, cutoffs, regions: RegionClassifier):
#         super().__init__(cutoffs, regions)
#         self._data = {
#             (i, region): {"sum": 0, "sum_sq": 0.0, "count": 0, "frequency": Counter()}
#             for i in range(len(self.cutoffs))
#             for region in self.regions.bounds
#         }

#     def update(self, cutoff_index, region, values):
#         values = np.asarray(values, dtype=np.int64)
#         if values.size == 0:
#             return
#         item = self._data[(cutoff_index, region)]
#         item["sum"] += int(values.sum())
#         item["sum_sq"] += float(np.sum(values.astype(np.float64) ** 2))
#         item["count"] += int(values.size)
#         item["frequency"].update(values.tolist())
