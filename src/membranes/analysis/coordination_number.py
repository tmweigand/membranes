import numpy as np
from MDAnalysis.lib.distances import self_capped_distance


class CoordinationNumber:
    """Computes coordination numbers frame by frame and feeds the stats."""

    def __init__(
        self,
        cutoffs: float | np.ndarray,
    ):
        self.cutoffs = cutoffs

    def find_pairs(self, positions: np.ndarray, cutoff: float, eps: float = 1e-6):
        """Neighbor search only: all pairs within the cutoff."""
        return self_capped_distance(
            positions,
            max_cutoff=cutoff + eps,
            return_distances=True,
        )

    def _get_coordination_numbers(self, positions, pairs, distances, cutoff):
        """Calculate coordination numbers for specified cutoff and return the array."""

        selected_pairs = pairs[distances <= cutoff]

        if len(selected_pairs):
            coordination = np.bincount(selected_pairs.ravel(), minlength=len(positions))
        else:
            coordination = np.zeros(len(positions), dtype=np.int64)

        return coordination

    def get_coordination_numbers(self, positions):
        if isinstance(self.cutoffs, (float, int)):
            cutoffs = [self.cutoffs]
        else:
            cutoffs = self.cutoffs

        pairs, distances = self.find_pairs(positions, max(cutoffs))

        coordination_numbers = {}

        for cutoff in cutoffs:
            coordination_numbers[cutoff] = self._get_coordination_numbers(
                positions, pairs, distances, cutoff
            )

        return coordination_numbers
