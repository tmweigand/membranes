from ..box import Box
from ..regions import RegionClassifier
from .coordination_stats import CoordinationStats
from ..input_config import InputConfig

import numpy as np
from MDAnalysis.lib.distances import self_capped_distance


def coordination_numbers(pair_indices, pair_distances, cutoff_nm, atom_count):
    cutoff_angstrom = np.nextafter(float(cutoff_nm) * 10.0, np.inf)
    selected_pairs = pair_indices[pair_distances <= cutoff_angstrom]
    coordination = np.zeros(atom_count, dtype=np.int64)

    if len(selected_pairs):
        np.add.at(coordination, selected_pairs[:, 0], 1)
        np.add.at(coordination, selected_pairs[:, 1], 1)

    return coordination


class CoordinationNumber:
    """Computes coordination numbers frame by frame and feeds the stats."""

    def __init__(self, config: InputConfig):
        self.config = config
        sim = config.sim
        self.oxygen_type = sim.oxygen_type
        self.region_bounds = sim.region_bounds
        self.frames_processed = 0
        self.box = Box(domain=config.box)
        self.regions = RegionClassifier(self.region_bounds, self.box)

        start, stop = sim.cutoffs
        step = sim.cutoff_step
        if start <= 0 or stop < start:
            raise ValueError("Invalid coordination cutoff range.")
        if step <= 0:
            raise ValueError("cutoff_step must be positive.")

        cutoffs = np.arange(start, stop + step / 2, step)
        self.stats = CoordinationStats(cutoffs, self.region_bounds)

    @property
    def max_cutoff_angstrom(self):
        return float(self.stats.cutoffs_angstrom.max())

    def find_pairs(self, positions):
        """Neighbor search only: all pairs within the largest cutoff."""
        # Small pad so pairs exactly at the max cutoff are not lost to round-off.
        return self_capped_distance(
            positions,
            max_cutoff=self.max_cutoff_angstrom + 1e-6,
            box=self.box.MD_box_from_domain(),
            return_distances=True,
        )

    def accumulate_frame(self, pairs, distances, regions, atom_count):
        """Compute CN at every cutoff and add it to the stats, per region."""
        for i, cutoff in enumerate(self.stats.cutoffs_angstrom):
            cn = coordination_numbers(pairs, distances, cutoff, atom_count)
            for region in self.region_bounds:
                self.stats.update(i, region, cn[regions == region])

    def process_frame(self, atom_data):
        """Full pipeline for one frame. Returns the number of pairs found."""
        positions = atom_data["coords"]

        labels = self.regions.classify_regions(positions[:, 2])

        pairs, distances = self.find_pairs(positions)
        self.accumulate_frame(pairs, distances, labels, len(positions))
        self.frames_processed += 1
        return len(pairs)
