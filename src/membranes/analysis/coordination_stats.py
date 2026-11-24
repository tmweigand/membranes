from collections import Counter
from ..regions import RegionClassifier
from .coordination_number import CoordinationNumber
import numpy as np


class CoordinationStats(CoordinationNumber):
    """Accumulates coordination-number observations per (cutoff, region)."""

    def __init__(self, cutoffs, regions: RegionClassifier | None = None):
        super().__init__(cutoffs)
        self.counter = 0
        self.regions = regions
        if isinstance(self.cutoffs, (float, int)):
            cutoffs = [self.cutoffs]
        else:
            cutoffs = self.cutoffs

        if self.regions is None:
            region_names = ["all"]
        else:
            region_names = list(self.regions.bounds) + ["unassigned"]

        self._data = {
            (cutoff, region): {
                "sum": 0,
                "sum_sq": 0.0,
                "count": 0,
                "frequency": Counter(),
            }
            for cutoff in cutoffs
            for region in region_names
        }

    def update(self, positions, atom_ids):

        atom_ids = np.asarray(atom_ids)
        values = self.get_coordination_numbers(positions)

        if self.regions is None:
            region_labels = np.full(len(positions), "all", dtype=object)
        else:
            region_labels = self.regions.classify_regions(positions[:, 2])

        for cutoff, cn_values in values.items():
            for region in self._data:
                current_cutoff, region_name = region
                if current_cutoff != cutoff:
                    continue
                if region_name == "unassigned":
                    region_mask = ~np.isin(region_labels, list(self.regions.bounds))
                else:
                    region_mask = region_labels == region_name

                region_atom_ids = atom_ids[region_mask]
                region_cn_values = cn_values[region_mask]

                self._update(self._data[region], region_cn_values)
        self.counter += 1

    def _update(self, data, values):
        data["sum"] += int(values.sum())
        data["sum_sq"] += float(np.sum(values.astype(np.float64) ** 2))
        data["count"] += int(values.size)
        data["frequency"].update(values.tolist())

    def summary(self):
        return {key: self._summary(data) for key, data in self._data.items()}

    def _summary(self, data):
        count = data["count"]
        if count == 0:
            return {"count": 0, "mean": np.nan, "std": np.nan}
        mean = data["sum"] / count
        variance = max(data["sum_sq"] / count - mean**2, 0.0)
        return {"count": count, "mean": mean, "std": float(np.sqrt(variance))}

    def distribution(self):
        return {key: self._distribution(data) for key, data in self._data.items()}

    def _distribution(self, data):

        frequency = data["frequency"]
        total = sum(frequency.values())

        return {
            cn: {
                "count": count,
                "probability": count / total if total else np.nan,
            }
            for cn, count in sorted(frequency.items())
        }
