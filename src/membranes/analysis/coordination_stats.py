from collections import Counter
from ..regions import RegionClassifier
from .coordination_number import CoordinationNumber
import numpy as np


class CoordinationStats(CoordinationNumber):
    """Accumulates coordination-number observations per (cutoff, region)."""

    def __init__(self, cutoffs, regions: RegionClassifier | None = None):
        super().__init__(cutoffs, regions)
        self.counter = 0
        if isinstance(self.cutoffs, (float, int)):
            if self.regions is None:
                self._data = {
                    "sum": 0,
                    "sum_sq": 0.0,
                    "count": 0,
                    "frequency": Counter(),
                }
            else:
                self._data = {
                    region: {
                        "sum": 0,
                        "sum_sq": 0.0,
                        "count": 0,
                        "frequency": Counter(),
                    }
                    for region in self.regions.bounds
                }
        else:
            if self.regions is None:
                self._data = {
                    i: {"sum": 0, "sum_sq": 0.0, "count": 0, "frequency": Counter()}
                    for i in self.cutoffs
                }
            else:
                self._data = {
                    (cutoff, region): {
                        "sum": 0,
                        "sum_sq": 0.0,
                        "count": 0,
                        "frequency": Counter(),
                    }
                    for cutoff in self.cutoffs
                    for region in self.regions.bounds
                }

    def update(self, positions):

        values = self.get_coordination_numbers(positions)
        # print(values.keys())
        if isinstance(self.cutoffs, (float, int)) and self.regions is None:
            self._update(self._data, values)
        else:
            for key in values.keys():
                self._update(self._data[key], values[key])

    def _update(self, data, values):
        data["sum"] += int(values.sum())
        data["sum_sq"] += float(np.sum(values.astype(np.float64) ** 2))
        data["count"] += int(values.size)
        data["frequency"].update(values.tolist())

    def summary(self):
        if isinstance(self.cutoffs, (float, int)) and self.regions is None:
            return self._summary(self._data)
        results = {}
        for key, data in self._data.items():
            results[key] = self._summary(data)
        return results

    def _summary(self, data):
        count = data["count"]
        if count == 0:
            return {"count": 0, "mean": np.nan, "std": np.nan}
        mean = data["sum"] / count
        variance = max(data["sum_sq"] / count - mean**2, 0.0)
        return {"count": count, "mean": mean, "std": float(np.sqrt(variance))}

    def distribution(self):

        if isinstance(self.cutoffs, (float, int)) and self.regions is None:
            return self._distribution(self._data)

        results = {}

        for key, data in self._data.items():
            results[key] = self._distribution(data)

        return results

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
