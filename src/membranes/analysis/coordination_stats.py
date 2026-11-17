from collections import Counter
import numpy as np

ANGSTROM_PER_NM = 10.0


class CoordinationStats:
    """Accumulates coordination-number observations per (cutoff, region).

    Cutoffs are identified by their integer index into `cutoffs_angstrom`,
    which avoids floating-point dictionary keys.
    """

    def __init__(self, cutoffs_angstrom, regions):
        self.cutoffs_angstrom = np.asarray(cutoffs_angstrom, dtype=np.float64)
        self.regions = list(regions)
        self._data = {
            (i, region): {"sum": 0, "sum_sq": 0.0, "count": 0, "frequency": Counter()}
            for i in range(len(self.cutoffs_angstrom))
            for region in self.regions
        }

    @property
    def cutoffs_nm(self):
        return self.cutoffs_angstrom / ANGSTROM_PER_NM

    def update(self, cutoff_index, region, values):
        values = np.asarray(values, dtype=np.int64)
        if values.size == 0:
            return
        item = self._data[(cutoff_index, region)]
        item["sum"] += int(values.sum())
        item["sum_sq"] += float(np.sum(values.astype(np.float64) ** 2))
        item["count"] += int(values.size)
        item["frequency"].update(values.tolist())

    def summary(self, cutoff_index, region):
        """Return count, mean, and (population) std for one cutoff/region."""
        item = self._data[(cutoff_index, region)]
        count = item["count"]
        if count == 0:
            return {"count": 0, "mean": np.nan, "std": np.nan}
        mean = item["sum"] / count
        variance = max(item["sum_sq"] / count - mean**2, 0.0)
        return {"count": count, "mean": mean, "std": float(np.sqrt(variance))}

    def distribution(self, cutoff_index, region):
        """Return {cn: {"count": int, "probability": float}} sorted by cn."""
        frequency = self._data[(cutoff_index, region)]["frequency"]
        total = sum(frequency.values())
        return {
            cn: {
                "count": count,
                "probability": count / total if total else np.nan,
            }
            for cn, count in sorted(frequency.items())
        }
