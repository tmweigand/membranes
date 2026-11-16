import numpy as np


def classify_regions(z_values, bounds):
    labels = np.full(len(z_values), "", dtype=object)

    if len(bounds) == 1:
        only_region = next(iter(bounds))
        low, high = bounds[only_region]
        if np.isneginf(low) and np.isposinf(high):
            return np.full(len(z_values), only_region, dtype=object)

    for region, (low, high) in bounds.items():
        if region == "membrane" and "feed" in bounds and "permeate" in bounds:
            mask = (z_values >= low) & (z_values <= high)
        elif region == "feed" and "membrane" in bounds:
            mask = z_values < high
        elif region == "permeate" and "membrane" in bounds:
            mask = z_values > low
        else:
            if np.isneginf(low):
                mask = z_values <= high
            elif np.isposinf(high):
                mask = z_values > low
            else:
                mask = (z_values > low) & (z_values <= high)

        # If any overlaps occur, raise an error. This ensures that each oxygen atom is assigned to exactly one region.
        overlap = mask & (labels != "")
        if np.any(overlap):
            raise ValueError(
                f"Region definitions overlap for {region!r}; check region_bounds."
            )
        labels[mask] = region

    # if np.any(labels == ""):
    #    examples = z_values[labels == ""][:5]
    #    raise ValueError(
    #        "Some oxygen atoms do not fall into any region. "
    #        f"Example z coordinates: {examples}. Check region_bounds."
    #    )
    return labels
