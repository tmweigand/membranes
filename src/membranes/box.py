import numpy as np

__all__ = ["MD_box_from_domain", "domain_length", "validate_box"]


def MD_box_from_domain(domain):
    """Get the domain lengths in each dimension."""
    lengths = domain_length(domain)
    return np.concatenate((lengths, [90.0, 90.0, 90.0]))


def domain_length(domain):
    """Get the domain lengths in each dimension."""
    return np.diff(domain, axis=1).flatten()


def validate_box(box, frame_label):
    box = np.asarray(box, dtype=np.float64)
    if box.size < 6 or not np.all(np.isfinite(box[:6])):
        raise ValueError(f"Invalid box dimensions in {frame_label}: {box}")
    if np.any(box[:3] <= 0):
        raise ValueError(f"Invalid box lengths in {frame_label}: {box[:3]}")
    return box[:6]
