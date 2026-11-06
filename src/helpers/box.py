import numpy as np


def validate_box(box, frame_label):
    box = np.asarray(box, dtype=np.float64)
    if box.size < 6 or not np.all(np.isfinite(box[:6])):
        raise ValueError(f"Invalid box dimensions in {frame_label}: {box}")
    if np.any(box[:3] <= 0):
        raise ValueError(f"Invalid box lengths in {frame_label}: {box[:3]}")
    return box[:6]
