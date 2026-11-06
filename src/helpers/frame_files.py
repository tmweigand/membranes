import os
import glob
import re
from membranes.helpers import helpers


def get_frame_files(
    pattern=None, use_helpers_get_files=False, machine=None, file_type=None
):
    """Get frame files using the original helper by default, or a custom glob."""
    if use_helpers_get_files:
        files = helpers.get_files(machine=machine, file_type=file_type)
    else:
        if not pattern:
            raise ValueError(
                "Set frame_files_pattern when use_helpers_get_files=False."
            )
        files = [path for path in glob.glob(pattern) if os.path.isfile(path)]

    # Sort numerically by the final number in the filename (e.g. timestep).
    def sort_key(path):
        name = os.path.basename(path)
        numbers = re.findall(r"\d+", name)
        if numbers:
            return (0, int(numbers[-1]), name)
        return (1, 0, name)

    files = sorted(files, key=sort_key)
    if not files:
        source = (
            f"helpers.get_files(machine={machine!r}, file_type={file_type!r})"
            if use_helpers_get_files
            else f"glob pattern {pattern!r}"
        )
        raise FileNotFoundError(f"No frame files found using {source}.")
    return files


def safe_label(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("_")
