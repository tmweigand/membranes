import numpy as np


class Box:
    def __init__(self, domain):
        self.domain = np.asarray(domain, dtype=np.float64)
        if self.domain.shape != (3, 2):
            raise ValueError("Domain must be a 3x2 array representing the box corners.")

    def domain_lengths(self):
        """Get the domain lengths in each dimension."""
        return np.diff(self.domain, axis=1).flatten()

    def MD_box_from_domain(self):
        """Get the domain lengths in each dimension and append angles."""
        lengths = self.domain_lengths()
        return np.concatenate((lengths, [90.0, 90.0, 90.0]))
