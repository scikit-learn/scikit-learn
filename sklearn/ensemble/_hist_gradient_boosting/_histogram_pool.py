import numpy as np
from .common import HISTOGRAM_DTYPE


class HistogramsPool:
    """Histograms pool to be used by the growers."""
    def __init__(self, n_features, n_bins):
        self.n_features = n_features
        self.n_bins = n_bins
        self.pool = []
        self.used_indices = set()
        self.available_indices = set()

    def reset(self):
        """Reset the pool."""
        # Only set histograms that were used to zero
        for used_idx in self.used_indices:
            self.pool[used_idx].fill(0)

        self.used_indices = set()
        self.available_indices = set(range(len(self.pool)))

    def get_new_histograms(self):
        """Return a histograms and its location in the pool."""
        if self.available_indices:
            idx = self.available_indices.pop()
            histograms = self.pool[idx]
        else:
            histograms = np.zeros(
                shape=(self.n_features, self.n_bins), dtype=HISTOGRAM_DTYPE
            )
            self.pool.append(histograms)
            idx = len(self.pool) - 1

        self.used_indices.add(idx)
        return idx, histograms

    def __getitem__(self, idx):
        """Get element from the pool."""
        return self.pool[idx]
