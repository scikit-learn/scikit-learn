import numpy as np
from .common import HISTOGRAM_DTYPE


class HistogramsPool:
    """Histograms pool to be used by the growers."""
    def __init__(self, n_features, n_bins):
        self.n_features = n_features
        self.n_bins = n_bins
        self.pool = []
        self.unused_histograms_indices = set()

    def reset(self):
        """Reset the pool."""
        self.unused_histograms_indices = set(range(len(self.pool)))
        for histograms in self.pool:
            histograms.fill(0)

    def get_new_histograms(self):
        """Return a histograms and its location in the cache."""
        if self.unused_histograms_indices:
            index = self.unused_histograms_indices.pop()
            return index, self.pool[index]

        # all histograms are used
        histograms = np.zeros(
            shape=(self.n_features, self.n_bins), dtype=HISTOGRAM_DTYPE
        )
        self.pool.append(histograms)
        return len(self.pool) - 1, histograms

    def __getitem__(self, idx):
        """Get element from the pool."""
        return self.pool[idx]
