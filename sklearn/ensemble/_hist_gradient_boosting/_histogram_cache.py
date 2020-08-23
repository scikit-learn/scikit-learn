import numpy as np
from .common import HISTOGRAM_DTYPE


class HistogramsCache:
    """Histograms cache to be used by the growers."""
    def __init__(self, n_features, n_bins):
        self.n_features = n_features
        self.n_bins = n_bins
        self.cache = []
        self.unused_histograms_indices = set()

    def reset(self):
        """Reset the cache."""
        self.unused_histograms_indices = set(range(len(self.cache)))
        for histograms in self.cache:
            histograms.fill(0)

    def get_new_histograms(self):
        """Return a histograms and its location in the cache."""
        if self.unused_histograms_indices:
            index = self.unused_histograms_indices.pop()
            return index, self.cache[index]

        # all histograms are used
        histograms = np.zeros(
            shape=(self.n_features, self.n_bins), dtype=HISTOGRAM_DTYPE
        )
        self.cache.append(histograms)
        return len(self.cache) - 1, histograms

    def __getitem__(self, idx):
        """Get element from the cache."""
        return self.cache[idx]
