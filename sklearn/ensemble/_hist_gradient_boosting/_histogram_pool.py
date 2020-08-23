import weakref

import numpy as np
from .common import HISTOGRAM_DTYPE


class HistogramsPool:
    """Histograms pool to be used by the growers."""
    def __init__(self, n_features, n_bins):
        self.n_features = n_features
        self.n_bins = n_bins
        self.avaliable_pool = []
        self.used_pool = []

    def reset(self):
        """Reset the pool."""
        # Only set histograms that were used to zero
        for histograms in self.used_pool:
            histograms.fill(0)
        self.avaliable_pool.extend(self.used_pool)
        self.used_pool = []

    def get_new_histograms(self):
        """Return a weakref to a histograms."""
        if self.avaliable_pool:
            histograms = self.avaliable_pool.pop()
        else:
            histograms = np.zeros(
                shape=(self.n_features, self.n_bins), dtype=HISTOGRAM_DTYPE
            )
        self.used_pool.append(histograms)
        return weakref.ref(histograms)
