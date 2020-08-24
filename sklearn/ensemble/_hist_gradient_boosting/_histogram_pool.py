import numpy as np
from .common import HISTOGRAM_DTYPE


class HistogramsPool:
    """Histograms pool to be used by the growers.

    The pool allocates and returns histograms to the caller. When the `reset`
    method is called, all the previously allocated histograms will be available
    for the next grower to use. New histograms will be allocated when there are
    no histograms left in the available_pool. HistogramsPool is used for memory
    allocation/management only. The computation of the histograms is done in
    HistogramBuilder.
    """
    def __init__(self, n_features, n_bins):
        self.n_features = n_features
        self.n_bins = n_bins
        self.available_pool = []
        self.used_pool = []

    def reset(self):
        """Reset the pool."""
        self.available_pool.extend(self.used_pool)
        self.used_pool = []

    def get(self):
        """Return a non-initialized histogram object."""
        if self.available_pool:
            histograms = self.available_pool.pop()
        else:
            histograms = np.empty(
                shape=(self.n_features, self.n_bins), dtype=HISTOGRAM_DTYPE
            )
        self.used_pool.append(histograms)
        return histograms
