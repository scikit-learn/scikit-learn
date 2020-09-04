import numpy as np
from .common import HISTOGRAM_DTYPE


class HistogramPool:
    """Histogram pool to be used by the growers.

    The pool allocates and returns histograms to the caller. When the `reset`
    method is called, all the previously allocated histograms will be available
    for the next grower to use. New histograms will be allocated when there are
    no histograms left in the available_pool. HistogramPool is used for memory
    allocation/management only. The computation of the histograms is done in
    HistogramBuilder.

    Empirically, this strategy has proven to be more memory efficient (on macOS
    in particular) than allocating new histograms each time and relying on the
    Python GC to keep the overall python process memory low when fitting GBRT
    on datasets with many features and target classes.
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
        """Return a non-initialized array of histogram for one grower node."""
        if self.available_pool:
            histograms = self.available_pool.pop()
        else:
            histograms = np.empty(
                shape=(self.n_features, self.n_bins), dtype=HISTOGRAM_DTYPE
            )
        self.used_pool.append(histograms)
        return histograms

    def release(self, histograms):
        """Move a specific histogram array to the available pool"""
        try:
            idx = next(idx for idx, h in enumerate(self.used_pool)
                       if h is histograms)
        except StopIteration as e:
            raise ValueError("Could not find histograms in used_pool") from e
        self.used_pool.pop(idx)
        self.available_pool.append(histograms)
