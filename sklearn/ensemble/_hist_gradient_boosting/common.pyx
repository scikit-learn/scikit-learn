import numpy as np

cimport cython
from ...utils._typedefs cimport uint32_t

# Y_DYTPE is the dtype to which the targets y are converted to. This is also
# dtype for leaf values, gains, and sums of gradients / hessians. The gradients
# and hessians arrays are stored as floats to avoid using too much memory.
Y_DTYPE = np.float64
X_DTYPE = np.float64
X_BINNED_DTYPE = np.uint8  # hence max_bins == 256
# dtype for gradients and hessians arrays
G_H_DTYPE = np.float32
X_BITSET_INNER_DTYPE = np.uint32

# Note that we use Y_DTYPE=float64 to avoid issues with floating point precision when
# summing gradients and hessians (both float32). Those are difficult to protect via
# tools like (Kahan-) Neumaier summation as in CPython, see
# https://github.com/python/cpython/issues/100425, or pairwise summation as numpy, see
# https://github.com/numpy/numpy/pull/3685, due to the way histograms are summed
# (number of additions per bin is not known in advance). See also comment in
# _subtract_histograms.
HISTOGRAM_DTYPE = np.dtype([
    ('sum_gradients', Y_DTYPE),  # sum of sample gradients in bin
    ('sum_hessians', Y_DTYPE),  # sum of sample hessians in bin
    ('count', np.uint32),  # number of samples in bin
])

PREDICTOR_RECORD_DTYPE = np.dtype([
    ('value', Y_DTYPE),
    ('count', np.uint32),
    ('feature_idx', np.intp),
    ('num_threshold', X_DTYPE),
    ('missing_go_to_left', np.uint8),
    ('left', np.uint32),
    ('right', np.uint32),
    ('gain', Y_DTYPE),
    ('depth', np.uint32),
    ('is_leaf', np.uint8),
    ('bin_threshold', X_BINNED_DTYPE),
    ('is_categorical', np.uint8),
    # The index of the corresponding bitsets in the Predictor's bitset arrays.
    # Only used if is_categorical is True
    ('bitset_idx', np.uint32)
])

ALMOST_INF = 1e300  # see LightGBM AvoidInf()


@cython.final
cdef class Histograms:
    """An extension type (class) for histograms with variable bins.

    This class only allocates the smallest possible 1d ndarray and provides access
    to it like a 2d jagged array. It is comparable to a jagged/ragged array.

    Attributes
    ----------
    Public, i.e. accessible from Python:

        bin_offsets : ndarray of shape (n_features + 1), dtype=np.uint32
            The bin offsets specify which partition of the histograms ndarray belongs
            to which features: feature j goes from `histograms[bin_offsets[j]]` until
            `histograms[bin_offsets[j + 1] - 1]`. `bin_offsets[n_features + 1]` gives
            the total number of bins over all features.
        histograms : ndarray of shape (bin_offsets[n_features + 1],), \
                dtype=HISTOGRAM_DTYPE
            The 1-dimensional array of all histograms for all features.

    Private, i.e. only accessible from Cython:

        bin_offsets_view : memoryview of `bin_offsets`, dtype=uint32_t
        histograms_view : memoryview of `histograms`, dtype=hist_stucts
        n_features : int
    """

    def __init__(self, n_features, bin_offsets):
        self.n_features = n_features
        if isinstance(bin_offsets, int):
            self.bin_offsets = np.zeros(shape=self.n_features + 1, dtype=np.uint32)
            self.bin_offsets[1:] = np.cumsum(np.full(shape=n_features, fill_value=bin_offsets))
        else:
            self.bin_offsets = bin_offsets
        self.bin_offsets_view = self.bin_offsets

        self.histograms = np.empty(
            shape=self.bin_offsets[self.bin_offsets.shape[0] - 1],  # bin_offsets[-1]
            dtype=HISTOGRAM_DTYPE,
        )
        self.histograms_view = self.histograms

    cdef inline uint32_t n_bins(self, int feature_idx) noexcept nogil:
        return self.bin_offsets_view[feature_idx + 1] - self.bin_offsets_view[feature_idx]

    cdef inline hist_struct* at(self, int feature_idx, uint32_t bin_idx) noexcept nogil:
        return &self.histograms_view[self.bin_offsets_view[feature_idx] + bin_idx]

    # From here on, only for better testing.
    def fill_zeros(self):
        cdef:
            uint32_t i
            n = self.histograms_view.shape[0]
        for i in range(n):
            self.histograms_view[i].sum_gradients = 0.
            self.histograms_view[i].sum_hessians = 0.
            self.histograms_view[i].count = 0
        return self

    def copy(self):
        """Return a copy of self."""
        h = Histograms(
            n_features=self.n_features,
            bin_offsets=self.bin_offsets.copy(),
        )
        np.copyto(h.histograms, self.histograms)
        return h

    def feature_histo(self, feature_idx):
        start = self.bin_offsets[feature_idx]
        end = self.bin_offsets[feature_idx + 1]
        return self.histograms[start:end]

    def feature_sum(self, key):
        return np.asarray(
            [np.sum(self.feature_histo(f)[key]) for f in range(self.n_features)]
        )

    def value_at(self, feature_idx, bin_idx):
        return self.histograms[self.bin_offsets[feature_idx] + bin_idx]
