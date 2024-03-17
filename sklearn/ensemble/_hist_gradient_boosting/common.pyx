import numpy as np

cimport cython
from ...utils._typedefs cimport uint32_t

# Y_DYTPE is the dtype to which the targets y are converted to. This is also
# dtype for leaf values, gains, and sums of gradients / hessians. The gradients
# and hessians arrays are stored as floats to avoid using too much memory.
Y_DTYPE = np.float64
X_DTYPE = np.float64
# Potential mix of uint8 (max_bins = 256) and uint16, therefore we set it to uint16,
# see BinnedData.
X_BINNED_DTYPE = np.uint16
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


@cython.final
cdef class Bitsets:
    """An extension type (class) for the bitsets of all features together.

    This class only allocates the smallest possible 1d ndarray and provides access
    to it like a 2d jagged array. It is comparable to a jagged/ragged array.

    Attributes
    ----------
    Public, i.e. accessible from Python:

        offsets : ndarray of shape (n_features + 1), dtype=np.uint32
            The offsets specify which partition of the bitsets ndarray belongs
            to which features: feature j goes from `bitsets[offsets[j]]` until
            `bitsets[offsets[j + 1] - 1]`. `offsets[n_features + 1]` gives
            the total number of base bitsets (X_BITSET_INNER_DTYPE) over all features.
        bitsets : ndarray of shape (offsets[n_features + 1],), \
                dtype=X_BITSET_INNER_DTYPE
            The 1-dimensional array of all bitsets for all features.

    Private, i.e. only accessible from Cython:

        offsets_view : memoryview of `offsets`, dtype=uint32_t
        bitsets_view : memoryview of `bitsets`, dtype=X_BITSET_INNER_DTYPE
        n_features : int
    """

    def __init__(self, offsets):
        self.offsets = offsets
        self.offsets_view = self.offsets

        self.bitsets = np.zeros(
            shape=self.offsets[self.offsets.shape[0] - 1],  # offsets[-1]
            dtype=X_BITSET_INNER_DTYPE,
        )
        self.bitsets_view = self.bitsets

    cdef inline uint32_t n_inner_bitsets(self, int feature_idx) noexcept nogil:
        return self.offsets_view[feature_idx + 1] - self.offsets_view[feature_idx]

    cdef inline BITSET_INNER_DTYPE_C* at(self, int feature_idx) noexcept nogil:
        return &self.bitsets_view[self.offsets_view[feature_idx]]

    def __reduce__(self):
        """Reduce method used for pickling."""
        return (Bitsets, (self.offsets,))


@cython.final
cdef class BinnedData:
    """An extension type as data container for mixed uint8 and uint16 columns.

    Parameters
    ----------
    n_samples : int
        The number of samples of the data.
    n_bins : ndarray of shape (n_features,)
        For each feature the number of bins including the bin for missing values.

    Attributes
    ----------
    Public, i.e. accessible from Python:

        X8 : ndarray of shape (n_samples, n_features_8bit), dtype=np.uint8, F-aligned
            The binned data, all features for which 256 bins are enough.
        X16 : ndarray of shape (n_samples, n_features_bit16), dtype=np.uint16
            The binned data, all features which exceed 256 bins; up to 2**16 = 65536
            bins are supported.
        feature_is_8bit : ndarray of shape (n_features,), dtype=bool (np.uint8)
            Array to indicate if feature is in X8, if false feature is in X16.
        feature_index : ndarray of shape (n_features,), dtype=np.uint32
            Array to map from feature index to corresponding column in either X8 or
            X16.

    Private, i.e. only accessible from Cython:

        X8_view : memoryview of `X8`, dtype=uint8_t
        X16_view : memoryview of `X16`, dtype=uint16_t
        feature_is_8bit_view : memoryview of `feature_is_8bit`, dtype=uint8_t
        feature_index_view : memoryview of `feature_index`, dtype=uint32_t
    """
    def __init__(self, n_samples, n_bins):
        n_features = n_bins.shape[0]
        self.feature_is_8bit = n_bins <= 256
        n_features8 = np.sum(self.feature_is_8bit)
        n_features16 = n_features - n_features8
        self.X8 = np.empty((n_samples, n_features8), dtype=np.uint8, order="F")
        self.X16 = np.empty((n_samples, n_features16), dtype=np.uint16, order="F")
        self.feature_index = np.empty(n_features, dtype=np.int32)
        self.feature_index[self.feature_is_8bit] = np.arange(n_features8)
        self.feature_index[~self.feature_is_8bit] = np.arange(n_features16)

        self.X8_view = self.X8
        self.X16_view = self.X16
        self.feature_is_8bit_view = self.feature_is_8bit
        self.feature_index_view = self.feature_index

    @property
    def shape(self):
        return (self.X8.shape[0], self.feature_is_8bit.shape[0])

    def __array__(self, dtype=None, copy=None):
        if np.all(self.feature_is_8bit):
            dtype = np.uint8
        else:
            dtype = np.uint16
        X = np.empty(self.shape, dtype=dtype)
        for j in range(X.shape[1]):
            X[:, j] = self[:, j]
        return X

    def __eq__(self, x):
        x = np.asarray(x)
        r = np.empty(self.shape, dtype=bool, order="F")
        if x.ndim == 0:
            for j in range(self.shape[1]):
                r[:, j] = (self[:, j] == x)
        elif x.ndim == 1:
            for j in range(self.shape[1]):
                r[:, j] = (self[:, j] == x[j])
        else:
            for j in range(self.shape[1]):
                r[:, j] = (self[:, j] == x[:, j])
        return r

    def __getitem__(self, x):
        if isinstance(x, tuple):
            i, j = x
            f_idx = self.feature_index[j]
            if self.feature_is_8bit[j]:
                return self.X8[i, f_idx]
            else:
                return self.X16[i, f_idx]
        else:
            from exception import NotImplementedError
            raise NotImplementedError()

    @staticmethod
    def from_array(array):
        """For testing only."""
        if array.dtype == np.uint8:
            b = BinnedData(array.shape[0], np.max(array, axis=0))
            b.X8[:] = array
        elif array.dtype == np.uint16:
            b = BinnedData(array.shape[0], np.maximum(256 + 1, np.max(array, axis=0)))
            b.X16[:] = array
        else:
            raise ValueError("Wrong dtype of array.")
        return b

    cdef inline uint8_t[::1] get_feature_view8(self, int feature_idx) noexcept nogil:
        return self.X8_view[:, self.feature_index_view[feature_idx]]

    cdef inline uint16_t[::1] get_feature_view16(self, int feature_idx) noexcept nogil:
        return self.X16_view[:, self.feature_index_view[feature_idx]]

    cdef inline uint8_t get_item8(self, int i, int feature_idx) noexcept nogil:
        return self.X8_view[i, self.feature_index_view[feature_idx]]

    cdef inline uint16_t get_item16(self, int i, int feature_idx) noexcept nogil:
        return self.X16_view[i, self.feature_index_view[feature_idx]]
