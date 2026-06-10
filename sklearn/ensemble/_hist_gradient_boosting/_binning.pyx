# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from cython.parallel import prange
from libc.math cimport isnan
from libc.stdint cimport uint8_t

from sklearn.ensemble._hist_gradient_boosting.common cimport X_DTYPE_C, X_BINNED_DTYPE_C
from sklearn.utils._typedefs cimport uint8_t


def _map_to_bins(const X_DTYPE_C [:, :] data,
                 list binning_thresholds,
                 const uint8_t[::1] is_categorical,
                 const uint8_t missing_values_bin_idx,
                 int n_threads,
                 X_BINNED_DTYPE_C [::1, :] binned):
    """Bin continuous and categorical values to discrete integer-coded levels.

    A given value x is mapped into bin value i iff
    thresholds[i - 1] < x <= thresholds[i]

    Parameters
    ----------
    data : ndarray, shape (n_samples, n_features)
        The data to bin.
    binning_thresholds : list of arrays
        For each feature, stores the increasing numeric values that are
        used to separate the bins.
    is_categorical : ndarray of uint8_t of shape (n_features,)
        Indicates categorical features.
    n_threads : int
        Number of OpenMP threads to use.
    binned : ndarray, shape (n_samples, n_features)
        Output array, must be fortran aligned.
    """
    cdef:
        int feature_idx

    for feature_idx in range(data.shape[1]):
        _map_col_to_bins(
            data[:, feature_idx],
            binning_thresholds[feature_idx],
            is_categorical[feature_idx],
            missing_values_bin_idx,
            n_threads,
            binned[:, feature_idx]
        )


cdef void _map_col_to_bins(
    const X_DTYPE_C [:] data,
    const X_DTYPE_C [::1] binning_thresholds,
    const uint8_t is_categorical,
    const uint8_t missing_values_bin_idx,
    int n_threads,
    X_BINNED_DTYPE_C [::1] binned
):
    """Binary search to find the bin index for each value in the data."""
    cdef:
        int i

    # If this limit increases, log2ceil will need to be changed to support
    # higher numbers:
    assert len(binning_thresholds) < 256

    for i in prange(data.shape[0], schedule='static', nogil=True,
                    num_threads=n_threads):
        if (
            isnan(data[i]) or
            # To follow LightGBM's conventions, negative values for
            # categorical features are considered as missing values.
            (is_categorical and data[i] < 0)
        ):
            binned[i] = missing_values_bin_idx
        else:
            # for known values, use binary search
            binned[i] = _binary_search(
                data[i],
                binning_thresholds,
                len(binning_thresholds)
            )


cdef inline size_t _binary_search(
    X_DTYPE_C value,
    const X_DTYPE_C [::1] binning_thresholds,
    size_t size,
) noexcept nogil:
    cdef:
        size_t left
        size_t half
        size_t middle
        size_t remaining_size

    # This implementation is designed to minimize branch mispredictions. See:
    # https://pvk.ca/Blog/2012/07/03/binary-search-star-eliminates-star-branch-mispredictions/
    # https://pvk.ca/Blog/2015/11/29/retrospective-on-binary-search-and-on-compression-slash-compilation/
    left = 0
    remaining_size = size

    # Fixed number of loops, instead of less-predictable while loop:
    for _ in range(log2ceil(size)):
        half = remaining_size / 2
        middle = left + half
        # Try for cmov instead of branch; see
        # https://en.algorithmica.org/hpc/pipelining/branchless/ for details:
        left = middle if (binning_thresholds[middle] < value) else left
        remaining_size -= half

    # Try for cmov instead of branch:
    left = left + 1 if (
        (left < size) and (value > binning_thresholds[left])
    ) else left
    return left


# Created with:
#
#  int_to_log2ceil[0] = 0
#  for i in range(1, 256):
#      int_to_log2ceil[i] = int(math.ceil(math.log2(i)))
cdef uint8_t[256] int_to_log2ceil = [
    0, 0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4,
    4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8
]

cdef inline uint8_t log2ceil(uint8_t x) noexcept nogil:
    # Using a lookup table is slightly faster than calculating on demand:
    return int_to_log2ceil[x]
