# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import math

from cython.parallel import prange
from libc.math cimport isnan

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

    assert len(binning_thresholds) < 256
    for i in prange(data.shape[0], schedule='static', nogil=True,
                    num_threads=1):
        if (
            isnan(data[i]) or
            # To follow LightGBM's conventions, negative values for
            # categorical features are considered as missing values.
            (is_categorical and data[i] <= 0)
        ):
            binned[i] = missing_values_bin_idx
        else:
            # for known values, use binary search
            binned[i] = _binary_search(data[i], binning_thresholds, len(binning_thresholds))


# For testing
def binary_search(
    X_DTYPE_C value,
    const X_DTYPE_C [::1] binning_thresholds,
    int size,
):
    return _binary_search(value, binning_thresholds, size)


cdef inline int _binary_search(
    X_DTYPE_C value,
    const X_DTYPE_C [::1] binning_thresholds,
    int size,
) nogil:
    cdef:
        int left
        unsigned int i
        int initial_middle
    left = 0

    # Do one special-cased round, to handle sizes that aren't power of 2:
    i = log2ceil(size) - 1
    initial_middle = size - (1 << (log2ceil(size) - 1))
    if binning_thresholds[initial_middle - 1] < value:
        left += initial_middle

    size = (1 << (log2ceil(size) - 1))
    # with gil:
    #     print(value, list(binning_thresholds), initial_middle, left, size)

    # Do the rest with assumption of power of 2:
    while i != 0:
        i -= 1
        size /= 2
        if binning_thresholds[left + size - 1] < value:
            left += size
        # with gil:
        #     print(value, list(binning_thresholds), left, size)

    if value > binning_thresholds[len(binning_thresholds) - 1]:
        left = len(binning_thresholds)
    return left  # min(left, len(binning_thresholds) - 1)


cdef int bins_to_ints[256]
bins_to_ints[0] = 0
for i in range(1, 256):
    bins_to_ints[i] = int(math.ceil(math.log2(i)))


cdef inline unsigned int log2ceil(unsigned int x) nogil:
    return bins_to_ints[x]
