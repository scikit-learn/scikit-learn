# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False
# cython: language_level=3

# Author: Nicolas Hug

cimport cython

import numpy as np
cimport numpy as np
from cython.parallel import prange
from libc.math cimport isnan

from .types cimport X_DTYPE_C, X_BINNED_DTYPE_C

cpdef _map_to_bins(const X_DTYPE_C [:, :] data,
                   list binning_thresholds,
                   const unsigned char [:] support_missing_values,
                   X_BINNED_DTYPE_C [::1, :] binned):
    """Bin numerical values to discrete integer-coded levels.

    Parameters
    ----------
    data : ndarray, shape (n_samples, n_features)
        The numerical data to bin.
    binning_thresholds : list of arrays
        For each feature, stores the increasing numeric values that are
        used to separate the bins.
    support_missing_values: ndarray, shape (n_features,)
        For each feature, indicate whether the first bin is reserved for
        missing values.
    binned : ndarray, shape (n_samples, n_features)
        Output array, must be fortran aligned.
    """
    cdef:
        int feature_idx

    for feature_idx in range(data.shape[1]):

        _map_num_col_to_bins(data[:, feature_idx],
                             binning_thresholds[feature_idx],
                             support_missing_values[feature_idx],
                             binned[:, feature_idx])


cpdef void _map_num_col_to_bins(const X_DTYPE_C [:] data,
                                const X_DTYPE_C [:] binning_thresholds,
                                const unsigned char support_missing_values,
                                X_BINNED_DTYPE_C [:] binned):
    """Binary search to find the bin index for each value in the data."""
    cdef:
        int i
        int left
        int right
        int middle

    for i in prange(data.shape[0], schedule='static', nogil=True):

        if isnan(data[i]):
            # unkown values are mapped to first bin. For this to be correct,
            # support_missing_values must have been correctly set at fit time.
            binned[i] = 0
        else:
            # for known values, use binary search
            left, right = support_missing_values, binning_thresholds.shape[0]
            while left < right:
                middle = (right + left - 1) // 2
                if data[i] <= binning_thresholds[middle]:
                    right = middle
                else:
                    left = middle + 1
            binned[i] = left
