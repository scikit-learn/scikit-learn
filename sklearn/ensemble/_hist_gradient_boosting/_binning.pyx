# Author: Nicolas Hug

from cython.parallel import prange
from libc.math cimport isnan

from .common cimport X_DTYPE_C, X_BINNED_DTYPE_C
from ...utils._typedefs cimport uint8_t


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
    const X_DTYPE_C [:] binning_thresholds,
    const uint8_t is_categorical,
    const uint8_t missing_values_bin_idx,
    int n_threads,
    X_BINNED_DTYPE_C [:] binned
):
    """Binary search to find the bin index for each value in the data."""
    cdef:
        int i
        int left
        int right
        int middle

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
            left, right = 0, binning_thresholds.shape[0]
            while left < right:
                # equal to (right + left - 1) // 2 but avoids overflow
                middle = left + (right - left - 1) // 2
                if data[i] <= binning_thresholds[middle]:
                    right = middle
                else:
                    left = middle + 1

            binned[i] = left
