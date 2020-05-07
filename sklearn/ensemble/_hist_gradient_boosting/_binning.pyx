# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False
# cython: language_level=3

# Author: Nicolas Hug

cimport cython

import numpy as np
cimport numpy as np
from numpy.math cimport INFINITY
from cython.parallel import prange
from libc.math cimport isnan

from .common cimport X_DTYPE_C, X_BINNED_DTYPE_C

np.import_array()


def _map_num_to_bins(const X_DTYPE_C [:, :] data,
                     list binning_thresholds,
                     const unsigned char missing_values_bin_idx,
                     X_BINNED_DTYPE_C [::1, :] binned):
    """Bin numerical values to discrete integer-coded levels.

    Parameters
    ----------
    data : ndarray, shape (n_samples, n_features)
        The numerical data to bin.
    binning_thresholds : list of arrays
        For each feature, stores the increasing numeric values that are
        used to separate the bins.
    binned : ndarray, shape (n_samples, n_features)
        Output array, must be fortran aligned.
    """
    cdef:
        int feature_idx
        X_DTYPE_C [:] binning_threshold

    for feature_idx in range(data.shape[1]):
        binning_threshold = binning_thresholds[feature_idx]
        # binned_threshold is None when the feature is categorical
        if binning_threshold is not None:
            _map_num_col_to_bins(data[:, feature_idx],
                                binning_threshold,
                                missing_values_bin_idx,
                                binned[:, feature_idx])


cdef void _map_num_col_to_bins(const X_DTYPE_C [:] data,
                               const X_DTYPE_C [:] binning_thresholds,
                               const unsigned char missing_values_bin_idx,
                               X_BINNED_DTYPE_C [:] binned):
    """Binary search to find the bin index for each value in the data."""
    cdef:
        int i
        int left
        int right
        int middle

    for i in prange(data.shape[0], schedule='static', nogil=True):

        if isnan(data[i]):
            binned[i] = missing_values_bin_idx
        else:
            # for known values, use binary search
            left, right = 0, binning_thresholds.shape[0]
            while left < right:
                middle = (right + left - 1) // 2
                if data[i] <= binning_thresholds[middle]:
                    right = middle
                else:
                    left = middle + 1
            binned[i] = left


def _map_cat_to_bins(const X_DTYPE_C [:, :] data,
                     const long [:] categorical_indicies,
                     list bin_catgories,
                     const unsigned char missing_values_bin_idx,
                     X_BINNED_DTYPE_C [::1, :] binned):
    """Encode categories.

    Missing values and unknown values are mapped to the missing bin.

    Parameters
    ----------
    data : ndarray of shape (n_samples, n_features)
        data to encoded.
    categorical_indices : list of int
        columns in ``data`` that are categorical.
    bin_categories : list of arrays
        categories learned during training that corresponds to
        ``categorical_indices``.
    missing_values_bin_idx : uint8
        The index of the bin where missing values are mapped.
    binned : ndarray, shape (n_samples, n_features)
        Output array
    """
    cdef:
        long feature_idx
        X_DTYPE_C [:] categories

    for i, feature_idx in enumerate(categorical_indicies):
        categories = bin_catgories[i]
        _map_cat_col_to_bins(data[:, feature_idx], categories,
                             missing_values_bin_idx, binned[:, feature_idx])


cdef void _map_cat_col_to_bins(const X_DTYPE_C [:] data,
                               const X_DTYPE_C [:] categories,
                               const unsigned char missing_values_bin_idx,
                               X_BINNED_DTYPE_C [:] binned):
        """Binary search for categories."""
        cdef:
            int i, left, right, middle
            unsigned char found
            X_DTYPE_C middle_value, current_value

        for i in prange(data.shape[0], schedule='static', nogil=True):
            if isnan(data[i]) or data[i] < 0:
                binned[i] = missing_values_bin_idx
            else:
                current_value = data[i]
                found = False
                left, right = 0, categories.shape[0] - 1
                while left <= right:
                    middle = (left + right) // 2
                    middle_value = categories[middle]
                    if middle_value < current_value:
                        left = middle + 1
                    elif middle_value > current_value:
                        right = middle - 1
                    else:
                        binned[i] = middle
                        found = True
                        break
                # unknown
                if not found:
                    binned[i] = missing_values_bin_idx


