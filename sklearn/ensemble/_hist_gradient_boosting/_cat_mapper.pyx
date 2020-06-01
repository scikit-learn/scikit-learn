# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
# cython: nonecheck=False

from .common cimport X_DTYPE_C
from libc.math cimport isnan

cdef class CategoryMapper:

    def __init__(self, X_BINNED_DTYPE_C missing_values_bin_idx):
        self.missing_values_bin_idx = missing_values_bin_idx

    def insert(self, int feature_idx, X_DTYPE_C[:] category_bins):
        cdef int i

        for i in range(category_bins.shape[0]):
            self.raw_category_to_bin[feature_idx][<int>(category_bins[i])] = i

    cdef X_BINNED_DTYPE_C map_to_bin(self, int feature_idx,
                                     X_DTYPE_C raw_category) nogil:
        # negative values and nans are mapped to missing value
        if isnan(raw_category) or raw_category < 0:
            return self.missing_values_bin_idx

        # This should never happen, but to be safe we check for feature_idx in
        # raw_category_to_bin
        if self.raw_category_to_bin.count(feature_idx) == 0:
            return self.missing_values_bin_idx

        cdef:
            int int_value = <int>(raw_category)
            map[int, X_BINNED_DTYPE_C] category_to_bin = \
                self.raw_category_to_bin[feature_idx]

        if category_to_bin.count(int_value) == 0:
            return self.missing_values_bin_idx

        return category_to_bin[int_value]
