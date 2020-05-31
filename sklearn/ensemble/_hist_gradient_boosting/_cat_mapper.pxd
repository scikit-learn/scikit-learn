# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
# cython: nonecheck=False
# distutils: language=c++

from libcpp.map cimport map
from libcpp.vector cimport vector
from .common cimport X_DTYPE_C
from .common cimport X_BINNED_DTYPE_C

cdef class CategoryMapper:
    cdef:
        # raw_category_to_bin[feature_idx] is a map from a raw category to bin
        map[int, map[int, X_BINNED_DTYPE_C]] raw_category_to_bin
        X_BINNED_DTYPE_C missing_values_bin_idx

    cdef X_BINNED_DTYPE_C map_to_bin(self, int feature_idx,
                                     X_DTYPE_C raw_category) nogil
