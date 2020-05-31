# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
# cython: nonecheck=False
# distutils: language=c++
from libcpp.set cimport set
from libcpp.map cimport map
from libcpp.vector cimport vector
from .common cimport BITSET_INNER_DTYPE_C
from .common cimport X_BINNED_DTYPE_C
from .common cimport X_DTYPE_C

cdef class PredictorBitSet:
    cdef map[int, set[int]] node_to_raw_bitset
    cdef map[int, vector[BITSET_INNER_DTYPE_C]] node_to_binned_bitset
    cdef map[int, set[int]] feature_idx_raw_cats

    cdef unsigned char raw_category_in_bitset(self, unsigned int node_idx,
                                              X_DTYPE_C category) nogil

    cdef unsigned char binned_category_in_bitset(self, unsigned int node_idx,
                                                 X_BINNED_DTYPE_C category) nogil

    cdef unsigned char is_known_category(self, unsigned int feature_idx,
                                         X_DTYPE_C category) nogil
