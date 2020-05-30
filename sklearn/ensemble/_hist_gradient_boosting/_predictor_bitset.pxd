# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
# cython: nonecheck=False
# distutils: language=c++
from cython cimport floating
from libcpp.set cimport set
from libcpp.map cimport map
from .common cimport BITSET_INNER_DTYPE_C
from .common cimport X_BINNED_DTYPE_C

cdef class PredictorBitSet:
    cdef map[int, set[int]] node_to_raw_bitset
    cdef map[int, BITSET_INNER_DTYPE_C] node_to_binned_bitset

    cdef unsigned char raw_category_in_bitset(self, unsigned int node_idx,
                                              floating category) nogil

    cdef unsigned char binned_category_in_bitset(self, unsigned int node_idx, X_BINNED_DTYPE_C category) nogil
