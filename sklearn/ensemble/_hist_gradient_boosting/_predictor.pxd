# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
# distutils: language=c++

from .common cimport node_struct
from libcpp.vector cimport vector

cdef class TreePredictor:
    cdef vector[node_struct] _nodes
    cdef node_struct* get(self, int node_idx) nogil
    cdef int get_size(self) nogil
