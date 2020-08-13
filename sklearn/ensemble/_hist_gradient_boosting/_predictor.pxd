# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
# distutils: language=c++

from .common cimport node_struct
from libcpp.vector cimport vector

cdef class TreePredictor:
    # Python access to nodes will create a new container and copy the data into
    # it.
    # https://cython.readthedocs.io/en/latest/src/userguide/wrapping_CPlusPlus.html#standard-library
    cdef public vector[node_struct] nodes
