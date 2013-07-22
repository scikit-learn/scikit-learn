cimport numpy as cnp
ctypedef cnp.float64_t FLOAT_t
ctypedef cnp.intp_t INT_t
ctypedef cnp.ulong_t INDEX_t
ctypedef cnp.uint8_t BOOL_t
from _basis cimport Basis
from _record cimport PruningPassRecord

cdef class PruningPasser:
    cdef cnp.ndarray X
    cdef cnp.ndarray B
    cdef cnp.ndarray y
    cdef cnp.ndarray weights
    cdef INDEX_t m
    cdef INDEX_t n
    cdef Basis basis
    cdef FLOAT_t penalty
    cdef FLOAT_t sst
    cdef PruningPassRecord record
    
    cpdef run(PruningPasser self)
    
    cpdef PruningPassRecord trace(PruningPasser self)