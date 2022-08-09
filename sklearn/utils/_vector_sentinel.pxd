cimport numpy as cnp

from libcpp.vector cimport vector

ctypedef fused vector_typed:
    vector[cnp.float64_t]
    vector[cnp.intp_t]
    vector[cnp.int32_t]
    vector[cnp.int64_t]

cdef cnp.ndarray vector_to_nd_array(vector_typed * vect_ptr)
