cimport numpy as cnp

from libcpp.vector cimport vector
from ..utils._typedefs cimport ITYPE_t, DTYPE_t, INT32TYPE_t, INT64TYPE_t

ctypedef fused vector_typed:
    vector[DTYPE_t]
    vector[ITYPE_t]
    vector[INT32TYPE_t]
    vector[INT64TYPE_t]

cdef cnp.ndarray vector_to_nd_array(vector_typed * vect_ptr)
