cimport numpy as np
from libcpp.vector cimport vector

from ..utils._typedefs cimport DTYPE_t, INT32TYPE_t, INT64TYPE_t, ITYPE_t

ctypedef fused vector_typed:
    vector[DTYPE_t]
    vector[ITYPE_t]
    vector[INT32TYPE_t]
    vector[INT64TYPE_t]

cdef np.ndarray vector_to_nd_array(vector_typed * vect_ptr)
