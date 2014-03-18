import numpy as np
cimport numpy as np

ctypedef np.npy_float32 DTYPE_t

cdef class CompiledPredictor:
   cdef void* handle
   cdef void* func
