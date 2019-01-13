import numpy as np
cimport numpy as np


ctypedef np.npy_float64 X_DTYPE_C
ctypedef np.npy_uint8 X_BINNED_DTYPE_C
ctypedef np.npy_float64 Y_DTYPE_C

# Same as histogram dtype but we need a struct to declare views. It needs to be
# packed since by default numpy dtypes aren't aligned
cdef packed struct hist_struct:
    Y_DTYPE_C sum_gradients
    Y_DTYPE_C sum_hessians
    unsigned int count
