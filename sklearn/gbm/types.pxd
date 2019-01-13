import numpy as np
cimport numpy as np


ctypedef np.npy_float64 NPY_X_DTYPE
ctypedef np.npy_uint8 NPY_X_BINNED_DTYPE
ctypedef np.npy_float32 NPY_Y_DTYPE

# Same as histogram dtype but we need a struct to declare views. It needs to be
# packed since by default numpy dtypes aren't aligned
cdef packed struct hist_struct:
    float sum_gradients
    float sum_hessians
    unsigned int count
