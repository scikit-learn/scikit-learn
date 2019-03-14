# Authors: Arnaud Joly
#
# License: BSD 3 clause
#
# cython: language_level=3


import numpy as np
cimport numpy as np
ctypedef np.npy_uint32 UINT32_t 

cdef enum:
    # Max value for our rand_r replacement (near the bottom).
    # We don't use RAND_MAX because it's different across platforms and
    # particularly tiny on Windows/MSVC.
    RAND_R_MAX = 0x7FFFFFFF

cpdef sample_without_replacement(np.int_t n_population,
                                 np.int_t n_samples,
                                 method=*,
                                 random_state=*)

cdef UINT32_t our_rand_r(UINT32_t* seed) nogil
