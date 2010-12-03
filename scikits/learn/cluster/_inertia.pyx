
import numpy as np
cimport numpy as np
cimport cython
ctypedef np.float64_t DOUBLE


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def compute_inertia(np.ndarray[DOUBLE, ndim=1] mi_1,\
                    np.ndarray[DOUBLE, ndim=1] mj_1,\
                    np.ndarray[DOUBLE, ndim=2] mi_2,\
                    np.ndarray[DOUBLE, ndim=2] mj_2,\
                    np.ndarray[DOUBLE, ndim=2] mi_3,\
                    np.ndarray[DOUBLE, ndim=2] mj_3,\
                    np.ndarray[DOUBLE, ndim=1] res):
    cdef int sizemax = mj_3.shape[0]
    cdef int nfeat = mj_3.shape[1]
    cdef int i, j
    cdef double pa, n
    for i in range(sizemax):
        n = mi_1[i] + mj_1[i]
        pa = 0.
        for j in range(nfeat):
            pa += mi_3[i,j] + mj_3[i,j]
            pa -= ((mi_2[i,j] + mj_2[i,j])**2)/n
        res[i] = pa
    return res

