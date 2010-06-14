"""
Wrapper for minilearn
http://fabianp.github.com/minilearn/
"""

import numpy as np
cimport numpy as np

cdef extern from "minilearn.h":
    void lars_fit (int, int, double *, double *, double *, int *, double *, int)
    # int RowMajorStrg

def lars_fit_wrap (np.ndarray[np.float64_t, ndim=2, mode='c'] X,
                   np.ndarray[np.float64_t, ndim=1, mode='c'] res,
                   np.ndarray[np.float64_t, ndim=1, mode='c'] beta,
                   np.ndarray[np.int32_t,   ndim=1, mode='c'] ind,
                   np.ndarray[np.float64_t, ndim=1, mode='c'] chol,
                   np.npy_intp niter):
    """
    res will be rewritten with the residual at output
    """
    lars_fit (<int> X.shape[1], <int> X.shape[0], <double *> X.data,
              <double *> res.data, <double *> beta.data, <int *>
              ind.data, <double *> chol.data, <int> niter)

