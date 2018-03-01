# Synopsis: Some fundamental utilities
# Author: Elvis Dohmatob <gmdopp@gmail.Com>

from cython cimport floating

cdef extern from "math.h" nogil:
    double fabs(double x)
    float fabsf(float x)

cdef floating fmax(floating x, floating y) nogil
cdef floating arr_max(int n, floating *X, int incX) nogil
cdef floating abs_max(int n, floating *X, int incX) nogil
cdef floating diff_abs_max(int n, floating* X, int incX, floating* Y,
                           int incY) nogil
cdef void relu(int n, floating *X, int incX) nogil
cdef floating fused_nrm2_squared(int N, floating *X, int incX) nogil
cdef floating fsign(floating x) nogil

cdef floating compute_dual_gap(int n_samples, int n_features, int n_tasks,
                               floating *W_ptr,  # C-order 2d array
                               floating reg, floating l2_reg,
                               floating *X_or_Gram_ptr,  # F-order 2d
                               floating *Y_or_Cov_ptr,  # F-order 2d
                               floating *R_or_RCov_ptr,  # F-order 2d
                               floating *XtA_ptr,  # F-order 2d
                               floating Y_norm2, bint precomputed,
                               int penalty_model, bint positive) nogil
