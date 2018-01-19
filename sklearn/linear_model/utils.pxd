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
cdef inline floating fsign(floating x) nogil
