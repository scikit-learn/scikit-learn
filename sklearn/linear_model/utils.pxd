# Synopsis: Some fundamental utilities
# Author: Elvis Dohmatob <gmdopp@gmail.Com>

from types cimport floating, complexing

cdef extern from "math.h" nogil:
    double fabs(double x)
    float fabsf(float x)

# complex absolute-value and real-part functions
cdef extern from "complex.h" nogil:
   double cabs(double complex)
   float cabsf(float complex)
   double creal(double complex)
   float crealf(float complex)

cdef floating fmax(floating x, floating y) nogil
cdef floating real_max(int n, floating *X, int incX) nogil
cdef double abs_max(int n, complexing *X, int incX) nogil
cdef double diff_abs_max(int n, complexing* X, int incX, complexing* Y,
                         int incY) nogil
cdef void relu(int n, floating *X, int incX) nogil
cdef complexing real_part(complexing X) nogil
cdef complexing fused_nrm2_squared(int N, complexing *X, int incX) nogil
