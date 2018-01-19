from libc.math cimport fabs
from types cimport floating, complexing

# complex absolute-value and real-part functions
cdef extern from "complex.h" nogil:
   double cabs(double complex)
   float cabsf(float complex)

cdef floating fmax(floating x, floating y) nogil
cdef floating fmax_arr(int n, floating* a) nogil
cdef double abs_max(int n, complexing* a) nogil
cdef double diff_abs_max(int n, complexing* a, complexing* b) nogil
cdef void relu(int n, floating *x) nogil
cdef double l1_norm(int n, complexing *w, int inc) nogil


