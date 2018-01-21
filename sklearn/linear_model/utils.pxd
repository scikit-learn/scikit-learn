from libc.math cimport fabs
from types cimport floating, complexing

# complex absolute-value and real-part functions
cdef extern from "complex.h" nogil:
   double cabs(double complex)
   float cabsf(float complex)
   double creal(double complex)
   float crealf(float complex)

cdef floating fmax(floating x, floating y) nogil
cdef floating fmax_arr(int n, floating* a, int inca) nogil
cdef void abs_max(int n, complexing *a, int inca, floating *b) nogil
cdef void diff_abs_max(int n, complexing* a, int inca, complexing* b, int incb, floating *c) nogil
cdef void relu(int n, floating *x, int incx) nogil
cdef void real_part(complexing x, floating *y) nogil


