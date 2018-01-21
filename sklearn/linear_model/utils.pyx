from cython cimport floating
from blas_api cimport fused_dot


cdef void real_part(complexing x, floating *y) nogil:
    if complexing is float or complexing is double:
        y[0] = x
    elif complexing is complex:
        y[0] = creal(x)
    else:
        y[0] = crealf(x)


cdef inline floating fmax(floating x, floating y) nogil:
    """Max of two real numbers"""
    if x < y:
        return y
    else:
        return x


cdef floating arr_max(int n, floating *a, int inca) nogil:
    cdef int i
    cdef floating m = a[0]
    for i in range(1, n):
        i *= inca
        m = fmax(m, a[i])
    return m


cdef floating abs_max(int n, floating *a, int inca) nogil:
    """b = np.max(np.abs(a))"""
    cdef int i
    cdef floating m = fabs(a[0])
    for i in range(1, n):
        i *= inca
        m = fmax(m, fabs(a[i]))
    return m


cdef floating diff_abs_max(int n, floating *a, int inca, floating *b,
                           int incb) nogil:
    """c = np.max(np.abs(a - b))"""
    cdef int i
    cdef double m = fabs(a[0] - b[0])
    for i in range(1, n):
        m = fmax(m, fabs(a[i * inca] - b[i * incb]))
    return m


cdef inline void relu(int n, floating *x, int incx) nogil:
    cdef int i
    for i in range(n):
        i *= incx
        if x[i] < 0.:
            x[i] = 0.


cdef floating fused_nrm2_squared(int N, floating *X, int incX) nogil:
    """Computes squared L2 norm of X"""
    return fused_dot(N, X, incX, X, incX)


cdef inline floating fsign(floating x) nogil:
    if x == 0:
        return 0
    elif x > 0:
        return 1.
    else:
        return -1.
