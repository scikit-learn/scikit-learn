from types cimport floating, complexing
from blas_api cimport fused_dotc


cdef complexing real_part(complexing x) nogil:
    if complexing is float or complexing is double:
        return x
    elif complexing is complex:
        return creal(x)
    else:
        return crealf(x)


cdef inline floating fmax(floating x, floating y) nogil:
    if x < y:
        return y
    else:
        return x


cdef double abs_max(int n, complexing *a, int inca) nogil:
    """b = np.max(np.abs(a))"""
    if complexing is float or complexing is double:
        abs = fabs
    else:
        abs = cabs
    cdef int i
    cdef double m = abs(a[0])
    for i in range(1, n):
        i *= inca
        m = fmax(m, abs(a[i]))
    return m


cdef floating real_max(int n, floating *a, int inca) nogil:
    """np.max(a.real)"""
    cdef int i
    cdef floating m = a[0]
    for i in range(1, n):
        i *= inca
        m = fmax(m, a[i])
    return m


cdef double diff_abs_max(int n, complexing *a, int inca, complexing *b,
                         int incb) nogil:
    """c = np.max(np.abs(a - b))"""
    if complexing is float or complexing is double:
        abs = fabs
    else:
        abs = cabs
    cdef int i
    cdef double m = abs(a[0] - b[0])
    for i in range(1, n):
        m = fmax(m, abs(a[i * inca] - b[i * incb]))
    return m


cdef inline void relu(int n, floating *x, int incx) nogil:
    cdef int i
    for i in range(n):
        i *= incx
        if x[i] < 0.:
            x[i] = 0.


cdef complexing fused_nrm2_squared(int N, complexing *X, int incX) nogil:
    """Computes squared L2 norm of X"""
    return fused_dotc(N, X, incX, X, incX)
