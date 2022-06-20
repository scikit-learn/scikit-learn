# Author: jakirkham, Meekail Zain, Thomas Fan

from libc.math cimport isnan, isinf
cimport cython

cpdef enum FiniteStatus:
    all_finite = 0
    has_nan = 1
    has_infinite = 2

def cy_isfinite(cython.floating[::1] a, bint allow_nan=False):
    cdef int result
    with nogil:
        result = _isfinite(a, allow_nan)
    return result

cdef inline int _isfinite(cython.floating[::1] a, bint allow_nan) nogil:
    if allow_nan:
        return _isfinite_allow_nan(a)
    else:
        return _isfinite_disable_nan(a)

cdef inline int _isfinite_allow_nan(cython.floating[::1] a) nogil:
    for i in range(len(a)):
        if isinf(a[i]) != 0:
            return FiniteStatus.has_infinite
    return FiniteStatus.all_finite

cdef inline int _isfinite_disable_nan(cython.floating[::1] a) nogil:
    for i in range(len(a)):
        if isnan(a[i]) != 0:
            return FiniteStatus.has_nan
        if isinf(a[i]) != 0:
            return FiniteStatus.has_infinite
    return FiniteStatus.all_finite
