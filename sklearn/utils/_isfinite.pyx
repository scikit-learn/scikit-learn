# Author: jakirkham, Meekail Zain, Thomas Fan

from libc.math cimport isnan, isinf
cimport cython
from cython cimport floating

cpdef enum FiniteStatus:
    all_finite = 0
    has_nan = 1
    has_infinite = 2

def cy_isfinite(floating[::1] a, bint allow_nan=False):
    cdef int result
    with nogil:
        result = _isfinite(a, allow_nan)
    return result

cdef inline int _isfinite(floating[::1] a, bint allow_nan) nogil:
    cdef floating* a_ptr = &a[0]
    cdef Py_ssize_t length = len(a)
    if allow_nan:
        return _isfinite_allow_nan(a_ptr, length)
    else:
        return _isfinite_disable_nan(a_ptr, length)

cdef inline int _isfinite_allow_nan(floating* a_ptr, Py_ssize_t length) nogil:
    cdef Py_ssize_t i
    for i in range(length):
        if isinf(a_ptr[i]) != 0:
            return FiniteStatus.has_infinite
    return FiniteStatus.all_finite

cdef inline int _isfinite_disable_nan(floating* a_ptr, Py_ssize_t length) nogil:
    cdef Py_ssize_t i
    for i in range(length):
        if isnan(a_ptr[i]) != 0:
            return FiniteStatus.has_nan
        if isinf(a_ptr[i]) != 0:
            return FiniteStatus.has_infinite
    return FiniteStatus.all_finite
