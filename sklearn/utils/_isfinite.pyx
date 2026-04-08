# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from libc.math cimport isnan, isinf
from cython cimport floating


cpdef enum FiniteStatus:
    all_finite = 0
    has_nan = 1
    has_infinite = 2


def cy_isfinite(floating[::1] a, bint allow_nan=False):
    cdef FiniteStatus result
    with nogil:
        result = _isfinite(a, allow_nan)
    return result


cdef inline FiniteStatus _isfinite(floating[::1] a, bint allow_nan) noexcept nogil:
    cdef floating* a_ptr = &a[0]
    cdef Py_ssize_t length = len(a)
    if allow_nan:
        return _isfinite_allow_nan(a_ptr, length)
    else:
        return _isfinite_disable_nan(a_ptr, length)


cdef inline FiniteStatus _isfinite_allow_nan(floating* a_ptr,
                                             Py_ssize_t length) noexcept nogil:
    cdef Py_ssize_t i
    cdef floating v
    for i in range(length):
        v = a_ptr[i]
        if isinf(v):
            return FiniteStatus.has_infinite
    return FiniteStatus.all_finite


cdef inline FiniteStatus _isfinite_disable_nan(floating* a_ptr,
                                               Py_ssize_t length) noexcept nogil:
    cdef Py_ssize_t i
    cdef floating v
    for i in range(length):
        v = a_ptr[i]
        if isnan(v):
            return FiniteStatus.has_nan
        elif isinf(v):
            return FiniteStatus.has_infinite
    return FiniteStatus.all_finite
