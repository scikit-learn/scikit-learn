# Author: jakirkham, Meekail Zain, Thomas Fan

from libc.math cimport isnan as c_isnan, isinf as c_isinf

cimport cython
cimport numpy as cnp
import numpy as np
cnp.import_array()


cdef fused fprecision:
    cnp.npy_float32
    cnp.npy_float64

cdef fused const_fprecision:
    const cnp.npy_float32
    const cnp.npy_float64

def cy_isfinite(cnp.ndarray a, bint allow_nan=False):
    cdef:
        char* a_data
        cnp.NPY_TYPES a_type
        Py_ssize_t a_size = a.size
        int result
        bint err

    if a.dtype.type not in {np.float32, np.float64}:
        raise TypeError("Unsupported array type: %s" % repr(a.dtype))

    with nogil:
        a_data = a.data
        a_type = <cnp.NPY_TYPES>a.descr.type_num

        if a_type == cnp.NPY_TYPES.NPY_FLOAT:
            result = _isfinite(<const float*>a_data, a_size, allow_nan)
        elif a_type == cnp.NPY_TYPES.NPY_DOUBLE:
            result = _isfinite(<const double*>a_data, a_size, allow_nan)

    return result

cdef inline int _isfinite(const_fprecision* a_ptr, Py_ssize_t size, bint allow_nan) nogil:
    if allow_nan:
        return _isfinite_allow_nan(a_ptr, size)
    else:
        return _isfinite_disable_nan(a_ptr, size)

cpdef enum FiniteStatus:
    all_finite = 0
    has_nan = 1
    has_infinite = 2

cdef int _isfinite_allow_nan(const_fprecision* a_ptr, Py_ssize_t size) nogil:
    for i in range(size):
        if c_isinf(a_ptr[i]) != 0:
            return FiniteStatus.has_infinite
    return FiniteStatus.all_finite

cdef int _isfinite_disable_nan(const_fprecision* a_ptr, Py_ssize_t size) nogil:
    for i in range(size):
        if c_isnan(a_ptr[i]) != 0:
            return FiniteStatus.has_nan
        if c_isinf(a_ptr[i]) != 0:
            return FiniteStatus.has_infinite
    return FiniteStatus.all_finite
