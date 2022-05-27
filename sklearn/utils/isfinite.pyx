# cython: profile=True, linetrace=True
# Author: jakirkham

cimport libc
cimport libc.math
from libc.math cimport isnan as _isnan, isinf as _isinf

cimport cython
from cython.parallel import prange
cimport numpy
import numpy
numpy.import_array()


cdef fused fprecision:
    numpy.npy_float32
    numpy.npy_float64

cdef fused const_fprecision:
    const numpy.npy_float32
    const numpy.npy_float64

@cython.cdivision(True)
def cy_isfinite(numpy.ndarray a, bint allow_nan=False):
    cdef char* a_data
    cdef numpy.NPY_TYPES a_type
    cdef Py_ssize_t a_size = a.size
    cdef bint disallow_nan

    cdef int result
    cdef bint err

    with nogil:
        a_data = a.data
        a_type = <numpy.NPY_TYPES>a.descr.type_num

        disallow_nan = not allow_nan

        err = False
        if a_type == numpy.NPY_TYPES.NPY_FLOAT:
            result = c_isfinite(<const float*>a_data, a_size, disallow_nan)
        elif a_type == numpy.NPY_TYPES.NPY_DOUBLE:
            result = c_isfinite(<const double*>a_data, a_size, disallow_nan)
        else:
            err = True

    if err == False:
        return result
    elif err == True:
        raise TypeError("Unsupported array type: %s" % repr(numpy.PyArray_TypeObjectFromType(a_type)))


cdef inline int c_isfinite(const_fprecision* a_ptr, Py_ssize_t size, bint disallow_nan) nogil:
    if disallow_nan:
        return c_isfinite_disable_nan(a_ptr, size)
    else:
        return c_isfinite_allow_nan(a_ptr, size)

cdef enum FiniteStatus:
    all_finite = 0
    has_nan = 1
    has_infinite = 2

cdef int c_isfinite_allow_nan(const_fprecision* a_ptr, Py_ssize_t size) nogil:
    for i in range(size):
        if _isinf(a_ptr[i]) != 0:
            return FiniteStatus.has_infinite
    return FiniteStatus.all_finite

cdef int c_isfinite_disable_nan(const_fprecision* a_ptr, Py_ssize_t size) nogil:
    for i in range(size):
        if _isnan(a_ptr[i]) != 0:
            return FiniteStatus.has_nan
        if _isinf(a_ptr[i]) != 0:
            return FiniteStatus.has_infinite
    return FiniteStatus.all_finite
