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


cdef enum bint_enum:
    false = False
    true = True


cdef struct bint_false_type:
    char empty

cdef struct bint_true_type:
    char empty

cdef fused bint_type:
    bint_false_type
    bint_true_type


@cython.cdivision(True)
def cy_isfinite(numpy.ndarray a, bint allow_nan=False):
    cdef char* a_data
    cdef numpy.NPY_TYPES a_type
    cdef Py_ssize_t a_size = a.size
    cdef bint disallow_nan

    cdef int result
    cdef bint_enum err

    with nogil:
        a_data = a.data
        a_type = <numpy.NPY_TYPES>a.descr.type_num

        disallow_nan = not allow_nan

        err = bint_enum.false
        if a_type == numpy.NPY_TYPES.NPY_FLOAT:
            result = c_isfinite(<const float*>a_data, a_size, <bint_enum>disallow_nan)
        elif a_type == numpy.NPY_TYPES.NPY_DOUBLE:
            result = c_isfinite(<const double*>a_data, a_size, <bint_enum>disallow_nan)
        else:
            err = bint_enum.true

    if err == bint_enum.false:
        return result
    elif err == bint_enum.true:
        raise TypeError("Unsupported array type: %s" % repr(numpy.PyArray_TypeObjectFromType(a_type)))


cdef inline int c_isfinite(const_fprecision* a_ptr, Py_ssize_t size, bint_enum disallow_nan) nogil:
    if disallow_nan == bint_enum.true:
        return c_isfinite_bint_type(a_ptr, size, <bint_true_type*>NULL)
    elif disallow_nan == bint_enum.false:
        return c_isfinite_bint_type(a_ptr, size, <bint_false_type*>NULL)


cdef int c_isfinite_bint_type(const_fprecision* a_ptr, Py_ssize_t size, bint_type* disallow_nan) nogil:
    cdef Py_ssize_t i
    cdef int out = 0

    for i from 0 <= i < size:
        out = c_isnonfinite(a_ptr[i], disallow_nan)
        if out !=0:
            return out
    return out

cdef inline int c_isnonfinite(fprecision v, bint_type* disallow_nan) nogil:
    # Returns 0 if finite
    # Returns 1 if NaN
    # Returns 2 if inf

    if bint_type is bint_true_type and _isnan(v) == 1:
        return 1
    return _isinf(v)*2
