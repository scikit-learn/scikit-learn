# Author: jakirkham

cimport libc
cimport libc.math
from libc.math cimport isfinite as _isfinite, isinf as _isinf

cimport cython
from cython.parallel import prange
cimport numpy
import numpy
numpy.import_array()


cdef fused fprecision:
    float
    double

cdef fused const_fprecision:
    const float
    const double


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
    cdef numpy.ndarray a_flat

    cdef char* a_data
    cdef numpy.NPY_TYPES a_type
    cdef Py_ssize_t a_step
    cdef Py_ssize_t a_size

    cdef bint disallow_nan

    cdef bint result
    cdef bint_enum err


    a_flat = numpy.PyArray_Reshape(a, -1)

    with nogil:
        a_data = a_flat.data
        a_type = <numpy.NPY_TYPES>a_flat.descr.type_num
        a_step = a_flat.strides[0] / a_flat.descr.itemsize
        a_size = a_flat.shape[0]

        disallow_nan = not allow_nan

        err = bint_enum.false

        if a_type == numpy.NPY_TYPES.NPY_BOOL:
            result = True
        elif a_type == numpy.NPY_TYPES.NPY_BYTE:
            result = True
        elif a_type == numpy.NPY_TYPES.NPY_UBYTE:
            result = True
        elif a_type == numpy.NPY_TYPES.NPY_SHORT:
            result = True
        elif a_type == numpy.NPY_TYPES.NPY_USHORT:
            result = True
        elif a_type == numpy.NPY_TYPES.NPY_INT:
            result = True
        elif a_type == numpy.NPY_TYPES.NPY_UINT:
            result = True
        elif a_type == numpy.NPY_TYPES.NPY_LONG:
            result = True
        elif a_type == numpy.NPY_TYPES.NPY_ULONG:
            result = True
        elif a_type == numpy.NPY_TYPES.NPY_LONGLONG:
            result = True
        elif a_type == numpy.NPY_TYPES.NPY_ULONGLONG:
            result = True
        elif a_type == numpy.NPY_TYPES.NPY_FLOAT:
            result = c_isfinite(<const float*>a_data, a_step, a_size, <bint_enum>disallow_nan)
        elif a_type == numpy.NPY_TYPES.NPY_DOUBLE:
            result = c_isfinite(<const double*>a_data, a_step, a_size, <bint_enum>disallow_nan)
        else:
            err = bint_enum.true

    if err == bint_enum.false:
        return result
    elif err == bint_enum.true:
        raise TypeError("Unsupported array type: %s" % repr(numpy.PyArray_TypeObjectFromType(a_type)))


cdef inline bint c_isfinite(const_fprecision* a_ptr, Py_ssize_t step, Py_ssize_t size, bint_enum disallow_nan) nogil:
    if disallow_nan == bint_enum.true:
        return c_isfinite_bint_type(a_ptr, step, size, <bint_true_type*>NULL)
    elif disallow_nan == bint_enum.false:
        return c_isfinite_bint_type(a_ptr, step, size, <bint_false_type*>NULL)


cdef bint c_isfinite_bint_type(const_fprecision* a_ptr, Py_ssize_t step, Py_ssize_t size, bint_type* disallow_nan) nogil:
    cdef Py_ssize_t i

    for i in prange(0, size, step, nogil=True):
        if c_isnonfinite(a_ptr[i], disallow_nan):
            return False

    return True

cdef inline bint c_isnonfinite(fprecision v, bint_type* disallow_nan) nogil:

    if bint_type is bint_true_type:
        return _isfinite(v) == 0
    else:
        return _isinf(v) != 0
