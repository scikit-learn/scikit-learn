# Author: jakirkham

cimport libc
cimport libc.math
from libc.math cimport isfinite as _isfinite, isinf as _isinf
from libc.stdint cimport uint16_t

cimport cython

cimport numpy
import numpy
numpy.import_array()


cdef extern from "numpy/halffloat.h":
    ctypedef uint16_t npy_half

    # conversion function
    float npy_half_to_float(npy_half h);


cdef fused fprecision:
    npy_half
    float
    double
    float complex
    double complex

cdef fused const_fprecision:
    const npy_half
    const float
    const double
    const float complex
    const double complex


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
def py_isfinite(numpy.ndarray a, bint allow_nan=False):
    cdef numpy.ndarray a_flat

    cdef char* a_data
    cdef numpy.NPY_TYPES a_type
    cdef size_t a_step
    cdef size_t a_size

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

        if a_type == numpy.NPY_BOOL:
            result = True
        elif a_type == numpy.NPY_BYTE:
            result = True
        elif a_type == numpy.NPY_UBYTE:
            result = True
        elif a_type == numpy.NPY_SHORT:
            result = True
        elif a_type == numpy.NPY_USHORT:
            result = True
        elif a_type == numpy.NPY_INT:
            result = True
        elif a_type == numpy.NPY_UINT:
            result = True
        elif a_type == numpy.NPY_LONG:
            result = True
        elif a_type == numpy.NPY_ULONG:
            result = True
        elif a_type == numpy.NPY_LONGLONG:
            result = True
        elif a_type == numpy.NPY_ULONGLONG:
            result = True
        elif a_type == numpy.NPY_HALF:
            result = c_isfinite(<const npy_half*>a_data, a_step, a_size, <bint_enum>disallow_nan)
        elif a_type == numpy.NPY_FLOAT:
            result = c_isfinite(<const float*>a_data, a_step, a_size, <bint_enum>disallow_nan)
        elif a_type == numpy.NPY_DOUBLE:
            result = c_isfinite(<const double*>a_data, a_step, a_size, <bint_enum>disallow_nan)
        elif a_type == numpy.NPY_CFLOAT:
            result = c_isfinite(<const float complex*>a_data, a_step, a_size, <bint_enum>disallow_nan)
        elif a_type == numpy.NPY_CDOUBLE:
            result = c_isfinite(<const double complex*>a_data, a_step, a_size, <bint_enum>disallow_nan)
        else:
            err = bint_enum.true

    if err == bint_enum.false:
        return result
    elif err == bint_enum.true:
        raise TypeError("Unsupported array type: %s" % repr(numpy.PyArray_TypeObjectFromType(a_type)))


cdef inline bint c_isfinite(const_fprecision* a_ptr, size_t step, size_t size, bint_enum disallow_nan) nogil:
    if disallow_nan == bint_enum.true:
        return c_isfinite_bint_type(a_ptr, step, size, <bint_true_type*>NULL)
    elif disallow_nan == bint_enum.false:
        return c_isfinite_bint_type(a_ptr, step, size, <bint_false_type*>NULL)


cdef bint c_isfinite_bint_type(const_fprecision* a_ptr, size_t step, size_t size, bint_type* disallow_nan) nogil:
    cdef size_t i

    for i from 0 <= i < size by step:
        if c_isnonfinite(a_ptr[i], disallow_nan):
            return False

    return True


cdef inline bint c_isnonfinite(fprecision v, bint_type* disallow_nan) nogil:

    if fprecision is floatcomplex or fprecision is doublecomplex:
        return c_isnonfinite(v.real, disallow_nan) or c_isnonfinite(v.imag, disallow_nan)
    elif fprecision is npy_half:
        with gil:
            if bint_type is bint_true_type:
                return _isfinite(npy_half_to_float(v)) == 0
            else:
                return _isinf(npy_half_to_float(v)) != 0
    elif fprecision is float or fprecision is double:
        if bint_type is bint_true_type:
            return _isfinite(v) == 0
        else:
            return _isinf(v) != 0
