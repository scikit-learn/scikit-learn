# cython: language_level=3

# -------------------------------------
# Readonly array wrapper
# -------------------------------------
# TODO: Remove with Cython >= 3.0 which supports const memoryviews for fused types.
#
# This class supports the buffer protocol, thus can wrap arrays and memoryviews.
# All it does is LIE about the readonly attribute: tell it's false!
# This way, we can use it on arrays that we don't touch.
# !!! USE CAREFULLY !!!


from cpython cimport Py_buffer
from cpython.buffer cimport PyObject_GetBuffer, PyBuffer_Release, PyBUF_WRITABLE

import numpy as np
cimport numpy as np


np.import_array()


ctypedef fused NUM_TYPES:
    np.npy_float64
    np.npy_float32
    np.npy_int64
    np.npy_int32


cdef class ReadonlyWrapper:
    cdef object wraps

    def __init__(self, wraps):
        self.wraps = wraps

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        request_for_writeable = False
        if flags & PyBUF_WRITABLE:
            flags ^= PyBUF_WRITABLE
            request_for_writeable = True
        PyObject_GetBuffer(self.wraps, buffer, flags)
        if request_for_writeable:
            buffer.readonly = False  # This is a lie!

    def __releasebuffer__(self, Py_buffer *buffer):
        PyBuffer_Release(buffer)


def _test_sum(NUM_TYPES[:] x):
    """This function is for testing only.

    As this function does not modify x, we would like to define it as
        _test_sum(const NUM_TYPES[:] x)
    which is not supported for fused types in Cython<3.0.
    """
    cdef:
        int i
        int n = x.shape[0]
        NUM_TYPES sum = 0

    for i in range(n):
        sum += x[i]
    return sum
