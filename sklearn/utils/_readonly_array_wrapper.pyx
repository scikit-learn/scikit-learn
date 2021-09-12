"""
ReadonlyArrayWrapper implements the buffer protocol to make the wraped buffer behave as if
writeable, even for readonly buffers. This way, even readonly arrays can be passed as
argument of type (non const) memoryview.
This is a workaround for the missing support for const fused-typed memoryviews in
Cython < 3.0.

Note: All it does is LIE about the readonly attribute: tell it's false!
This way, we can use it on arrays that we don't touch.
!!! USE CAREFULLY !!!
"""
# TODO: Remove with Cython >= 3.0 which supports const memoryviews for fused types.

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


cdef class ReadonlyArrayWrapper:
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
            # The following is a lie when self.wraps is readonly!
            buffer.readonly = False

    def __releasebuffer__(self, Py_buffer *buffer):
        PyBuffer_Release(buffer)


def _test_sum(NUM_TYPES[:] x):
    """This function is for testing only.

    As this function does not modify x, we would like to define it as

            _test_sum(const NUM_TYPES[:] x)

    which is not possible as fused typed const memoryviews aren't
    supported in Cython<3.0.
    """
    cdef:
        int i
        int n = x.shape[0]
        NUM_TYPES sum = 0

    for i in range(n):
        sum += x[i]
    return sum
