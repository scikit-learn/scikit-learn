cimport cython
cimport numpy as np
import numpy as np


cdef class ArrayBuilder(object):
    """Helper class to incrementally build a 1-d numpy.ndarray"""
    # Or: let's reinvent the wheel!

    GROWTH_FACTOR = 1.5

    cdef Py_ssize_t _nelems
    cdef object _arr    # object because we don't know the dtype statically

    def __init__(self, dtype, initial_capacity=256):
        assert self.GROWTH_FACTOR > 1
        assert initial_capacity >= 2
        self._arr = np.empty(initial_capacity, dtype=dtype)
        self._nelems = 0

    def __len__(self):
        return self._nelems

    @cython.boundscheck(False)
    def append(self, x):
        """Append a single value.

        Complexity: amortized O(1).
        """
        if self._nelems == self._arr.size:
            self._grow()
        self._arr[self._nelems] = x
        self._nelems += 1

    def get(self):
        """Return the constructed array.

        Don't use an ArrayBuilder after calling this method.
        """
        self._arr.resize(self._nelems)
        return self._arr

    cdef _grow(self):
        self._arr.resize(self._arr.size * self.GROWTH_FACTOR)
