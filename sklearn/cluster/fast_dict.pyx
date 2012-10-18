"""
Uses C++ map containers for fast dict-like behavior with keys being
integers, and values float.
"""
# Author: Gael Varoquaux
# License: BSD

cimport cython

# C++
from cython.operator cimport dereference as deref, preincrement as inc
from libcpp.map cimport map as cpp_map
from libc.math cimport fmin

import numpy as np

# Import the C-level symbols of numpy
cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int32
ctypedef np.int32_t ITYPE_t

###############################################################################
# An object to be used in Python

# Lookup is faster than dict (up to 10 times), and so is full traversal
# (up to 50 times), and assignement (up to 6 times), but creation is
# slower (up to 3 times)

cdef class IntFloatDict:
    cdef cpp_map[ITYPE_t, DTYPE_t] my_map

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __init__(self, np.ndarray[ITYPE_t, ndim=1] keys,
                       np.ndarray[DTYPE_t, ndim=1] values):
        cdef cpp_map[ITYPE_t,DTYPE_t] my_map = self.my_map
        cdef int i
        cdef int size = values.size
        # Should check that sizes for keys and values are equal, and
        # after should boundcheck(False)
        for i in range(size):
            my_map[keys[i]] = values[i]
        self.my_map = my_map

    def __len__(self):
        return self.my_map.size()

    def __getitem__(self, int key):
        cdef cpp_map[ITYPE_t, DTYPE_t].iterator it = self.my_map.find(key)
        if it == self.my_map.end():
            # The key is not in the dict
            raise KeyError('%i' % key)
        return deref(it).second

    def __setitem__(self, int key, float value):
        self.my_map[key] = value

    def __iter__(self):
        cdef cpp_map[ITYPE_t, DTYPE_t].iterator it = self.my_map.begin()
        cdef cpp_map[ITYPE_t, DTYPE_t].iterator end = self.my_map.end()
        while it != end:
            yield deref(it).first, deref(it).second
            inc(it)

    def update(self, IntFloatDict other):
        cdef cpp_map[ITYPE_t, DTYPE_t].iterator it = other.my_map.begin()
        cdef cpp_map[ITYPE_t, DTYPE_t].iterator end = other.my_map.end()
        while it != end:
            self.my_map[deref(it).first] = deref(it).second
            inc(it)

    def copy(self):
        cdef IntFloatDict out_obj = IntFloatDict.__new__(IntFloatDict)
        # The '=' operator is a copy operator for C++ maps
        out_obj.my_map = self.my_map
        return out_obj


###############################################################################
# An object with a specific merge strategy

@cython.boundscheck(False)
@cython.wraparound(False)
def min_merge(IntFloatDict a, IntFloatDict b,
              np.ndarray[ITYPE_t, ndim=1] mask):
    a = a.copy()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator my_it = a.my_map.begin()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator my_end = a.my_map.end()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator other_it = b.my_map.begin()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator other_end = b.my_map.end()
    cdef ITYPE_t other_key
    cdef DTYPE_t other_value
    while other_it != other_end:
        other_key = deref(other_it).first
        other_value = deref(other_it).second
        if mask[other_key]:
            my_it = a.my_map.find(other_key)
            if my_it == my_end:
                # Key not found
                a.my_map[other_key] = other_value
            else:
                deref(my_it).second = fmin(deref(my_it).second, other_value)
        inc(other_it)
    return a


