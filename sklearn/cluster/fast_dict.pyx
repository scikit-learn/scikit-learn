"""
Uses C++ map containers for fast dict-like behavior with keys being
integers, and values float.
"""
# Author: Gael Varoquaux
# License: BSD

cimport cython

# C++
from cython.operator cimport dereference as deref, preincrement as inc, \
    predecrement as dec
from libcpp.utility cimport pair
from libcpp.map cimport map as cpp_map
from libc.math cimport fmin, fmax

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

    def append(self, ITYPE_t key, DTYPE_t value):
        cdef cpp_map[ITYPE_t, DTYPE_t].iterator end = self.my_map.end()
        # Decrement the iterator
        dec(end)
        # Construct our arguments
        cdef pair[ITYPE_t, DTYPE_t] args
        args.first = key
        args.second = value
        self.my_map.insert(end, args)


###############################################################################
# An object with a specific merge strategy

@cython.boundscheck(False)
@cython.wraparound(False)
def min_merge(IntFloatDict a, IntFloatDict b,
              np.ndarray[ITYPE_t, ndim=1] mask):
    cdef IntFloatDict out_obj = IntFloatDict.__new__(IntFloatDict)
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator a_it = a.my_map.begin()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator a_end = a.my_map.end()
    cdef ITYPE_t key
    cdef DTYPE_t value
    # First copy a into out
    while a_it != a_end:
        key = deref(a_it).first
        if mask[key]:
            out_obj.my_map[key] = deref(a_it).second
        inc(a_it)

    # Then merge b into out
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator out_it = out_obj.my_map.begin()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator out_end = out_obj.my_map.end()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator b_it = b.my_map.begin()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator b_end = b.my_map.end()
    while b_it != b_end:
        key = deref(b_it).first
        value = deref(b_it).second
        if mask[key]:
            out_it = out_obj.my_map.find(key)
            if out_it == out_end:
                # Key not found
                out_obj.my_map[key] = value
            else:
                deref(out_it).second = fmin(deref(out_it).second, value)
        inc(b_it)
    return out_obj


@cython.boundscheck(False)
@cython.wraparound(False)
def max_merge(IntFloatDict a, IntFloatDict b,
              np.ndarray[ITYPE_t, ndim=1] mask):
    cdef IntFloatDict out_obj = IntFloatDict.__new__(IntFloatDict)
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator a_it = a.my_map.begin()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator a_end = a.my_map.end()
    cdef ITYPE_t key
    cdef DTYPE_t value
    # First copy a into out
    while a_it != a_end:
        key = deref(a_it).first
        if mask[key]:
            out_obj.my_map[key] = deref(a_it).second
        inc(a_it)

    # Then merge b into out
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator out_it = out_obj.my_map.begin()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator out_end = out_obj.my_map.end()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator b_it = b.my_map.begin()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator b_end = b.my_map.end()
    while b_it != b_end:
        key = deref(b_it).first
        value = deref(b_it).second
        if mask[key]:
            out_it = out_obj.my_map.find(key)
            if out_it == out_end:
                # Key not found
                out_obj.my_map[key] = value
            else:
                deref(out_it).second = fmax(deref(out_it).second, value)
        inc(b_it)
    return out_obj


###############################################################################
# An edge object for fast comparisons 

cdef class WeightedEdge:
    cdef public ITYPE_t a
    cdef public ITYPE_t b
    cdef public DTYPE_t weight
    
    def __init__(self, DTYPE_t weight, ITYPE_t a, ITYPE_t b):
        self.weight = weight
        self.a = a
        self.b = b

    @cython.nonecheck(False)
    def __cmp__(self, WeightedEdge other):
        """Return negative if self < other, zero if self == other,
           positive if self > other.
        """
        return self.weight > other.weight
        
    def __repr__(self):
        return "%s(weight=%f, a=%i, b=%i)" % (self.__class__.__name__,
                                              self.weight,
                                              self.a, self.b)

