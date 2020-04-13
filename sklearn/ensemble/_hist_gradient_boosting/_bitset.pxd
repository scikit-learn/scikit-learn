# cython: language_level=3
from .common cimport X_BITSET_DTYPE_C
from .common cimport X_BINNED_DTYPE_C


cdef void init_bitset(X_BITSET_DTYPE_C bitset) nogil

cdef void insert_bitset(X_BINNED_DTYPE_C val, X_BITSET_DTYPE_C bitset) nogil

cdef unsigned char in_bitset(X_BINNED_DTYPE_C val, X_BITSET_DTYPE_C bitset) nogil
