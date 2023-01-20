from libcpp.vector cimport vector

from ._tree cimport DTYPE_t          # Type of X
from ._tree cimport DOUBLE_t         # Type of y, sample_weight
from ._tree cimport SIZE_t           # Type for indices and counters
from ._tree cimport INT32_t          # Signed 32 bit integer
from ._tree cimport UINT32_t         # Unsigned 32 bit integer
from ._tree cimport Tree

cdef class QuantileTree(Tree):
    cdef vector[SIZE_t] leaf_node_ids   # the ids of each leaf node
    cdef vector[vector[DOUBLE_t]] leaf_samples  # (node_count, start:end) array of 'y' samples that fall into the leaf node
    

