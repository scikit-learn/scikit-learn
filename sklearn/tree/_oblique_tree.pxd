# distutils: language = c++

# Authors: Adam Li <adam2392@gmail.com>
#          Chester Huynh <chester.huynh924@gmail.com>
#          Parth Vora <pvora4@jhu.edu>
#
# License: BSD 3 clause

# See _oblique_tree.pyx for details.

import numpy as np
cimport numpy as cnp

from libcpp.vector cimport vector

from ._tree cimport DTYPE_t          # Type of X
from ._tree cimport DOUBLE_t         # Type of y, sample_weight
from ._tree cimport SIZE_t           # Type for indices and counters
from ._tree cimport INT32_t          # Signed 32 bit integer
from ._tree cimport UINT32_t         # Unsigned 32 bit integer
from ._tree cimport Tree, Node, TreeBuilder

from ._splitter cimport SplitRecord
from ._oblique_splitter cimport ObliqueSplitRecord

cdef class ObliqueTree(Tree):
    cdef vector[vector[DTYPE_t]] proj_vec_weights # (capacity, n_features) array of projection vectors
    cdef vector[vector[SIZE_t]] proj_vec_indices  # (capacity, n_features) array of projection vectors

    cdef int _resize_c(self, SIZE_t capacity=*) nogil except -1
    cdef int _set_node_values(self, SplitRecord* split_node, Node *node)  nogil except -1
    cdef DTYPE_t _compute_feature(self, const DTYPE_t[:] X_ndarray, Node *node, SIZE_t node_id) nogil
    cdef void _compute_feature_importances(self, DOUBLE_t* importance_data,
                                Node* node, SIZE_t node_id) nogil

    cpdef DTYPE_t compute_feature_value(self, object X, SIZE_t node_id)
    cpdef cnp.ndarray get_projection_matrix(self)
