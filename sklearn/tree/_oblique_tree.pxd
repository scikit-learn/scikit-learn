# distutils: language = c++

# Authors: Gilles Louppe <g.louppe@gmail.com>
#          Peter Prettenhofer <peter.prettenhofer@gmail.com>
#          Brian Holt <bdholt1@gmail.com>
#          Joel Nothman <joel.nothman@gmail.com>
#          Arnaud Joly <arnaud.v.joly@gmail.com>
#          Jacob Schreiber <jmschreiber91@gmail.com>
#          Nelson Liu <nelson@nelsonliu.me>
#
# License: BSD 3 clause

# See _oblique_tree.pyx for details.

import numpy as np
cimport numpy as np

from libcpp.vector cimport vector

from ._tree cimport DTYPE_t          # Type of X
from ._tree cimport DOUBLE_t         # Type of y, sample_weight
from ._tree cimport SIZE_t           # Type for indices and counters
from ._tree cimport INT32_t          # Signed 32 bit integer
from ._tree cimport UINT32_t         # Unsigned 32 bit integer
from ._tree cimport Tree, Node
from ._oblique_splitter cimport ObliqueSplitter
from ._oblique_splitter cimport ObliqueSplitRecord


cdef class ObliqueTree(Tree):
    # The Tree object is a binary tree structure constructed by the
    # TreeBuilder. The tree structure is used for predictions and
    # feature importances.

    # oblique forests
    # cdef Node* oblique_nodes               # Array of "oblique" nodes
    cdef vector[vector[DTYPE_t]] proj_vec_weights # (capacity, n_features) array of projection vectors
    cdef vector[vector[SIZE_t]] proj_vec_indices  # (capacity, n_features) array of projection vectors

    # Methods
    cdef SIZE_t _add_oblique_node(self, SIZE_t parent, bint is_left, bint is_leaf,
                          SIZE_t feature, double threshold, double impurity,
                          SIZE_t n_node_samples,
                          double weighted_n_samples,
                          vector[DTYPE_t]* proj_vec_weights,
                          vector[SIZE_t]* proj_vec_indices) nogil except -1


# =============================================================================
# Tree builder
# =============================================================================

cdef class ObliqueTreeBuilder:
    # The TreeBuilder recursively builds a Tree object from training samples,
    # using a Splitter object for splitting internal nodes and assigning
    # values to leaves.
    #
    # This class controls the various stopping criteria and the node splitting
    # evaluation order, e.g. depth-first or best-first.

    cdef ObliqueSplitter splitter   # Splitting algorithm

    cdef SIZE_t min_samples_split       # Minimum number of samples in an internal node
    cdef SIZE_t min_samples_leaf        # Minimum number of samples in a leaf
    cdef double min_weight_leaf         # Minimum weight in a leaf
    cdef SIZE_t max_depth               # Maximal tree depth
    cdef double min_impurity_split
    cdef double min_impurity_decrease   # Impurity threshold for early stopping

    cpdef build(self, ObliqueTree tree, object X, np.ndarray y,
                np.ndarray sample_weight=*,
                np.ndarray X_idx_sorted=*)
    cdef _check_input(self, object X, np.ndarray y, np.ndarray sample_weight)
