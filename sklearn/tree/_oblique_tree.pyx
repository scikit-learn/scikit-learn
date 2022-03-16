# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False

# Authors: Adam Li <adam2392@gmail.com>
#          Chester Huynh <chester.huynh924@gmail.com>
#          Parth Vora <pvora4@jhu.edu>
#
# License: BSD 3 clause

from cpython cimport Py_INCREF, PyObject, PyTypeObject

from libc.stdlib cimport malloc, free
from libc.math cimport fabs
from libc.string cimport memcpy
from libc.string cimport memset
from libc.stdint cimport SIZE_MAX
from libc.math cimport round

import numpy as np
cimport numpy as np
np.import_array()

from scipy.sparse import issparse
from scipy.sparse import csr_matrix

from cython.operator cimport dereference as deref

from ._utils cimport safe_realloc
from ._utils cimport sizet_ptr_to_ndarray

# Gets Node dtype exposed inside oblique_tree.
# See "_tree.pyx" for more details.
cdef Node dummy;
NODE_DTYPE = np.asarray(<Node[:1]>(&dummy)).dtype

# Mitigate precision differences between 32 bit and 64 bit
cdef DTYPE_t FEATURE_THRESHOLD = 1e-7

# =============================================================================
# ObliqueTree
# =============================================================================

cdef class ObliqueTree(Tree):
    """Array-based representation of a binary oblique decision tree.

    The oblique decision tree is represented as a number of parallel arrays. The i-th
    element of each array holds information about the node `i`. Node 0 is the
    tree's root. You can find a detailed description of all arrays in
    `_tree.pxd`. NOTE: Some of the arrays only apply to either leaves or split
    nodes, resp. In this case the values of nodes of the other type are
    arbitrary!

    Attributes
    ----------
    node_count : int
        The number of nodes (internal nodes + leaves) in the tree.

    capacity : int
        The current capacity (i.e., size) of the arrays, which is at least as
        great as `node_count`.

    max_depth : int
        The depth of the tree, i.e. the maximum depth of its leaves.

    children_left : array of int, shape [node_count]
        children_left[i] holds the node id of the left child of node i.
        For leaves, children_left[i] == TREE_LEAF. Otherwise,
        children_left[i] > i. This child handles the case where
        X[:, feature[i]] <= threshold[i].

    children_right : array of int, shape [node_count]
        children_right[i] holds the node id of the right child of node i.
        For leaves, children_right[i] == TREE_LEAF. Otherwise,
        children_right[i] > i. This child handles the case where
        X[:, feature[i]] > threshold[i].

    feature : array of int, shape [node_count]
        feature[i] holds the feature to split on, for the internal node i.

    threshold : array of double, shape [node_count]
        threshold[i] holds the threshold for the internal node i.

    value : array of double, shape [node_count, n_outputs, max_n_classes]
        Contains the constant prediction value of each node.

    impurity : array of double, shape [node_count]
        impurity[i] holds the impurity (i.e., the value of the splitting
        criterion) at node i.

    n_node_samples : array of int, shape [node_count]
        n_node_samples[i] holds the number of training samples reaching node i.

    weighted_n_node_samples : array of int, shape [node_count]
        weighted_n_node_samples[i] holds the weighted number of training samples
        reaching node i.
    """
    def __cinit__(self, int n_features, np.ndarray[SIZE_t, ndim=1] n_classes,
                int n_outputs):
        """Constructor."""
        # Input/Output layout
        self.n_features = n_features
        self.n_outputs = n_outputs
        self.n_classes = NULL
        safe_realloc(&self.n_classes, n_outputs)

        self.max_n_classes = np.max(n_classes)
        self.value_stride = n_outputs * self.max_n_classes

        cdef SIZE_t k
        for k in range(n_outputs):
            self.n_classes[k] = n_classes[k]

        # Inner structures
        self.max_depth = 0
        self.node_count = 0
        self.capacity = 0
        self.value = NULL
        self.nodes = NULL

        self.proj_vec_weights = vector[vector[DTYPE_t]](self.capacity)
        self.proj_vec_indices = vector[vector[SIZE_t]](self.capacity)

    def __reduce__(self):
        """Reduce re-implementation, for pickling."""
        return (ObliqueTree, (self.n_features,
                       sizet_ptr_to_ndarray(self.n_classes, self.n_outputs),
                       self.n_outputs), self.__getstate__())

    def __getstate__(self):
        """Getstate re-implementation, for pickling."""
        d = {}
        # capacity is inferred during the __setstate__ using nodes
        d["max_depth"] = self.max_depth
        d["node_count"] = self.node_count
        d["nodes"] = self._get_node_ndarray()
        d["values"] = self._get_value_ndarray()

        proj_vecs = self.get_projection_matrix()
        d['proj_vecs'] = proj_vecs
        return d

    def __setstate__(self, d):
        """Setstate re-implementation, for unpickling."""
        self.max_depth = d["max_depth"]
        self.node_count = d["node_count"]

        if 'nodes' not in d:
            raise ValueError('You have loaded ObliqueTree version which '
                            'cannot be imported')

        node_ndarray = d['nodes']
        value_ndarray = d['values']

        value_shape = (node_ndarray.shape[0], self.n_outputs,
                    self.max_n_classes)
        if (node_ndarray.ndim != 1 or
                node_ndarray.dtype != NODE_DTYPE or
                not node_ndarray.flags.c_contiguous or
                value_ndarray.shape != value_shape or
                not value_ndarray.flags.c_contiguous or
                value_ndarray.dtype != np.float64):
            raise ValueError('Did not recognise loaded array layout')

        self.capacity = node_ndarray.shape[0]
        if self._resize_c(self.capacity) != 0:
            raise MemoryError("resizing tree to %d" % self.capacity)

        # now set the projection vector weights and indices
        proj_vecs = d['proj_vecs']
        self.n_features = proj_vecs.shape[1]
        for i in range(0, self.node_count):
            for j in range(0, self.n_features):
                weight = proj_vecs[i, j]
                if weight == 0:
                    continue
                self.proj_vec_weights[i].push_back(weight)
                self.proj_vec_indices[i].push_back(j)

        nodes = memcpy(self.nodes, (<np.ndarray> node_ndarray).data,
                    self.capacity * sizeof(Node))
        value = memcpy(self.value, (<np.ndarray> value_ndarray).data,
                    self.capacity * self.value_stride * sizeof(double))

    cpdef np.ndarray get_projection_matrix(self):
        """Get the projection matrix of shape (node_count, n_features)."""
        proj_vecs = np.zeros((self.node_count, self.n_features))
        for i in range(0, self.node_count):
            for j in range(0, self.proj_vec_weights[i].size()):
                weight = self.proj_vec_weights[i][j]
                feat = self.proj_vec_indices[i][j]
                proj_vecs[i, feat] = weight
        return proj_vecs

    cdef int _resize_c(self, SIZE_t capacity=SIZE_MAX) nogil except -1:
        """Guts of _resize.

        Additionally resizes the projection indices and weights.

        Returns -1 in case of failure to allocate memory (and raise MemoryError)
        or 0 otherwise.
        """
        if capacity == self.capacity and self.nodes != NULL:
            return 0

        if capacity == SIZE_MAX:
            if self.capacity == 0:
                capacity = 3  # default initial value
            else:
                capacity = 2 * self.capacity

        safe_realloc(&self.nodes, capacity)
        safe_realloc(&self.value, capacity * self.value_stride)

        # only thing added for oblique trees
        # TODO: this could possibly be removed if we can 
        # add projection indices and weights to Node
        self.proj_vec_weights.resize(capacity)
        self.proj_vec_indices.resize(capacity)

        # value memory is initialised to 0 to enable classifier argmax
        if capacity > self.capacity:
            memset(<void*>(self.value + self.capacity * self.value_stride), 0,
                (capacity - self.capacity) * self.value_stride *
                sizeof(double))

        # if capacity smaller than node_count, adjust the counter
        if capacity < self.node_count:
            self.node_count = capacity

        self.capacity = capacity
        return 0

    cdef int _set_node_values(self, SplitRecord* split_node, Node *node) nogil except -1:
        """Set node data.
        """
        cdef ObliqueSplitRecord* oblique_split_node = <ObliqueSplitRecord*>(split_node)
        cdef SIZE_t node_id = self.node_count

        node.feature = deref(oblique_split_node).feature
        node.threshold = deref(oblique_split_node).threshold

        # oblique trees store the projection indices and weights
        # inside the tree itself
        self.proj_vec_weights[node_id] = deref(deref(oblique_split_node).proj_vec_weights)
        self.proj_vec_indices[node_id] = deref(deref(oblique_split_node).proj_vec_indices)
        return 1

    cdef DTYPE_t _compute_feature(self, const DTYPE_t[:] X_ndarray, Node *node, SIZE_t node_id) nogil:
        """Compute feature from a given data matrix, X.

        In oblique-aligned trees, this is the projection of X.
        In this case, we take a simple linear combination of some columns of X.
        """
        cdef vector[DTYPE_t] proj_vec_weights
        cdef vector[SIZE_t] proj_vec_indices
        cdef DTYPE_t proj_feat = 0.0

        # compute projection of the data based on trained tree
        proj_vec_weights = self.proj_vec_weights[node_id]
        proj_vec_indices = self.proj_vec_indices[node_id]
        for j in range(proj_vec_indices.size()):
            # proj_feat += round(X_ndarray[proj_vec_indices[j]] * proj_vec_weights[j] * 1.0 / FEATURE_THRESHOLD)  * FEATURE_THRESHOLD
            proj_feat += X_ndarray[proj_vec_indices[j]] * proj_vec_weights[j]

        # round to handle numerical precision
        # proj_feat = round(proj_feat * 1.0 / FEATURE_THRESHOLD) * FEATURE_THRESHOLD

        return proj_feat