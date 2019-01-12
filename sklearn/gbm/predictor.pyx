# cython: profile=True
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
"""
This module contains the TreePredictor class which is used for prediction.
"""
import numpy as np
cimport numpy as np

from .types import X_DTYPE


PREDICTOR_RECORD_DTYPE = np.dtype([
    ('value', np.float32),
    ('count', np.uint32),
    ('feature_idx', np.uint32),
    ('threshold', X_DTYPE),
    ('left', np.uint32),
    ('right', np.uint32),
    ('gain', np.float32),
    ('depth', np.uint32),
    ('is_leaf', np.uint8),
    ('bin_threshold', np.uint8),
    # TODO: shrinkage in leaf for feature importance error bar?
])

ctypedef np.npy_float64 NPY_X_DTYPE

cdef packed struct node_struct:
    float value
    unsigned int count
    unsigned int feature_idx
    NPY_X_DTYPE threshold
    unsigned int left
    unsigned int right
    float gain
    unsigned int depth
    unsigned char is_leaf
    unsigned char bin_threshold


class TreePredictor:
    """Tree class used for predictions.

    Parameters
    ----------
    nodes : list of PREDICTOR_RECORD_DTYPE.
        The nodes of the tree.
    """
    def __init__(self, nodes):
        self.nodes = nodes

    def get_n_leaf_nodes(self):
        """Return number of leaves."""
        return int(self.nodes['is_leaf'].sum())

    def get_max_depth(self):
        """Return maximum depth among all leaves."""
        return int(self.nodes['depth'].max())

    def predict(self, X):
        """Predict raw values for non-binned data.

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_features)
            The input samples.

        Returns
        -------
        y : array, shape (n_samples,)
            The raw predicted values.
        """
        # TODO: change dtype of out (should be same as Y_DTYPE I think since
        # the value is grad/hess which are Y_DTYPE)
        out = np.empty(X.shape[0], dtype=np.float32)
        _predict_from_numeric_data(self.nodes, X, out)
        return out


cdef float _predict_one_from_numeric_data(
    node_struct [:] nodes,
    const NPY_X_DTYPE [:] numeric_data) nogil:

    cdef:
        node_struct node = nodes[0]

    while True:
        if node.is_leaf:
            return node.value
        if numeric_data[node.feature_idx] <= node.threshold:
            node = nodes[node.left]
        else:
            node = nodes[node.right]


cdef void _predict_from_numeric_data(
    node_struct [:] nodes,
    const NPY_X_DTYPE [:, :] numeric_data,
    float [:] out) nogil:

    cdef:
        unsigned int i

    for i in range(numeric_data.shape[0]):
        out[i] = _predict_one_from_numeric_data(nodes, numeric_data[i])
