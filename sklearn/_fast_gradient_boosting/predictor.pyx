# cython: profile=True
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
"""
This module contains the TreePredictor class which is used for prediction.
"""
cimport cython
from cython.parallel import prange
import numpy as np
cimport numpy as np

from .types import X_DTYPE
from .types cimport X_DTYPE_C
from .types import Y_DTYPE
from .types cimport Y_DTYPE_C
from .types import X_BINNED_DTYPE
from .types cimport X_BINNED_DTYPE_C


PREDICTOR_RECORD_DTYPE = np.dtype([
    ('value', Y_DTYPE),
    ('count', np.uint32),
    ('feature_idx', np.uint32),
    ('threshold', X_DTYPE),
    ('left', np.uint32),
    ('right', np.uint32),
    ('gain', Y_DTYPE),
    ('depth', np.uint32),
    ('is_leaf', np.uint8),
    ('bin_threshold', X_BINNED_DTYPE),
])


cdef packed struct node_struct:
    Y_DTYPE_C value
    unsigned int count
    unsigned int feature_idx
    X_DTYPE_C threshold
    unsigned int left
    unsigned int right
    Y_DTYPE_C gain
    unsigned int depth
    unsigned char is_leaf
    X_BINNED_DTYPE_C bin_threshold


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
        out = np.empty(X.shape[0], dtype=Y_DTYPE)
        _predict_from_numeric_data(self.nodes, X, out)
        return out

    def predict_binned(self, X):
        """Predict raw values for binned data.

        Parameters
        ----------
        X : array-like, shape=(n_samples, n_features)
            The input samples.

        Returns
        -------
        y : array, shape (n_samples,)
            The raw predicted values.
        """
        out = np.empty(X.shape[0], dtype=Y_DTYPE)
        _predict_from_binned_data(self.nodes, X, out)
        return out

cdef inline Y_DTYPE_C _predict_one_from_numeric_data(
    node_struct [:] nodes,
    const X_DTYPE_C [:, :] numeric_data,
    const int row
    ) nogil:
    # Need to pass the whole array, else prange won't work. See issue Cython
    # #2798

    cdef:
        node_struct node = nodes[0]

    while True:
        if node.is_leaf:
            return node.value
        if numeric_data[row, node.feature_idx] <= node.threshold:
            node = nodes[node.left]
        else:
            node = nodes[node.right]


cdef void _predict_from_numeric_data(
    node_struct [:] nodes,
    const X_DTYPE_C [:, :] numeric_data,
    Y_DTYPE_C [:] out) nogil:

    cdef:
        int i

    for i in prange(numeric_data.shape[0], schedule='static'):
        out[i] = _predict_one_from_numeric_data(nodes, numeric_data, i)


cdef inline Y_DTYPE_C _predict_one_from_binned_data(
    node_struct [:] nodes,
    const X_BINNED_DTYPE_C [:, :] binned_data,
    const int row
    ) nogil:
    # Need to pass the whole array, else prange won't work. See issue Cython
    # #2798

    cdef:
        node_struct node = nodes[0]

    while True:
        if node.is_leaf:
            return node.value
        if binned_data[row, node.feature_idx] <= node.bin_threshold:
            node = nodes[node.left]
        else:
            node = nodes[node.right]


cdef void _predict_from_binned_data(
    node_struct [:] nodes,
    const X_BINNED_DTYPE_C [:, :] binned_data,
    Y_DTYPE_C [:] out) nogil:

    cdef:
        int i

    for i in prange(binned_data.shape[0], schedule='static'):
        out[i] = _predict_one_from_binned_data(nodes, binned_data, i)
