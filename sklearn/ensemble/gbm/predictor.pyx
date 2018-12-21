# cython: profile=True
"""
This module contains the TreePredictor class which is used for prediction.
"""
import numpy as np


PREDICTOR_RECORD_DTYPE = np.dtype([
    ('value', np.float32),
    ('count', np.uint32),
    ('feature_idx', np.uint32),
    ('threshold', np.float32),
    ('left', np.uint32),
    ('right', np.uint32),
    ('gain', np.float32),
    ('depth', np.uint32),
    ('is_leaf', np.uint8),
    ('bin_threshold', np.uint8),
    # TODO: shrinkage in leaf for feature importance error bar?
])

ctypedef fused float_or_double:
    float
    double

cdef packed struct node_struct:
    float value
    unsigned int count
    unsigned int feature_idx
    float threshold
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

    def predict_binned(self, binned_data, out=None):
        """Predict raw values for binned data.

        Parameters
        ----------
        binned_data : array-like of np.uint8, shape=(n_samples, n_features)
            The binned input samples.
        out : array-like, shape=(n_samples,), optional (default=None)
            If not None, predictions will be written inplace in ``out``.

        Returns
        -------
        y : array, shape (n_samples,)
            The raw predicted values.
        """
        if out is None:
            out = np.empty(binned_data.shape[0], dtype=np.float32)
        _predict_binned(self.nodes, binned_data, out)
        return out

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
        # TODO: introspect X to dispatch to numerical or categorical data
        # (dense or sparse) on a feature by feature basis.
        out = np.empty(X.shape[0], dtype=np.float32)
        _predict_from_numeric_data(self.nodes, X, out)
        return out


def _predict_one_binned(nodes, binned_data):
    node = nodes[0]
    while True:
        if node['is_leaf']:
            return node['value']
        if binned_data[node['feature_idx']] <= node['bin_threshold']:
            node = nodes[node['left']]
        else:
            node = nodes[node['right']]


def _predict_binned(nodes, binned_data, out):
    for i in range(binned_data.shape[0]):
        out[i] = _predict_one_binned(nodes, binned_data[i])


cdef float _predict_one_from_numeric_data(node_struct [:] nodes, float_or_double [:] numeric_data) nogil:
    cdef node_struct node = nodes[0]
    while True:
        if node.is_leaf:
            return node.value
        if numeric_data[node.feature_idx] <= node.threshold:
            node = nodes[node.left]
        else:
            node = nodes[node.right]


# TODO: having a view on numeric_data (passed by user) may not be supported,
# see sklearn issue 10624
def _predict_from_numeric_data(node_struct [:] nodes, float_or_double [:, :] numeric_data, float [:] out):

    cdef int i

    for i in range(numeric_data.shape[0]):
        out[i] = _predict_one_from_numeric_data(nodes, numeric_data[i])
