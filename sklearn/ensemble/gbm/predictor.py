"""
This module contains the TreePredictor class which is used for prediction.
"""
import numpy as np


PREDICTOR_RECORD_DTYPE = np.dtype([
    ('is_leaf', np.uint8),
    ('value', np.float32),
    ('count', np.uint32),
    ('feature_idx', np.uint32),
    ('bin_threshold', np.uint8),
    ('threshold', np.float32),
    ('left', np.uint32),
    ('right', np.uint32),
    ('gain', np.float32),
    ('depth', np.uint32),
    # TODO: shrinkage in leaf for feature importance error bar?
])


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


def _predict_one_from_numeric_data(nodes, numeric_data):
    node = nodes[0]
    while True:
        if node['is_leaf']:
            return node['value']
        if numeric_data[node['feature_idx']] <= node['threshold']:
            node = nodes[node['left']]
        else:
            node = nodes[node['right']]


def _predict_from_numeric_data(nodes, numeric_data, out):
    for i in range(numeric_data.shape[0]):
        out[i] = _predict_one_from_numeric_data(nodes, numeric_data[i])
