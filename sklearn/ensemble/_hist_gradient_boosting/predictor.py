"""
This module contains the TreePredictor class which is used for prediction.
"""
# Author: Nicolas Hug

import numpy as np

from .types import X_DTYPE
from .types import Y_DTYPE
from .types import X_BINNED_DTYPE
from ._predictor import _predict_from_numeric_data
from ._predictor import _predict_from_binned_data


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


class TreePredictor:
    """Tree class used for predictions.

    Parameters
    ----------
    nodes : list of PREDICTOR_RECORD_DTYPE
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
        X : ndarray, shape (n_samples, n_features)
            The input samples.

        Returns
        -------
        y : ndarray, shape (n_samples,)
            The raw predicted values.
        """
        out = np.empty(X.shape[0], dtype=Y_DTYPE)
        _predict_from_numeric_data(self.nodes, X, out)
        return out

    def predict_binned(self, X):
        """Predict raw values for binned data.

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_features)
            The input samples.

        Returns
        -------
        y : ndarray, shape (n_samples,)
            The raw predicted values.
        """
        out = np.empty(X.shape[0], dtype=Y_DTYPE)
        _predict_from_binned_data(self.nodes, X, out)
        return out
