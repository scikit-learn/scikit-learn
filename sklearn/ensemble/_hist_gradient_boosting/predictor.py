"""
This module contains the TreePredictor class which is used for prediction.
"""
# Author: Nicolas Hug

import numpy as np

from .common import Y_DTYPE
from ._predictor import _predict_from_numeric_data
from ._predictor import _predict_from_binned_data
from ._predictor import _compute_partial_dependence


class TreePredictor:
    """Tree class used for predictions.

    Parameters
    ----------
    nodes : ndarray of PREDICTOR_RECORD_DTYPE
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

    def predict_binned(self, X, missing_values_bin_idx):
        """Predict raw values for binned data.

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_features)
            The input samples.
        missing_values_bin_idx : uint8
            Index of the bin that is used for missing values. This is the
            index of the last bin and is always equal to max_bins (as passed
            to the GBDT classes), or equivalently to n_bins - 1.

        Returns
        -------
        y : ndarray, shape (n_samples,)
            The raw predicted values.
        """
        out = np.empty(X.shape[0], dtype=Y_DTYPE)
        _predict_from_binned_data(self.nodes, X, missing_values_bin_idx, out)
        return out

    def compute_partial_dependence(self, grid, target_features, out):
        """Fast partial dependence computation.

        Parameters
        ----------
        grid : ndarray, shape (n_samples, n_target_features)
            The grid points on which the partial dependence should be
            evaluated.
        target_features : ndarray, shape (n_target_features)
            The set of target features for which the partial dependence
            should be evaluated.
        out : ndarray, shape (n_samples)
            The value of the partial dependence function on each grid
            point.
        """
        _compute_partial_dependence(self.nodes, grid, target_features, out)
