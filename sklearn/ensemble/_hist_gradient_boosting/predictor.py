"""
This module contains the TreePredictor class which is used for prediction.
"""
# Author: Nicolas Hug

import numpy as np

from .common import Y_DTYPE
from .common import X_BITSET_INNER_DTYPE
from ._bitset import set_bitset_mv
from ._predictor import _predict_from_data
from ._predictor import _predict_from_binned_data
from ._predictor import _compute_partial_dependence


class TreePredictor:
    """Tree class used for predictions.

    Parameters
    ----------
    nodes : ndarray of PREDICTOR_RECORD_DTYPE
        The nodes of the tree.
    binned_categorical_bitsets : ndarray of shape (n_categorical, 8), \
            dtype=uint32
        Bitset for binned categorical used in predict_binned.
    raw_categorical_bitsets : ndarray of shape (n_categorical, 8), \
            dtype=uint32
        Bitset for raw categorical used in predict.

    """
    def __init__(self, nodes, binned_categorical_bitsets,
                 raw_categorical_bitsets):
        self.nodes = nodes
        self.binned_categorical_bitsets = binned_categorical_bitsets
        self.raw_categorical_bitsets = raw_categorical_bitsets

    def get_n_leaf_nodes(self):
        """Return number of leaves."""
        return int(self.nodes['is_leaf'].sum())

    def get_max_depth(self):
        """Return maximum depth among all leaves."""
        return int(self.nodes['depth'].max())

    def predict(self, X, known_cat_bitset, orig_feat_to_known_cats_idx):
        """Predict raw values for non-binned numerical data and binned
        categorical data.

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_features)
            The input samples.

        known_cat_bitset : ndarray of shape (n_categorical, 8)
            Bitset for known categories.

        orig_idx_to_cat_idx : ndarray of shape (n_features,)
            Maps from original feature idx to the categorical idx in
            known_cat_bitset.

        Returns
        -------
        y : ndarray, shape (n_samples,)
            The raw predicted values.
        """
        out = np.empty(X.shape[0], dtype=Y_DTYPE)
        _predict_from_data(self.nodes, X, self.raw_categorical_bitsets,
                           known_cat_bitset, orig_feat_to_known_cats_idx, out)
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
        _predict_from_binned_data(self.nodes, X,
                                  self.binned_categorical_bitsets,
                                  missing_values_bin_idx, out)
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


def _make_known_categories(bin_thresholds, is_categorical):
    n_features = len(bin_thresholds)
    categorical_indices = np.flatnonzero(is_categorical)
    orig_feat_to_known_cats_idx = np.zeros(n_features, dtype=np.uint8)

    orig_feat_to_known_cats_idx[categorical_indices] = \
        np.arange(categorical_indices.size, dtype=np.uint8)

    known_cat_bitset = np.zeros((categorical_indices.size, 8),
                                dtype=X_BITSET_INNER_DTYPE)
    for idx, cat_idx in enumerate(categorical_indices):
        for num in bin_thresholds[cat_idx]:
            set_bitset_mv(known_cat_bitset[idx], num)

    return known_cat_bitset, orig_feat_to_known_cats_idx
