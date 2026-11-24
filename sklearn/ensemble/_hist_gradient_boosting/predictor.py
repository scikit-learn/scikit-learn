"""
This module contains the TreePredictor class which is used for prediction.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from sklearn.ensemble._hist_gradient_boosting._predictor import (
    _compute_partial_dependence,
    _predict_from_binned_data,
    _predict_from_binned_data_vector,
    _predict_from_raw_data,
    _predict_from_raw_data_vector,
)
from sklearn.ensemble._hist_gradient_boosting.common import (
    PREDICTOR_RECORD_DTYPE,
    Y_DTYPE,
)


class TreePredictor:
    """Tree class used for predictions.

    Parameters
    ----------
    nodes : ndarray of PREDICTOR_RECORD_DTYPE
        The nodes of the tree.
    binned_left_cat_bitsets : ndarray of shape (n_categorical_splits, 8), dtype=uint32
        Array of bitsets for binned categories used in predict_binned when a
        split is categorical.
    raw_left_cat_bitsets : ndarray of shape (n_categorical_splits, 8), dtype=uint32
        Array of bitsets for raw categories used in predict when a split is
        categorical.
    node_values : ndarray of shape (n_nodes, n_outputs), dtype=float64, \
            default=None
        The values of the nodes of a tree with vector leaves. If None, the
        tree has scalar leaves whose values are stored in ``nodes["value"]``.
    """

    def __init__(
        self, nodes, binned_left_cat_bitsets, raw_left_cat_bitsets, node_values=None
    ):
        self.nodes = nodes
        self.binned_left_cat_bitsets = binned_left_cat_bitsets
        self.raw_left_cat_bitsets = raw_left_cat_bitsets
        self.node_values = node_values

    @property
    def n_outputs(self):
        """Number of outputs of a tree with vector leaves, None otherwise."""
        if self.node_values is None:
            return None
        return self.node_values.shape[1]

    def get_n_leaf_nodes(self):
        """Return number of leaves."""
        return int(self.nodes["is_leaf"].sum())

    def get_max_depth(self):
        """Return maximum depth among all leaves."""
        return int(self.nodes["depth"].max())

    def predict(self, X, known_cat_bitsets, f_idx_map, n_threads, out=None):
        """Predict raw values for non-binned data.

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_features)
            The input samples.

        known_cat_bitsets : ndarray of shape (n_categorical_features, 8)
            Array of bitsets of known categories, for each categorical feature.

        f_idx_map : ndarray of shape (n_features,)
            Map from original feature index to the corresponding index in the
            known_cat_bitsets array.

        n_threads : int
            Number of OpenMP threads to use.

        out : ndarray of shape (n_samples, n_outputs), default=None
            Only for trees with vector leaves: if not None, the predicted
            values are added to ``out`` in place.

        Returns
        -------
        y : ndarray, shape (n_samples,) or (n_samples, n_outputs)
            The raw predicted values. For trees with vector leaves, this is
            ``out`` if it is not None.
        """
        if self.node_values is not None:
            if out is None:
                out = np.zeros((X.shape[0], self.n_outputs), dtype=Y_DTYPE)
            _predict_from_raw_data_vector(
                self.nodes,
                self.node_values,
                X,
                self.raw_left_cat_bitsets,
                known_cat_bitsets,
                f_idx_map,
                n_threads,
                out,
            )
            return out

        out = np.empty(X.shape[0], dtype=Y_DTYPE)

        _predict_from_raw_data(
            self.nodes,
            X,
            self.raw_left_cat_bitsets,
            known_cat_bitsets,
            f_idx_map,
            n_threads,
            out,
        )
        return out

    def predict_binned(self, X, missing_values_bin_idx, n_threads, out=None):
        """Predict raw values for binned data.

        Parameters
        ----------
        X : ndarray, shape (n_samples, n_features)
            The input samples.
        missing_values_bin_idx : uint8
            Index of the bin that is used for missing values. This is the
            index of the last bin and is always equal to max_bins (as passed
            to the GBDT classes), or equivalently to n_bins - 1.
        n_threads : int
            Number of OpenMP threads to use.
        out : ndarray of shape (n_samples, n_outputs), default=None
            Only for trees with vector leaves: if not None, the predicted
            values are added to ``out`` in place.

        Returns
        -------
        y : ndarray, shape (n_samples,) or (n_samples, n_outputs)
            The raw predicted values. For trees with vector leaves, this is
            ``out`` if it is not None.
        """
        if self.node_values is not None:
            if out is None:
                out = np.zeros((X.shape[0], self.n_outputs), dtype=Y_DTYPE)
            _predict_from_binned_data_vector(
                self.nodes,
                self.node_values,
                X,
                self.binned_left_cat_bitsets,
                missing_values_bin_idx,
                n_threads,
                out,
            )
            return out

        out = np.empty(X.shape[0], dtype=Y_DTYPE)
        _predict_from_binned_data(
            self.nodes,
            X,
            self.binned_left_cat_bitsets,
            missing_values_bin_idx,
            n_threads,
            out,
        )
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
        out : ndarray, shape (n_samples) or (n_outputs, n_samples)
            The value of the partial dependence function on each grid
            point. 2-dimensional for trees with vector leaves. Values are
            added to ``out`` in place.
        """
        if self.node_values is None:
            node_values = np.ascontiguousarray(self.nodes["value"]).reshape(-1, 1)
            out_2d = out.reshape(1, -1)
        else:
            node_values = self.node_values
            out_2d = out
        _compute_partial_dependence(
            self.nodes, node_values, grid, target_features, out_2d
        )

    def __setstate__(self, state):
        try:
            super().__setstate__(state)
        except AttributeError:
            self.__dict__.update(state)

        # The dtype of feature_idx is np.intp which is platform dependent. Here, we
        # make sure that saving and loading on different bitness systems works without
        # errors. For instance, on a 64 bit Python runtime, np.intp = np.int64,
        # while on 32 bit np.intp = np.int32.
        #
        # TODO: consider always using platform agnostic dtypes for fitted
        # estimator attributes. For this particular estimator, this would
        # mean replacing the intp field of PREDICTOR_RECORD_DTYPE by an int32
        # field. Ideally this should be done consistently throughout
        # scikit-learn along with a common test.
        if self.nodes.dtype != PREDICTOR_RECORD_DTYPE:
            self.nodes = self.nodes.astype(PREDICTOR_RECORD_DTYPE, casting="same_kind")

        # Predictors pickled before vector leaves were introduced.
        if not hasattr(self, "node_values"):
            self.node_values = None
