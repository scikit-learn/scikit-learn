"""
This module contains the TreePredictor class which is used for prediction.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from sklearn.ensemble._hist_gradient_boosting._predictor import (
    _compute_partial_dependence,
    _predict_from_binned_data,
    _predict_from_raw_data,
)
from sklearn.ensemble._hist_gradient_boosting.common import (
    PREDICTOR_RECORD_DTYPE,
    X_BITSET_LENGTH,
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
    """

    def __init__(self, nodes, binned_left_cat_bitsets, raw_left_cat_bitsets):
        self.nodes = nodes
        self.binned_left_cat_bitsets = binned_left_cat_bitsets
        self.raw_left_cat_bitsets = raw_left_cat_bitsets

    def get_n_leaf_nodes(self):
        """Return number of leaves."""
        return int(self.nodes["is_leaf"].sum())

    def get_max_depth(self):
        """Return maximum depth among all leaves."""
        return int(self.nodes["depth"].max())

    def predict(self, X, known_cat_bitsets, f_idx_map, n_threads):
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

        Returns
        -------
        y : ndarray, shape (n_samples,)
            The raw predicted values.
        """
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

    def predict_binned(self, X, missing_values_bin_idx, n_threads):
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

        Returns
        -------
        y : ndarray, shape (n_samples,)
            The raw predicted values.
        """
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
        out : ndarray, shape (n_samples)
            The value of the partial dependence function on each grid
            point.
        """
        _compute_partial_dependence(self.nodes, grid, target_features, out)

    def _check_feature_idx_bounds(self, n_features):
        """Validate that split nodes reference in-range feature columns.

        This complements the `left` / `right` / `bitset_idx` bounds validated
        at load time in `__setstate__`: `feature_idx` is dereferenced as a
        column index into `X` in the `nogil` traversal
        (`_predict_from_raw_data` / `_predict_from_binned_data`), which can
        cause a segfault if out of bounds. We use `_n_features` saved in `fit`
        time to validate `feature_idx` values, assuming it'd be the same during
        `predict`. Only split nodes dereference `feature_idx`, so only they are
        validated.
        """
        nodes = self.nodes
        is_split = ~nodes["is_leaf"].astype(bool)
        feature_idx = nodes["feature_idx"]
        if np.any(is_split & ((feature_idx < 0) | (feature_idx >= n_features))):
            raise ValueError(
                "predictor node array has out-of-bounds 'feature_idx' values: "
                f"expected each split node's to be in [0, {n_features}), got "
                f"{feature_idx}."
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

        # Bounds-check the node array so the nogil traversal in
        # `_predict_from_raw_data` / `_predict_from_binned_data` can't segfault
        # or loop forever on a malicious predictor. Leafness is the `is_leaf`
        # flag; only split nodes dereference `left` / `right` (node ids) and
        # `bitset_idx` (a row in the categorical bitset arrays). `feature_idx`
        # is checked separately in `_check_feature_idx_bounds`, which needs
        # `n_features`.
        nodes = self.nodes
        n_nodes = nodes.shape[0]
        node_index = np.arange(n_nodes)
        is_split = ~nodes["is_leaf"].astype(bool)

        # `left` / `right` must be node ids strictly greater than the own index:
        # this keeps them in `[0, n_nodes)` and rules out cycles.
        for field in ("left", "right"):
            child = nodes[field]
            if np.any(is_split & ((child <= node_index) | (child >= n_nodes))):
                raise ValueError(
                    f"predictor node array has out-of-bounds '{field}' values: "
                    f"expected each split node's to be a node id greater than "
                    f"the node's own index and less than {n_nodes}, got {child}."
                )

        # Together with the strict increase above, allowing a node to be the
        # child of at most one other node makes the array a tree rather than a
        # DAG. `_compute_partial_dependence` visits both children of a split and
        # sizes its stacks from the node count, so a node reachable through
        # several paths makes it do work exponential in the number of nodes.
        children = np.concatenate((nodes["left"][is_split], nodes["right"][is_split]))
        if children.size and np.bincount(children, minlength=n_nodes).max() > 1:
            raise ValueError(
                "predictor node array is inconsistent: expected each node to be "
                "the child of at most one node, got several nodes sharing a child."
            )

        # Categorical splits index the bitset arrays as
        # `bitsets[bitset_idx, binned_value // 32]`, where `binned_value` is a
        # uint8. Both the row (`bitset_idx`) and the implied column therefore
        # need to be in bounds, i.e. the arrays must have `X_BITSET_LENGTH`
        # columns and enough rows. Only touch those arrays when a categorical
        # split exists (a non-categorical predictor may carry empty ones).
        is_categorical_split = is_split & nodes["is_categorical"].astype(bool)
        if np.any(is_categorical_split):
            for name in ("binned_left_cat_bitsets", "raw_left_cat_bitsets"):
                bitsets = getattr(self, name, None)
                if (
                    bitsets is None
                    or bitsets.ndim != 2
                    or bitsets.shape[1] != X_BITSET_LENGTH
                ):
                    shape = None if bitsets is None else bitsets.shape
                    raise ValueError(
                        f"predictor '{name}' has an invalid shape: expected "
                        f"(n_categorical_splits, {X_BITSET_LENGTH}), got {shape}."
                    )

            bitset_idx = nodes["bitset_idx"]
            n_bitsets = min(
                self.binned_left_cat_bitsets.shape[0],
                self.raw_left_cat_bitsets.shape[0],
            )
            if np.any(is_categorical_split & (bitset_idx >= n_bitsets)):
                raise ValueError(
                    "predictor node array has out-of-bounds 'bitset_idx' values: "
                    f"expected each categorical split node's to be in "
                    f"[0, {n_bitsets}), got {bitset_idx}."
                )
