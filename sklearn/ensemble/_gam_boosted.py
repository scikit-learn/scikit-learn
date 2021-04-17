"""Generalized Additive Models (GAMs) where bagging and gradient boosted trees
modeling the shape functions."""
import numpy as np

from ..base import BaseEstimator
from ._hist_gradient_boosting.common import X_DTYPE, Y_DTYPE
from ._hist_gradient_boosting.loss import _LOSSES
from ._hist_gradient_boosting.binning import _BinMapper
from ._hist_gradient_boosting.grower import TreeGrower
from ._hist_gradient_boosting._gradient_boosting import _update_raw_predictions

from ..utils import check_random_state
from ..utils.validation import check_is_fitted


class GAMBoostingRegressor(BaseEstimator):
    """Generalized Additive Models (GAMs) with Bagged Gradient boosted trees.

    Parameters
    ----------
    max_iter : int, default=100

    max_leaf_nodes : int, default=3

    max_depth : int, default=None

    min_samples_leaf : int, default=20

    learning_rate: float, default=0.1

    random_state : int, RandomState instance or None, default=None

    """
    def __init__(self, max_iter=100, max_leaf_nodes=3, max_depth=None,
                 min_samples_leaf=20, learning_rate=0.1,
                 random_state=None):
        self.max_iter = max_iter
        self.max_leaf_nodes = max_leaf_nodes
        self.max_depth = max_depth
        self.min_samples_leaf = min_samples_leaf
        self.learning_rate = learning_rate
        self.random_state = random_state

    def _bin_data(self, X, is_training_data):
        if is_training_data:
            X_binned = self._bin_mapper.fit_transform(X)  # F-aligned array
        else:
            X_binned = self._bin_mapper.transform(X)  # F-aligned array
            # We convert the array to C-contiguous since predicting is faster
            # with this layout (training is faster on F-arrays though)
            X_binned = np.ascontiguousarray(X_binned)
        return X_binned

    def fit(self, X, y):

        # TODO: do not support missing values for now
        X, y = self._validate_data(X, y, dtype=X_DTYPE)
        # bin data

        rng = check_random_state(self.random_state)
        self._random_seed = rng.randint(np.iinfo(np.uint32).max,
                                        dtype='u8')

        # TODO: early stopping
        X_train, y_train = X, y

        # TODO: n_bins
        # TODO: Categorical support
        n_bins = 256
        self._bin_mapper = _BinMapper(
            n_bins=n_bins,
            is_categorical=None,
            known_categories=None,
            random_state=self._random_seed)
        X_binned_train = self._bin_data(X_train, is_training_data=True)

        n_samples, self.n_features_in_ = X_binned_train.shape

        sample_weight = np.ones(n_samples)
        self._loss = _LOSSES["squared_error"](sample_weight=sample_weight)

        self.n_trees_per_iteration_ = 1

        raw_predictions = np.zeros(
            shape=(self.n_trees_per_iteration_, n_samples),
            dtype=Y_DTYPE,
        )
        self._predictors = predictors = []

        gradients, hessians = self._loss.init_gradients_and_hessians(
            n_samples=n_samples,
            prediction_dim=self.n_trees_per_iteration_,
            sample_weight=None
        )

        indices_rng = check_random_state(self._random_seed)
        for iteration in range(self.max_iter):
            # sample weight for bagging
            sample_indices = indices_rng.randint(0, high=n_samples,
                                                 size=n_samples)
            sample_weight_train = (np.bincount(sample_indices,
                                               minlength=n_samples)
                                   .astype(Y_DTYPE))

            # Update gradients and hessians, inplace
            self._loss.update_gradients_and_hessians(gradients, hessians,
                                                     y_train, raw_predictions,
                                                     sample_weight_train)
            predictors.append([[] for _ in range(self.n_features_in_)])

            for feature_idx in range(self.n_features_in_):
                grower = TreeGrower(
                    X_binned_train, gradients[0, :], hessians[0, :],
                    n_bins=n_bins,
                    n_bins_non_missing=self._bin_mapper.n_bins_non_missing_,
                    has_missing_values=False,
                    is_categorical=None,
                    monotonic_cst=None,
                    max_leaf_nodes=self.max_leaf_nodes,
                    max_depth=self.max_depth,
                    min_samples_leaf=self.min_samples_leaf,
                    l2_regularization=0,
                    shrinkage=self.learning_rate,
                    feature_idx=feature_idx,
                )
                grower.grow()
                _update_raw_predictions(raw_predictions[0, :], grower)

                # Update for the next feature
                self._loss.update_gradients_and_hessians(
                    gradients, hessians, y_train, raw_predictions,
                    sample_weight_train)

                predictor = grower.make_predictor(
                    binning_thresholds=self._bin_mapper.bin_thresholds_
                )
                predictors[-1][feature_idx].append(predictor)

        return self

    def predict(self, X):
        check_is_fitted(self)
        X = self._validate_data(X, dtype=X_DTYPE, reset=False)
        n_samples = X.shape[0]
        raw_predictions = np.zeros(
            shape=(self.n_trees_per_iteration_, n_samples),
            dtype=Y_DTYPE,
        )
        self._predict_iterations(X, self._predictors, raw_predictions)
        return raw_predictions.ravel()

    def apply(self, X, feature_idx):
        check_is_fitted(self)

        X = self._validate_data(X, dtype=X_DTYPE, reset=False)
        n_samples = X.shape[0]
        raw_predictions = np.zeros(
            shape=(self.n_trees_per_iteration_, n_samples),
            dtype=Y_DTYPE
        )
        self._predict_iterations(X, self._predictors, raw_predictions,
                                 feature_idx)
        return raw_predictions.ravel()

    def _predict_iterations(self, X, predictors, raw_predictions,
                            feature_idx=None):
        known_cat_bitsets, f_idx_map = (
            self._bin_mapper.make_known_categories_bitsets())

        for predictors_of_iter in predictors:
            if feature_idx is None:
                predictor_iter = predictors_of_iter
            else:
                predictor_iter = [predictors_of_iter[feature_idx]]

            for predictors_per_feature in predictor_iter:
                for predictor_for_tree in predictors_per_feature:
                    raw_predictions[0, :] += predictor_for_tree.predict(
                        X,
                        known_cat_bitsets=known_cat_bitsets,
                        f_idx_map=f_idx_map)
