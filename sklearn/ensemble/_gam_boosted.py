"""Generalized Additive Models (GAMs) where bagging and gradient boosted trees
modeling the shape functions."""
from abc import ABC, abstractmethod
from itertools import product

import numpy as np

from ..preprocessing import LabelEncoder
from ..base import BaseEstimator, RegressorMixin, ClassifierMixin
from ._hist_gradient_boosting.common import X_DTYPE, Y_DTYPE
from ._hist_gradient_boosting.loss import _LOSSES
from ._hist_gradient_boosting.binning import _BinMapper
from ._hist_gradient_boosting.grower import TreeGrower
from ._hist_gradient_boosting._gradient_boosting import _update_raw_predictions

from ..utils import check_random_state
from ..utils.validation import check_is_fitted
from ..utils.multiclass import check_classification_targets


class BaseGAMBoosting(ABC, BaseEstimator):
    @abstractmethod
    def __init__(self, *, learning_rate=0.1,
                 max_iter=100, max_leaf_nodes=3, max_depth=2,
                 min_samples_leaf=2, random_state=None):
        self.max_iter = max_iter
        self.max_leaf_nodes = max_leaf_nodes
        self.max_depth = max_depth
        self.min_samples_leaf = min_samples_leaf
        self.learning_rate = learning_rate
        self.random_state = random_state

    @abstractmethod
    def _encode_y(self, y):
        """Encode target."""

    @abstractmethod
    def _get_loss(self, sample_weight):
        """Get loss used in boosting process."""

    def fit(self, X, y):

        # TODO: missing value support
        X, y = self._validate_data(X, y, dtype=X_DTYPE)
        y = self._encode_y(y)

        rng = check_random_state(self.random_state)
        self._random_seed = rng.randint(np.iinfo(np.uint32).max,
                                        dtype='u8')

        X_train, y_train = X, y

        n_bins = 256
        self._bin_mapper = _BinMapper(
            n_bins=n_bins,
            is_categorical=None,
            known_categories=None,
            random_state=self._random_seed)
        X_binned_train = self._bin_data(X_train, is_training_data=True)
        n_samples = X_binned_train.shape[0]

        # TODO: sample_weight support in fit. `sample_weight` is passed in
        # because the training loop will use sample weights for bagging.
        sample_weight = np.ones(n_samples)
        self._loss = self._get_loss(sample_weight=sample_weight)

        self._baseline_prediction = self._loss.get_baseline_prediction(
            y_train, None, self.n_trees_per_iteration_
        )
        raw_predictions = np.zeros(
            shape=(self.n_trees_per_iteration_, n_samples),
            dtype=Y_DTYPE,
        )
        raw_predictions += self._baseline_prediction

        gradients, hessians = self._loss.init_gradients_and_hessians(
            n_samples=n_samples,
            prediction_dim=self.n_trees_per_iteration_,
            sample_weight=sample_weight
        )

        # predictors are structured as nested list with the structure:
        # predictors[iteration][feature_idx][tree_idx]
        self._predictors = predictors = []

        if not hasattr(self, "_bagging_indices"):
            self._bagging_indices = check_random_state(self._random_seed)

        for _ in range(self.max_iter):
            # TODO: Out of bag early stopping?
            # sample weight for bagging
            sample_indices = self._bagging_indices.randint(0, high=n_samples,
                                                           size=n_samples)
            sample_weight_train = (np.bincount(sample_indices,
                                               minlength=n_samples)
                                   .astype(Y_DTYPE))

            # Update gradients and hessians, inplace
            self._loss.update_gradients_and_hessians(gradients, hessians,
                                                     y_train, raw_predictions,
                                                     sample_weight_train)
            predictors.append([[] for _ in range(self.n_features_in_)])

            tree_feature_iter = product(
                range(self.n_trees_per_iteration_), range(self.n_features_in_)
            )
            for tree_idx, feature_idx in tree_feature_iter:
                grower = TreeGrower(
                    X_binned_train, gradients[tree_idx, :],
                    hessians[tree_idx, :],
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
                _update_raw_predictions(raw_predictions[tree_idx, :], grower)

                # Update for the next feature
                self._loss.update_gradients_and_hessians(
                    gradients, hessians, y_train, raw_predictions,
                    sample_weight_train)

                predictor = grower.make_predictor(
                    binning_thresholds=self._bin_mapper.bin_thresholds_
                )
                predictors[-1][feature_idx].append(predictor)

        return self

    def _raw_predict(self, X, check_input=True, feature_idx=None):
        if check_input:
            X = self._validate_data(X, dtype=X_DTYPE, reset=False)

        n_samples = X.shape[0]
        raw_predictions = np.zeros(
            shape=(self.n_trees_per_iteration_, n_samples),
            dtype=Y_DTYPE,
        )

        # TODO: Does not actually change the shape, but this could be nicer?
        # Have the prediction contribute an equal amount to each feature
        if feature_idx is None:
            base_prediction = self._baseline_prediction
        else:
            base_prediction = self._baseline_prediction / self.n_features_in_

        raw_predictions += base_prediction
        self._predict_iterations(X, self._predictors, raw_predictions,
                                 feature_idx=feature_idx)
        return raw_predictions

    def _predict_iterations(self, X, predictors, raw_predictions,
                            feature_idx=None):
        known_cat_bitsets, f_idx_map = (
            self._bin_mapper.make_known_categories_bitsets())

        for predictors_for_iter in predictors:
            if feature_idx is None:
                predictor_iter = predictors_for_iter
            else:
                predictor_iter = [predictors_for_iter[feature_idx]]

            for predictors_for_feature in predictor_iter:
                for tree_idx, predictor in enumerate(predictors_for_feature):
                    raw_predictions[tree_idx, :] += predictor.predict(
                        X,
                        known_cat_bitsets=known_cat_bitsets,
                        f_idx_map=f_idx_map)

    def _bin_data(self, X, is_training_data):
        if is_training_data:
            X_binned = self._bin_mapper.fit_transform(X)  # F-aligned array
        else:
            X_binned = self._bin_mapper.transform(X)  # F-aligned array
            # We convert the array to C-contiguous since predicting is faster
            # with this layout (training is faster on F-arrays though)
            X_binned = np.ascontiguousarray(X_binned)
        return X_binned


class GAMBoostingRegressor(RegressorMixin, BaseGAMBoosting):
    """Generalized Additive Models (GAMs) with Bagged Histogram-based Gradient
    Boosting Regression Trees.

    Read more in the :ref:`User Guide <gam_boosted_trees>`.

    .. versionadded:: 1.0

    Parameters
    ----------
    learning_rate: float, default=0.1
        The learning rate, also known as *shrinkage*. This is used as a
        multiplicative factor for the leaves values. Use `1` for no
        shrinkage.
    max_iter : int, default=100
        The maximum number of iterations of the boosting process, i.e. the
        maximum number of trees.
    max_leaf_nodes : int, default=3
        The maximum number of leaves for each tree. Must be strictly greater
        than 1. If None, there is no maximum limit.
    max_depth : int, default=2
        The maximum depth of each tree. The depth of a tree is the number of
        edges to go from the root to the deepest leaf.
    min_samples_leaf : int, default=2
        The minimum number of samples per leaf. For small datasets with less
        than a few hundred samples, it is recommended to lower this value
        since only very shallow trees would be built.
    random_state : int, RandomState instance or None, default=None
        Pseudo-random number generator to control the subsampling in the
        binning process, and bagging for each iteration.
        Pass an int for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    Attributes
    ----------
    n_trees_per_iteration_ : int
        The number of tree that are built at each iteration. For regressors,
        this is always 1.
    n_features_in_ : int
        Number of features seen during :term:`fit`.

    See Also
    --------
    GAMBoostingClassifier
    """
    def __init__(self, *, learning_rate=0.1,
                 max_iter=100, max_leaf_nodes=3, max_depth=2,
                 min_samples_leaf=2, random_state=None):
        super().__init__(
            max_iter=max_iter,
            max_leaf_nodes=max_leaf_nodes,
            max_depth=max_depth,
            min_samples_leaf=min_samples_leaf,
            learning_rate=learning_rate,
            random_state=random_state,
        )

    def predict(self, X):
        """Predict values for X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The input samples.

        Returns
        -------
        y : ndarray, shape (n_samples,)
            The predicted values.
        """
        check_is_fitted(self)
        return self._loss.inverse_link_function(self._raw_predict(X).ravel())

    def apply(self, X, feature_idx):
        """Predict only for the shape function at `feature_idx`.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The input samples.
        feature_idx : int
            Index of shape function.

        Returns
        -------
        y : ndarray, shape (n_samples,)
            The predicted values for shape function.
        """
        check_is_fitted(self)
        return self._raw_predict(X, feature_idx=feature_idx).ravel()

    def _encode_y(self, y):
        self.n_trees_per_iteration_ = 1
        return y.astype(Y_DTYPE, copy=False)

    def _get_loss(self, sample_weight):
        return _LOSSES["squared_error"](sample_weight=sample_weight)


class GAMBoostingClassifier(ClassifierMixin, BaseGAMBoosting):
    """Generalized Additive Models (GAMs) with Bagged Histogram-based Gradient
    Boosting Regression Trees.

    Read more in the :ref:`User Guide <gam_boosted_trees>`.

    .. versionadded:: 1.0

    Parameters
    ----------
    learning_rate: float, default=0.1
        The learning rate, also known as *shrinkage*. This is used as a
        multiplicative factor for the leaves values. Use `1` for no
        shrinkage.
    max_iter : int, default=100
        The maximum number of iterations of the boosting process, i.e. the
        maximum number of trees.
    max_leaf_nodes : int, default=3
        The maximum number of leaves for each tree. Must be strictly greater
        than 1. If None, there is no maximum limit.
    max_depth : int, default=2
        The maximum depth of each tree. The depth of a tree is the number of
        edges to go from the root to the deepest leaf.
    min_samples_leaf : int, default=2
        The minimum number of samples per leaf. For small datasets with less
        than a few hundred samples, it is recommended to lower this value
        since only very shallow trees would be built.
    random_state : int, RandomState instance or None, default=None
        Pseudo-random number generator to control the subsampling in the
        binning process, and bagging for each iteration.
        Pass an int for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    Attributes
    ----------
    n_trees_per_iteration_ : int
        The number of tree that are built at each iteration. For regressors,
        this is always 1.
    n_features_in_ : int
        Number of features seen during :term:`fit`.

    See Also
    --------
    GAMBoostingRegressor
    """

    def __init__(self, *, learning_rate=0.1,
                 max_iter=100, max_leaf_nodes=3, max_depth=2,
                 min_samples_leaf=2, random_state=None):
        super().__init__(
            max_iter=max_iter,
            max_leaf_nodes=max_leaf_nodes,
            max_depth=max_depth,
            min_samples_leaf=min_samples_leaf,
            learning_rate=learning_rate,
            random_state=random_state,
        )

    def predict(self, X):
        """Predict classes for X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The input samples.

        Returns
        -------
        y : ndarray, shape (n_samples,)
            The predicted classes.
        """
        encoded_classes = np.argmax(self.predict_proba(X), axis=1)
        return self.classes_[encoded_classes]

    def predict_proba(self, X):
        """Predict class probabilities for X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The input samples.

        Returns
        -------
        p : ndarray, shape (n_samples, n_classes)
            The class probabilities of the input samples.
        """
        check_is_fitted(self)
        return self._loss.predict_proba(self._raw_predict(X))

    def _encode_y(self, y):
        # encode classes into 0 ... n_classes - 1 and sets attributes classes_
        # and n_trees_per_iteration_
        check_classification_targets(y)

        label_encoder = LabelEncoder()
        encoded_y = label_encoder.fit_transform(y)
        self.classes_ = label_encoder.classes_
        n_classes = self.classes_.shape[0]
        # only 1 tree for binary classification. For multiclass classification,
        # we build 1 tree per class.
        self.n_trees_per_iteration_ = 1 if n_classes <= 2 else n_classes
        encoded_y = encoded_y.astype(Y_DTYPE, copy=False)
        return encoded_y

    def _get_loss(self, sample_weight):
        if self.n_trees_per_iteration_ == 1:
            return _LOSSES['binary_crossentropy'](
                sample_weight=sample_weight)
        else:
            return _LOSSES['categorical_crossentropy'](
                sample_weight=sample_weight)
