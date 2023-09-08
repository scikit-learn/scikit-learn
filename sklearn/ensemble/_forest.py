"""
Forest of trees-based ensemble methods.

Those methods include random forests and extremely randomized trees.

The module structure is the following:

- The ``BaseForest`` base class implements a common ``fit`` method for all
  the estimators in the module. The ``fit`` method of the base ``Forest``
  class calls the ``fit`` method of each sub-estimator on random samples
  (with replacement, a.k.a. bootstrap) of the training set.

  The init of the sub-estimator is further delegated to the
  ``BaseEnsemble`` constructor.

- The ``ForestClassifier`` and ``ForestRegressor`` base classes further
  implement the prediction logic by computing an average of the predicted
  outcomes of the sub-estimators.

- The ``RandomForestClassifier`` and ``RandomForestRegressor`` derived
  classes provide the user with concrete implementations of
  the forest ensemble method using classical, deterministic
  ``DecisionTreeClassifier`` and ``DecisionTreeRegressor`` as
  sub-estimator implementations.

- The ``ExtraTreesClassifier`` and ``ExtraTreesRegressor`` derived
  classes provide the user with concrete implementations of the
  forest ensemble method using the extremely randomized trees
  ``ExtraTreeClassifier`` and ``ExtraTreeRegressor`` as
  sub-estimator implementations.

Single and multi-output problems are both handled.
"""

# Authors: Gilles Louppe <g.louppe@gmail.com>
#          Brian Holt <bdholt1@gmail.com>
#          Joly Arnaud <arnaud.v.joly@gmail.com>
#          Fares Hedayati <fares.hedayati@gmail.com>
#
# License: BSD 3 clause


import threading
from abc import ABCMeta, abstractmethod
from numbers import Integral, Real
from time import time
from warnings import catch_warnings, simplefilter, warn

import numpy as np
from scipy.sparse import hstack as sparse_hstack
from scipy.sparse import issparse

from sklearn.base import (
    ClassifierMixin,
    MultiOutputMixin,
    RegressorMixin,
    TransformerMixin,
    _fit_context,
    is_classifier,
)
from sklearn.ensemble._base import BaseEnsemble, _partition_estimators
from sklearn.ensemble._hist_gradient_boosting.binning import _BinMapper
from sklearn.exceptions import DataConversionWarning
from sklearn.metrics import accuracy_score, r2_score
from sklearn.preprocessing import OneHotEncoder
from sklearn.utils import check_random_state, compute_sample_weight
from sklearn.utils._openmp_helpers import _openmp_effective_n_threads
from sklearn.utils._param_validation import Interval, RealNotInt, StrOptions
from sklearn.utils._tags import _safe_tags
from sklearn.utils.multiclass import (
    _check_partial_fit_first_call,
    check_classification_targets,
    type_of_target,
)
from sklearn.utils.parallel import Parallel, delayed
from sklearn.utils.validation import (
    _check_feature_names_in,
    _check_sample_weight,
    _num_samples,
    check_is_fitted,
)

from ..tree import (
    BaseDecisionTree,
    DecisionTreeClassifier,
    DecisionTreeRegressor,
    ExtraTreeClassifier,
    ExtraTreeRegressor,
)
from ..tree._tree import DOUBLE, DTYPE

__all__ = [
    "RandomForestClassifier",
    "RandomForestRegressor",
    "ExtraTreesClassifier",
    "ExtraTreesRegressor",
    "RandomTreesEmbedding",
]

MAX_INT = np.iinfo(np.int32).max


def _get_n_samples_bootstrap(n_samples, max_samples):
    """
    Get the number of samples in a bootstrap sample.

    Parameters
    ----------
    n_samples : int
        Number of samples in the dataset.
    max_samples : int or float
        The maximum number of samples to draw from the total available:
            - if float, this indicates a fraction of the total and should be
              the interval `(0.0, 1.0]`;
            - if int, this indicates the exact number of samples;
            - if None, this indicates the total number of samples.

    Returns
    -------
    n_samples_bootstrap : int
        The total number of samples to draw for the bootstrap sample.
    """
    if max_samples is None:
        return n_samples

    if isinstance(max_samples, Integral):
        if max_samples > n_samples:
            msg = "`max_samples` must be <= n_samples={} but got value {}"
            raise ValueError(msg.format(n_samples, max_samples))
        return max_samples

    if isinstance(max_samples, Real):
        return max(round(n_samples * max_samples), 1)


def _generate_sample_indices(random_state, n_samples, n_samples_bootstrap):
    """
    Private function used to _parallel_build_trees function."""

    random_instance = check_random_state(random_state)
    sample_indices = random_instance.randint(0, n_samples, n_samples_bootstrap)

    return sample_indices


def _generate_unsampled_indices(random_state, n_samples, n_samples_bootstrap):
    """
    Private function used to forest._set_oob_score function."""
    sample_indices = _generate_sample_indices(
        random_state, n_samples, n_samples_bootstrap
    )
    sample_counts = np.bincount(sample_indices, minlength=n_samples)
    unsampled_mask = sample_counts == 0
    indices_range = np.arange(n_samples)
    unsampled_indices = indices_range[unsampled_mask]

    return unsampled_indices


def _parallel_build_trees(
    tree,
    bootstrap,
    X,
    y,
    sample_weight,
    tree_idx,
    n_trees,
    verbose=0,
    class_weight=None,
    n_samples_bootstrap=None,
    missing_values_in_feature_mask=None,
    classes=None,
):
    """
    Private function used to fit a single tree in parallel."""
    if verbose > 1:
        print("building tree %d of %d" % (tree_idx + 1, n_trees))

    if bootstrap:
        n_samples = X.shape[0]
        if sample_weight is None:
            curr_sample_weight = np.ones((n_samples,), dtype=np.float64)
        else:
            curr_sample_weight = sample_weight.copy()

        indices = _generate_sample_indices(
            tree.random_state, n_samples, n_samples_bootstrap
        )
        sample_counts = np.bincount(indices, minlength=n_samples)
        curr_sample_weight *= sample_counts

        if class_weight == "subsample":
            with catch_warnings():
                simplefilter("ignore", DeprecationWarning)
                curr_sample_weight *= compute_sample_weight("auto", y, indices=indices)
        elif class_weight == "balanced_subsample":
            curr_sample_weight *= compute_sample_weight("balanced", y, indices=indices)

        tree._fit(
            X,
            y,
            sample_weight=curr_sample_weight,
            check_input=False,
            missing_values_in_feature_mask=missing_values_in_feature_mask,
            classes=classes,
        )
    else:
        tree._fit(
            X,
            y,
            sample_weight=sample_weight,
            check_input=False,
            missing_values_in_feature_mask=missing_values_in_feature_mask,
            classes=classes,
        )

    return tree


def _parallel_update_trees(
    tree,
    bootstrap,
    X,
    y,
    sample_weight,
    tree_idx,
    n_trees,
    verbose=0,
    class_weight=None,
    n_samples_bootstrap=None,
    classes=None,
):
    """
    Private function used to fit a single tree in parallel."""
    if verbose > 1:
        print("Updating tree %d of %d" % (tree_idx + 1, n_trees))

    if bootstrap:
        n_samples = X.shape[0]
        indices = _generate_sample_indices(
            tree.random_state, n_samples, n_samples_bootstrap
        )

        tree.partial_fit(
            X[indices, :],
            y[indices],
            sample_weight=sample_weight,
            check_input=False,
            classes=classes,
        )
    else:
        tree.partial_fit(
            X,
            y,
            sample_weight=sample_weight,
            check_input=False,
            classes=classes,
        )

    return tree


class BaseForest(MultiOutputMixin, BaseEnsemble, metaclass=ABCMeta):
    """
    Base class for forests of trees.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """

    _parameter_constraints: dict = {
        "n_estimators": [Interval(Integral, 1, None, closed="left")],
        "bootstrap": ["boolean"],
        "oob_score": ["boolean", callable],
        "n_jobs": [Integral, None],
        "random_state": ["random_state"],
        "verbose": ["verbose"],
        "warm_start": ["boolean"],
        "max_samples": [
            None,
            Interval(RealNotInt, 0.0, 1.0, closed="right"),
            Interval(Integral, 1, None, closed="left"),
        ],
        "max_bins": [
            None,
            Interval(Integral, 1, None, closed="left"),
        ],
        "store_leaf_values": ["boolean"],
    }

    @abstractmethod
    def __init__(
        self,
        estimator,
        n_estimators=100,
        *,
        estimator_params=tuple(),
        bootstrap=False,
        oob_score=False,
        n_jobs=None,
        random_state=None,
        verbose=0,
        warm_start=False,
        class_weight=None,
        max_samples=None,
        base_estimator="deprecated",
        max_bins=None,
        store_leaf_values=False,
    ):
        super().__init__(
            estimator=estimator,
            n_estimators=n_estimators,
            estimator_params=estimator_params,
            base_estimator=base_estimator,
        )

        self.bootstrap = bootstrap
        self.oob_score = oob_score
        self.n_jobs = n_jobs
        self.random_state = random_state
        self.verbose = verbose
        self.warm_start = warm_start
        self.class_weight = class_weight
        self.max_samples = max_samples
        self.max_bins = max_bins
        self.store_leaf_values = store_leaf_values

    def apply(self, X):
        """
        Apply trees in the forest to X, return leaf indices.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The input samples. Internally, its dtype will be converted to
            ``dtype=np.float32``. If a sparse matrix is provided, it will be
            converted into a sparse ``csr_matrix``.

        Returns
        -------
        X_leaves : ndarray of shape (n_samples, n_estimators)
            For each datapoint x in X and for each tree in the forest,
            return the index of the leaf x ends up in.
        """
        X = self._validate_X_predict(X)

        # if we trained a binning tree, then we should re-bin the data
        # XXX: this is inefficient and should be improved to be in line with what
        # the Histogram Gradient Boosting Tree does, where the binning thresholds
        # are passed into the tree itself, thus allowing us to set the node feature
        # value thresholds within the tree itself.
        if self.max_bins is not None:
            X = self._bin_data(X, is_training_data=False).astype(DTYPE)

        results = Parallel(
            n_jobs=self.n_jobs,
            verbose=self.verbose,
            prefer="threads",
        )(delayed(tree.apply)(X, check_input=False) for tree in self.estimators_)

        return np.array(results).T

    def decision_path(self, X):
        """
        Return the decision path in the forest.

        .. versionadded:: 0.18

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The input samples. Internally, its dtype will be converted to
            ``dtype=np.float32``. If a sparse matrix is provided, it will be
            converted into a sparse ``csr_matrix``.

        Returns
        -------
        indicator : sparse matrix of shape (n_samples, n_nodes)
            Return a node indicator matrix where non zero elements indicates
            that the samples goes through the nodes. The matrix is of CSR
            format.

        n_nodes_ptr : ndarray of shape (n_estimators + 1,)
            The columns from indicator[n_nodes_ptr[i]:n_nodes_ptr[i+1]]
            gives the indicator value for the i-th estimator.
        """
        X = self._validate_X_predict(X)
        indicators = Parallel(
            n_jobs=self.n_jobs,
            verbose=self.verbose,
            prefer="threads",
        )(
            delayed(tree.decision_path)(X, check_input=False)
            for tree in self.estimators_
        )

        n_nodes = [0]
        n_nodes.extend([i.shape[1] for i in indicators])
        n_nodes_ptr = np.array(n_nodes).cumsum()

        return sparse_hstack(indicators).tocsr(), n_nodes_ptr

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y, sample_weight=None, classes=None):
        """
        Build a forest of trees from the training set (X, y).

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The training input samples. Internally, its dtype will be converted
            to ``dtype=np.float32``. If a sparse matrix is provided, it will be
            converted into a sparse ``csc_matrix``.

        y : array-like of shape (n_samples,) or (n_samples, n_outputs)
            The target values (class labels in classification, real numbers in
            regression).

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights. If None, then samples are equally weighted. Splits
            that would create child nodes with net zero or negative weight are
            ignored while searching for a split in each node. In the case of
            classification, splits are also ignored if they would result in any
            single class carrying a negative weight in either child node.

        classes : array-like of shape (n_classes,), default=None
            List of all the classes that can possibly appear in the y vector.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        # Validate or convert input data
        if issparse(y):
            raise ValueError("sparse multilabel-indicator for y is not supported.")

        X, y = self._validate_data(
            X,
            y,
            multi_output=True,
            accept_sparse="csc",
            dtype=DTYPE,
            force_all_finite=False,
        )
        # _compute_missing_values_in_feature_mask checks if X has missing values and
        # will raise an error if the underlying tree base estimator can't handle missing
        # values. Only the criterion is required to determine if the tree supports
        # missing values.
        estimator = type(self.estimator)(criterion=self.criterion)
        missing_values_in_feature_mask = (
            estimator._compute_missing_values_in_feature_mask(
                X, estimator_name=self.__class__.__name__
            )
        )

        if sample_weight is not None:
            sample_weight = _check_sample_weight(sample_weight, X)

        if issparse(X):
            # Pre-sort indices to avoid that each individual tree of the
            # ensemble sorts the indices.
            X.sort_indices()

        y = np.atleast_1d(y)
        if y.ndim == 2 and y.shape[1] == 1:
            warn(
                (
                    "A column-vector y was passed when a 1d array was"
                    " expected. Please change the shape of y to "
                    "(n_samples,), for example using ravel()."
                ),
                DataConversionWarning,
                stacklevel=2,
            )

        if y.ndim == 1:
            # reshape is necessary to preserve the data contiguity against vs
            # [:, np.newaxis] that does not.
            y = np.reshape(y, (-1, 1))

        if self.criterion == "poisson":
            if np.any(y < 0):
                raise ValueError(
                    "Some value(s) of y are negative which is "
                    "not allowed for Poisson regression."
                )
            if np.sum(y) <= 0:
                raise ValueError(
                    "Sum of y is not strictly positive which "
                    "is necessary for Poisson regression."
                )

        self.n_outputs_ = y.shape[1]

        y, expanded_class_weight = self._validate_y_class_weight(y, classes=classes)

        if getattr(y, "dtype", None) != DOUBLE or not y.flags.contiguous:
            y = np.ascontiguousarray(y, dtype=DOUBLE)

        if expanded_class_weight is not None:
            if sample_weight is not None:
                sample_weight = sample_weight * expanded_class_weight
            else:
                sample_weight = expanded_class_weight

        if not self.bootstrap and self.max_samples is not None:
            raise ValueError(
                "`max_sample` cannot be set if `bootstrap=False`. "
                "Either switch to `bootstrap=True` or set "
                "`max_sample=None`."
            )
        elif self.bootstrap:
            n_samples_bootstrap = _get_n_samples_bootstrap(
                n_samples=X.shape[0], max_samples=self.max_samples
            )
        else:
            n_samples_bootstrap = None

        self._validate_estimator()

        if not self.bootstrap and self.oob_score:
            raise ValueError("Out of bag estimation only available if bootstrap=True")

        random_state = check_random_state(self.random_state)

        if not self.warm_start or not hasattr(self, "estimators_"):
            # Free allocated memory, if any
            self.estimators_ = []

        n_more_estimators = self.n_estimators - len(self.estimators_)

        if self.max_bins is not None:
            # `_openmp_effective_n_threads` is used to take cgroups CPU quotes
            # into account when determine the maximum number of threads to use.
            n_threads = _openmp_effective_n_threads()

            # Bin the data
            # For ease of use of the API, the user-facing GBDT classes accept the
            # parameter max_bins, which doesn't take into account the bin for
            # missing values (which is always allocated). However, since max_bins
            # isn't the true maximal number of bins, all other private classes
            # (binmapper, histbuilder...) accept n_bins instead, which is the
            # actual total number of bins. Everywhere in the code, the
            # convention is that n_bins == max_bins + 1
            n_bins = self.max_bins + 1  # + 1 for missing values
            self._bin_mapper = _BinMapper(
                n_bins=n_bins,
                # is_categorical=self.is_categorical_,
                known_categories=None,
                random_state=random_state,
                n_threads=n_threads,
            )

            # XXX: in order for this to work with the underlying tree submodule's Cython
            # code, we need to convert this into the original data's DTYPE because
            # the Cython code assumes that `DTYPE` is used.
            # The proper implementation will be a lot more complicated and should be
            # tackled once scikit-learn has finalized their inclusion of missing data
            # and categorical support for decision trees
            X = self._bin_data(X, is_training_data=True)  # .astype(DTYPE)
        else:
            self._bin_mapper = None

        if n_more_estimators < 0:
            raise ValueError(
                "n_estimators=%d must be larger or equal to "
                "len(estimators_)=%d when warm_start==True"
                % (self.n_estimators, len(self.estimators_))
            )

        elif n_more_estimators == 0:
            warn(
                "Warm-start fitting without increasing n_estimators does not "
                "fit new trees."
            )
        else:
            if self.warm_start and len(self.estimators_) > 0:
                # We draw from the random state to get the random state we
                # would have got if we hadn't used a warm_start.
                random_state.randint(MAX_INT, size=len(self.estimators_))

            trees = [
                self._make_estimator(append=False, random_state=random_state)
                for i in range(n_more_estimators)
            ]

            # Parallel loop: we prefer the threading backend as the Cython code
            # for fitting the trees is internally releasing the Python GIL
            # making threading more efficient than multiprocessing in
            # that case. However, for joblib 0.12+ we respect any
            # parallel_backend contexts set at a higher level,
            # since correctness does not rely on using threads.
            trees = Parallel(
                n_jobs=self.n_jobs,
                verbose=self.verbose,
                prefer="threads",
            )(
                delayed(_parallel_build_trees)(
                    t,
                    self.bootstrap,
                    X,
                    y,
                    sample_weight,
                    i,
                    len(trees),
                    verbose=self.verbose,
                    class_weight=self.class_weight,
                    n_samples_bootstrap=n_samples_bootstrap,
                    missing_values_in_feature_mask=missing_values_in_feature_mask,
                    classes=classes,
                )
                for i, t in enumerate(trees)
            )

            # Collect newly grown trees
            self.estimators_.extend(trees)

        if self.oob_score and (
            n_more_estimators > 0 or not hasattr(self, "oob_score_")
        ):
            y_type = type_of_target(y)
            if y_type in ("multiclass-multioutput", "unknown"):
                # FIXME: we could consider to support multiclass-multioutput if
                # we introduce or reuse a constructor parameter (e.g.
                # oob_score) allowing our user to pass a callable defining the
                # scoring strategy on OOB sample.
                raise ValueError(
                    "The type of target cannot be used to compute OOB "
                    f"estimates. Got {y_type} while only the following are "
                    "supported: continuous, continuous-multioutput, binary, "
                    "multiclass, multilabel-indicator."
                )

            if callable(self.oob_score):
                self._set_oob_score_and_attributes(
                    X, y, scoring_function=self.oob_score
                )
            else:
                self._set_oob_score_and_attributes(X, y)

        # Decapsulate classes_ attributes
        if hasattr(self, "classes_") and self.n_outputs_ == 1:
            self.n_classes_ = self.n_classes_[0]
            self.classes_ = self.classes_[0]

        return self

    @abstractmethod
    def _set_oob_score_and_attributes(self, X, y, scoring_function=None):
        """Compute and set the OOB score and attributes.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data matrix.
        y : ndarray of shape (n_samples, n_outputs)
            The target matrix.
        scoring_function : callable, default=None
            Scoring function for OOB score. Default depends on whether
            this is a regression (R2 score) or classification problem
            (accuracy score).
        """

    def _compute_oob_predictions(self, X, y):
        """Compute and set the OOB score.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data matrix.
        y : ndarray of shape (n_samples, n_outputs)
            The target matrix.

        Returns
        -------
        oob_pred : ndarray of shape (n_samples, n_classes, n_outputs) or \
                (n_samples, 1, n_outputs)
            The OOB predictions.
        """
        # Prediction requires X to be in CSR format
        if issparse(X):
            X = X.tocsr()

        n_samples = y.shape[0]
        n_outputs = self.n_outputs_
        if is_classifier(self) and hasattr(self, "n_classes_"):
            # n_classes_ is a ndarray at this stage
            # all the supported type of target will have the same number of
            # classes in all outputs
            oob_pred_shape = (n_samples, self.n_classes_[0], n_outputs)
        else:
            # for regression, n_classes_ does not exist and we create an empty
            # axis to be consistent with the classification case and make
            # the array operations compatible with the 2 settings
            oob_pred_shape = (n_samples, 1, n_outputs)

        oob_pred = np.zeros(shape=oob_pred_shape, dtype=np.float64)
        n_oob_pred = np.zeros((n_samples, n_outputs), dtype=np.int64)

        n_samples_bootstrap = _get_n_samples_bootstrap(
            n_samples,
            self.max_samples,
        )
        for estimator in self.estimators_:
            unsampled_indices = _generate_unsampled_indices(
                estimator.random_state,
                n_samples,
                n_samples_bootstrap,
            )

            y_pred = self._get_oob_predictions(estimator, X[unsampled_indices, :])
            oob_pred[unsampled_indices, ...] += y_pred
            n_oob_pred[unsampled_indices, :] += 1

        for k in range(n_outputs):
            if (n_oob_pred == 0).any():
                warn(
                    (
                        "Some inputs do not have OOB scores. This probably means "
                        "too few trees were used to compute any reliable OOB "
                        "estimates."
                    ),
                    UserWarning,
                )
                n_oob_pred[n_oob_pred == 0] = 1
            oob_pred[..., k] /= n_oob_pred[..., [k]]

        return oob_pred

    def _validate_y_class_weight(self, y, classes=None):
        # Default implementation
        return y, None

    def _validate_X_predict(self, X):
        """
        Validate X whenever one tries to predict, apply, predict_proba."""
        check_is_fitted(self)
        if self.estimators_[0]._support_missing_values(X):
            force_all_finite = "allow-nan"
        else:
            force_all_finite = True

        X = self._validate_data(
            X,
            dtype=DTYPE,
            accept_sparse="csr",
            reset=False,
            force_all_finite=force_all_finite,
        )
        if issparse(X) and (X.indices.dtype != np.intc or X.indptr.dtype != np.intc):
            raise ValueError("No support for np.int64 index based sparse matrices")
        return X

    @property
    def feature_importances_(self):
        """
        The impurity-based feature importances.

        The higher, the more important the feature.
        The importance of a feature is computed as the (normalized)
        total reduction of the criterion brought by that feature.  It is also
        known as the Gini importance.

        Warning: impurity-based feature importances can be misleading for
        high cardinality features (many unique values). See
        :func:`sklearn.inspection.permutation_importance` as an alternative.

        Returns
        -------
        feature_importances_ : ndarray of shape (n_features,)
            The values of this array sum to 1, unless all trees are single node
            trees consisting of only the root node, in which case it will be an
            array of zeros.
        """
        check_is_fitted(self)

        all_importances = Parallel(n_jobs=self.n_jobs, prefer="threads")(
            delayed(getattr)(tree, "feature_importances_")
            for tree in self.estimators_
            if tree.tree_.node_count > 1
        )

        if not all_importances:
            return np.zeros(self.n_features_in_, dtype=np.float64)

        all_importances = np.mean(all_importances, axis=0, dtype=np.float64)
        return all_importances / np.sum(all_importances)

    def _bin_data(self, X, is_training_data):
        """Bin data X.

        If is_training_data, then fit the _bin_mapper attribute.
        Else, the binned data is converted to a C-contiguous array.
        """
        description = "training" if is_training_data else "validation"
        if self.verbose:
            print(
                "Binning {:.3f} GB of {} data: ".format(X.nbytes / 1e9, description),
                end="",
                flush=True,
            )
        tic = time()
        if is_training_data:
            X_binned = self._bin_mapper.fit_transform(X)  # F-aligned array
        else:
            X_binned = self._bin_mapper.transform(X)  # F-aligned array
            # We convert the array to C-contiguous since predicting is faster
            # with this layout (training is faster on F-arrays though)
            X_binned = np.ascontiguousarray(X_binned)
        toc = time()
        if self.verbose:
            duration = toc - tic
            print("{:.3f} s".format(duration))

        return X_binned

    def predict_quantiles(self, X, quantiles=0.5, method="nearest"):
        """Predict class or regression value for X at given quantiles.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Input data.
        quantiles : float, optional
            The quantiles at which to evaluate, by default 0.5 (median).
        method : str, optional
            The method to interpolate, by default 'linear'. Can be any keyword
            argument accepted by :func:`~np.quantile`.

        Returns
        -------
        y : ndarray of shape (n_samples, n_quantiles, [n_outputs])
            The predicted values. The ``n_outputs`` dimension is present only
            for multi-output regressors.
        """
        if not self.store_leaf_values:
            raise RuntimeError(
                "Quantile prediction is not available when store_leaf_values=False"
            )
        check_is_fitted(self)
        # Check data
        X = self._validate_X_predict(X)

        if not isinstance(quantiles, (np.ndarray, list)):
            quantiles = np.array([quantiles])

        # if we trained a binning tree, then we should re-bin the data
        # XXX: this is inefficient and should be improved to be in line with what
        # the Histogram Gradient Boosting Tree does, where the binning thresholds
        # are passed into the tree itself, thus allowing us to set the node feature
        # value thresholds within the tree itself.
        if self.max_bins is not None:
            X = self._bin_data(X, is_training_data=False).astype(DTYPE)

        # Assign chunk of trees to jobs
        # n_jobs, _, _ = _partition_estimators(self.n_estimators, self.n_jobs)

        # avoid storing the output of every estimator by summing them here
        if self.n_outputs_ > 1:
            y_hat = np.zeros(
                (X.shape[0], len(quantiles), self.n_outputs_), dtype=np.float64
            )
        else:
            y_hat = np.zeros((X.shape[0], len(quantiles)), dtype=np.float64)

        # get (n_samples, n_estimators) indicator of leaf nodes
        X_leaves = self.apply(X)

        # we now want to aggregate all leaf samples across all trees for each sample
        for idx in range(X.shape[0]):
            # get leaf nodes for this sample
            leaf_nodes = X_leaves[idx, :]

            # (n_total_leaf_samples, n_outputs)
            leaf_node_samples = np.vstack(
                [
                    est.tree_.leaf_nodes_samples[leaf_nodes[jdx]]
                    for jdx, est in enumerate(self.estimators_)
                ]
            )

            # get quantiles across all leaf node samples
            try:
                y_hat[idx, ...] = np.quantile(
                    leaf_node_samples, quantiles, axis=0, method=method
                )
            except TypeError:
                y_hat[idx, ...] = np.quantile(
                    leaf_node_samples, quantiles, axis=0, interpolation=method
                )

            if is_classifier(self):
                if self.n_outputs_ == 1:
                    for i in range(len(quantiles)):
                        class_pred_per_sample = y_hat[idx, i, :].squeeze().astype(int)
                        y_hat[idx, ...] = self.classes_.take(
                            class_pred_per_sample, axis=0
                        )
                else:
                    for k in range(self.n_outputs_):
                        for i in range(len(quantiles)):
                            class_pred_per_sample = (
                                y_hat[idx, i, k].squeeze().astype(int)
                            )
                            y_hat[idx, i, k] = self.classes_[k].take(
                                class_pred_per_sample, axis=0
                            )
        return y_hat

    def get_leaf_node_samples(self, X):
        """For each datapoint x in X, get the training samples in the leaf node.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Dataset to apply the forest to.

        Returns
        -------
        leaf_node_samples : a list of array-like
            Each sample is represented by the indices of the training samples that
            reached the leaf node. The ``n_leaf_node_samples`` may vary between
            samples, since the number of samples that fall in a leaf node is
            variable. Each array-like has shape (n_leaf_node_samples, n_outputs).
        """
        if not self.store_leaf_values:
            raise RuntimeError(
                "Leaf node samples are not available when store_leaf_values=False"
            )

        check_is_fitted(self)
        # Check data
        X = self._validate_X_predict(X)

        # if we trained a binning tree, then we should re-bin the data
        # XXX: this is inefficient and should be improved to be in line with what
        # the Histogram Gradient Boosting Tree does, where the binning thresholds
        # are passed into the tree itself, thus allowing us to set the node feature
        # value thresholds within the tree itself.
        if self.max_bins is not None:
            X = self._bin_data(X, is_training_data=False).astype(DTYPE)

        # Assign chunk of trees to jobs
        n_jobs, _, _ = _partition_estimators(self.n_estimators, self.n_jobs)

        # avoid storing the output of every estimator by summing them here
        result = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
            delayed(_accumulate_leaf_nodes_samples)(e.get_leaf_node_samples, X)
            for e in self.estimators_
        )
        leaf_nodes_samples = result[0]
        for result_ in result[1:]:
            for i, node_samples in enumerate(result_):
                leaf_nodes_samples[i] = np.vstack((leaf_nodes_samples[i], node_samples))
        return leaf_nodes_samples

    def _more_tags(self):
        # Only the criterion is required to determine if the tree supports
        # missing values
        estimator = type(self.estimator)(criterion=self.criterion)
        return {"allow_nan": _safe_tags(estimator, key="allow_nan")}


def _accumulate_prediction(predict, X, out, lock):
    """
    This is a utility function for joblib's Parallel.

    It can't go locally in ForestClassifier or ForestRegressor, because joblib
    complains that it cannot pickle it when placed there.
    """
    prediction = predict(X, check_input=False)
    with lock:
        if len(out) == 1:
            out[0] += prediction
        else:
            for i in range(len(out)):
                out[i] += prediction[i]


def _accumulate_leaf_nodes_samples(func, X):
    """
    This is a utility function for joblib's Parallel.

    It can't go locally in ForestClassifier or ForestRegressor, because joblib
    complains that it cannot pickle it when placed there.
    """
    leaf_nodes_samples = func(X, check_input=False)
    return leaf_nodes_samples


class ForestClassifier(ClassifierMixin, BaseForest, metaclass=ABCMeta):
    """
    Base class for forest of trees-based classifiers.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """

    @abstractmethod
    def __init__(
        self,
        estimator,
        n_estimators=100,
        *,
        estimator_params=tuple(),
        bootstrap=False,
        oob_score=False,
        n_jobs=None,
        random_state=None,
        verbose=0,
        warm_start=False,
        class_weight=None,
        max_samples=None,
        base_estimator="deprecated",
        max_bins=None,
        store_leaf_values=False,
    ):
        super().__init__(
            estimator=estimator,
            n_estimators=n_estimators,
            estimator_params=estimator_params,
            bootstrap=bootstrap,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose,
            warm_start=warm_start,
            class_weight=class_weight,
            max_samples=max_samples,
            base_estimator=base_estimator,
            max_bins=max_bins,
            store_leaf_values=store_leaf_values,
        )

    @staticmethod
    def _get_oob_predictions(tree, X):
        """Compute the OOB predictions for an individual tree.

        Parameters
        ----------
        tree : DecisionTreeClassifier object
            A single decision tree classifier.
        X : ndarray of shape (n_samples, n_features)
            The OOB samples.

        Returns
        -------
        y_pred : ndarray of shape (n_samples, n_classes, n_outputs)
            The OOB associated predictions.
        """
        y_pred = tree.predict_proba(X, check_input=False)
        y_pred = np.array(y_pred, copy=False)
        if y_pred.ndim == 2:
            # binary and multiclass
            y_pred = y_pred[..., np.newaxis]
        else:
            # Roll the first `n_outputs` axis to the last axis. We will reshape
            # from a shape of (n_outputs, n_samples, n_classes) to a shape of
            # (n_samples, n_classes, n_outputs).
            y_pred = np.rollaxis(y_pred, axis=0, start=3)
        return y_pred

    def _set_oob_score_and_attributes(self, X, y, scoring_function=None):
        """Compute and set the OOB score and attributes.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data matrix.
        y : ndarray of shape (n_samples, n_outputs)
            The target matrix.
        scoring_function : callable, default=None
            Scoring function for OOB score. Defaults to `accuracy_score`.
        """
        self.oob_decision_function_ = super()._compute_oob_predictions(X, y)
        if self.oob_decision_function_.shape[-1] == 1:
            # drop the n_outputs axis if there is a single output
            self.oob_decision_function_ = self.oob_decision_function_.squeeze(axis=-1)

        if scoring_function is None:
            scoring_function = accuracy_score

        self.oob_score_ = scoring_function(
            y, np.argmax(self.oob_decision_function_, axis=1)
        )

    def _validate_y_class_weight(self, y, classes=None):
        check_classification_targets(y)

        y = np.copy(y)
        expanded_class_weight = None

        if self.class_weight is not None:
            y_original = np.copy(y)

        self.classes_ = []
        self.n_classes_ = []

        y_store_unique_indices = np.zeros(y.shape, dtype=int)
        if classes is not None:
            classes = np.atleast_1d(classes)
            if classes.ndim == 1:
                classes = np.array([classes])

            for k in classes:
                self.classes_.append(np.array(k))
                self.n_classes_.append(np.array(k).shape[0])

            for i in range(y.shape[0]):
                for j in range(self.n_outputs_):
                    y_store_unique_indices[i, j] = np.where(
                        self.classes_[j] == y[i, j]
                    )[0][0]
        else:
            for k in range(self.n_outputs_):
                classes_k, y_store_unique_indices[:, k] = np.unique(
                    y[:, k], return_inverse=True
                )
                self.classes_.append(classes_k)
                self.n_classes_.append(classes_k.shape[0])

        y = y_store_unique_indices

        if self.class_weight is not None:
            valid_presets = ("balanced", "balanced_subsample")
            if isinstance(self.class_weight, str):
                if self.class_weight not in valid_presets:
                    raise ValueError(
                        "Valid presets for class_weight include "
                        '"balanced" and "balanced_subsample".'
                        'Given "%s".'
                        % self.class_weight
                    )
                if self.warm_start:
                    warn(
                        'class_weight presets "balanced" or '
                        '"balanced_subsample" are '
                        "not recommended for warm_start if the fitted data "
                        "differs from the full dataset. In order to use "
                        '"balanced" weights, use compute_class_weight '
                        '("balanced", classes, y). In place of y you can use '
                        "a large enough sample of the full training set "
                        "target to properly estimate the class frequency "
                        "distributions. Pass the resulting weights as the "
                        "class_weight parameter."
                    )

            if self.class_weight != "balanced_subsample" or not self.bootstrap:
                if self.class_weight == "balanced_subsample":
                    class_weight = "balanced"
                else:
                    class_weight = self.class_weight
                expanded_class_weight = compute_sample_weight(class_weight, y_original)

        return y, expanded_class_weight

    def partial_fit(self, X, y, sample_weight=None, classes=None):
        """Update a decision tree classifier from the training set (X, y).

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The training input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csc_matrix``.

        y : array-like of shape (n_samples,) or (n_samples, n_outputs)
            The target values (class labels) as integers or strings.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights. If None, then samples are equally weighted. Splits
            that would create child nodes with net zero or negative weight are
            ignored while searching for a split in each node. Splits are also
            ignored if they would result in any single class carrying a
            negative weight in either child node.

        classes : array-like of shape (n_classes,), default=None
            List of all the classes that can possibly appear in the y vector.
            Must be provided at the first call to partial_fit, can be omitted
            in subsequent calls.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        self._validate_params()

        # validate input parameters
        first_call = _check_partial_fit_first_call(self, classes=classes)

        # Fit if no tree exists yet
        if first_call:
            self.fit(
                X,
                y,
                sample_weight=sample_weight,
                classes=classes,
            )
            return self

        X, y = self._validate_data(
            X,
            y,
            multi_output=True,
            accept_sparse="csc",
            dtype=DTYPE,
            force_all_finite=False,
            reset=first_call,
        )

        if issparse(y):
            raise ValueError("sparse multilabel-indicator for y is not supported.")

        if sample_weight is not None:
            sample_weight = _check_sample_weight(sample_weight, X)

        if issparse(X):
            # Pre-sort indices to avoid that each individual tree of the
            # ensemble sorts the indices.
            X.sort_indices()

        y = np.atleast_1d(y)
        if y.ndim == 2 and y.shape[1] == 1:
            warn(
                (
                    "A column-vector y was passed when a 1d array was"
                    " expected. Please change the shape of y to "
                    "(n_samples,), for example using ravel()."
                ),
                DataConversionWarning,
                stacklevel=2,
            )

        if y.ndim == 1:
            # reshape is necessary to preserve the data contiguity against vs
            # [:, np.newaxis] that does not.
            y = np.reshape(y, (-1, 1))

        if self.criterion == "poisson":
            if np.any(y < 0):
                raise ValueError(
                    "Some value(s) of y are negative which is "
                    "not allowed for Poisson regression."
                )
            if np.sum(y) <= 0:
                raise ValueError(
                    "Sum of y is not strictly positive which "
                    "is necessary for Poisson regression."
                )

        self.n_outputs_ = y.shape[1]

        classes = self.classes_
        if self.n_outputs_ == 1:
            classes = [classes]

        y, expanded_class_weight = self._validate_y_class_weight(y, classes)

        if getattr(y, "dtype", None) != DOUBLE or not y.flags.contiguous:
            y = np.ascontiguousarray(y, dtype=DOUBLE)

        if expanded_class_weight is not None:
            if sample_weight is not None:
                sample_weight = sample_weight * expanded_class_weight
            else:
                sample_weight = expanded_class_weight

        if not self.bootstrap and self.max_samples is not None:
            raise ValueError(
                "`max_sample` cannot be set if `bootstrap=False`. "
                "Either switch to `bootstrap=True` or set "
                "`max_sample=None`."
            )
        elif self.bootstrap:
            n_samples_bootstrap = _get_n_samples_bootstrap(
                n_samples=X.shape[0], max_samples=self.max_samples
            )
        else:
            n_samples_bootstrap = None

        self._validate_estimator()

        if not self.bootstrap and self.oob_score:
            raise ValueError("Out of bag estimation only available if bootstrap=True")

        random_state = check_random_state(self.random_state)

        if self.max_bins is not None:
            # `_openmp_effective_n_threads` is used to take cgroups CPU quotes
            # into account when determine the maximum number of threads to use.
            n_threads = _openmp_effective_n_threads()

            # Bin the data
            # For ease of use of the API, the user-facing GBDT classes accept the
            # parameter max_bins, which doesn't take into account the bin for
            # missing values (which is always allocated). However, since max_bins
            # isn't the true maximal number of bins, all other private classes
            # (binmapper, histbuilder...) accept n_bins instead, which is the
            # actual total number of bins. Everywhere in the code, the
            # convention is that n_bins == max_bins + 1
            n_bins = self.max_bins + 1  # + 1 for missing values
            self._bin_mapper = _BinMapper(
                n_bins=n_bins,
                # is_categorical=self.is_categorical_,
                known_categories=None,
                random_state=random_state,
                n_threads=n_threads,
            )

            # XXX: in order for this to work with the underlying tree submodule's Cython
            # code, we need to convert this into the original data's DTYPE because
            # the Cython code assumes that `DTYPE` is used.
            # The proper implementation will be a lot more complicated and should be
            # tackled once scikit-learn has finalized their inclusion of missing data
            # and categorical support for decision trees
            X = self._bin_data(X, is_training_data=True)  # .astype(DTYPE)
        else:
            self._bin_mapper = None

        # We draw from the random state to get the random state we
        # would have got if we hadn't used a warm_start.
        random_state.randint(MAX_INT, size=len(self.estimators_))

        # Parallel loop: we prefer the threading backend as the Cython code
        # for fitting the trees is internally releasing the Python GIL
        # making threading more efficient than multiprocessing in
        # that case. However, for joblib 0.12+ we respect any
        # parallel_backend contexts set at a higher level,
        # since correctness does not rely on using threads.
        Parallel(
            n_jobs=self.n_jobs,
            verbose=self.verbose,
            prefer="threads",
        )(
            delayed(_parallel_update_trees)(
                t,
                self.bootstrap,
                X,
                y,
                sample_weight,
                i,
                len(self.estimators_),
                verbose=self.verbose,
                class_weight=self.class_weight,
                n_samples_bootstrap=n_samples_bootstrap,
                classes=classes[0],
            )
            for i, t in enumerate(self.estimators_)
        )

        if self.oob_score:
            y_type = type_of_target(y)
            if y_type in ("multiclass-multioutput", "unknown"):
                # FIXME: we could consider to support multiclass-multioutput if
                # we introduce or reuse a constructor parameter (e.g.
                # oob_score) allowing our user to pass a callable defining the
                # scoring strategy on OOB sample.
                raise ValueError(
                    "The type of target cannot be used to compute OOB "
                    f"estimates. Got {y_type} while only the following are "
                    "supported: continuous, continuous-multioutput, binary, "
                    "multiclass, multilabel-indicator."
                )

            if callable(self.oob_score):
                self._set_oob_score_and_attributes(
                    X, y, scoring_function=self.oob_score
                )
            else:
                self._set_oob_score_and_attributes(X, y)

        # Decapsulate classes_ attributes
        if hasattr(self, "classes_") and self.n_outputs_ == 1:
            self.n_classes_ = self.n_classes_[0]
            self.classes_ = self.classes_[0]
        return self

    def predict(self, X):
        """
        Predict class for X.

        The predicted class of an input sample is a vote by the trees in
        the forest, weighted by their probability estimates. That is,
        the predicted class is the one with highest mean probability
        estimate across the trees.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The input samples. Internally, its dtype will be converted to
            ``dtype=np.float32``. If a sparse matrix is provided, it will be
            converted into a sparse ``csr_matrix``.

        Returns
        -------
        y : ndarray of shape (n_samples,) or (n_samples, n_outputs)
            The predicted classes.
        """
        proba = self.predict_proba(X)

        if self.n_outputs_ == 1:
            return self.classes_.take(np.argmax(proba, axis=1), axis=0)

        else:
            n_samples = proba[0].shape[0]
            # all dtypes should be the same, so just take the first
            class_type = self.classes_[0].dtype
            predictions = np.empty((n_samples, self.n_outputs_), dtype=class_type)

            for k in range(self.n_outputs_):
                predictions[:, k] = self.classes_[k].take(
                    np.argmax(proba[k], axis=1), axis=0
                )

            return predictions

    def predict_proba(self, X):
        """
        Predict class probabilities for X.

        The predicted class probabilities of an input sample are computed as
        the mean predicted class probabilities of the trees in the forest.
        The class probability of a single tree is the fraction of samples of
        the same class in a leaf.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The input samples. Internally, its dtype will be converted to
            ``dtype=np.float32``. If a sparse matrix is provided, it will be
            converted into a sparse ``csr_matrix``.

        Returns
        -------
        p : ndarray of shape (n_samples, n_classes), or a list of such arrays
            The class probabilities of the input samples. The order of the
            classes corresponds to that in the attribute :term:`classes_`.
        """
        check_is_fitted(self)
        # Check data
        X = self._validate_X_predict(X)

        # if we trained a binning tree, then we should re-bin the data
        # XXX: this is inefficient and should be improved to be in line with what
        # the Histogram Gradient Boosting Tree does, where the binning thresholds
        # are passed into the tree itself, thus allowing us to set the node feature
        # value thresholds within the tree itself.
        if self.max_bins is not None:
            X = self._bin_data(X, is_training_data=False).astype(DTYPE)

        # Assign chunk of trees to jobs
        n_jobs, _, _ = _partition_estimators(self.n_estimators, self.n_jobs)

        # avoid storing the output of every estimator by summing them here
        all_proba = [
            np.zeros((X.shape[0], j), dtype=np.float64)
            for j in np.atleast_1d(self.n_classes_)
        ]
        lock = threading.Lock()
        Parallel(n_jobs=n_jobs, verbose=self.verbose, require="sharedmem")(
            delayed(_accumulate_prediction)(e.predict_proba, X, all_proba, lock)
            for e in self.estimators_
        )

        for proba in all_proba:
            proba /= len(self.estimators_)

        if len(all_proba) == 1:
            return all_proba[0]
        else:
            return all_proba

    def predict_log_proba(self, X):
        """
        Predict class log-probabilities for X.

        The predicted class log-probabilities of an input sample is computed as
        the log of the mean predicted class probabilities of the trees in the
        forest.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The input samples. Internally, its dtype will be converted to
            ``dtype=np.float32``. If a sparse matrix is provided, it will be
            converted into a sparse ``csr_matrix``.

        Returns
        -------
        p : ndarray of shape (n_samples, n_classes), or a list of such arrays
            The class probabilities of the input samples. The order of the
            classes corresponds to that in the attribute :term:`classes_`.
        """
        proba = self.predict_proba(X)

        if self.n_outputs_ == 1:
            return np.log(proba)

        else:
            for k in range(self.n_outputs_):
                proba[k] = np.log(proba[k])

            return proba

    def _more_tags(self):
        return {"multilabel": True}


class ForestRegressor(RegressorMixin, BaseForest, metaclass=ABCMeta):
    """
    Base class for forest of trees-based regressors.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """

    @abstractmethod
    def __init__(
        self,
        estimator,
        n_estimators=100,
        *,
        estimator_params=tuple(),
        bootstrap=False,
        oob_score=False,
        n_jobs=None,
        random_state=None,
        verbose=0,
        warm_start=False,
        max_samples=None,
        base_estimator="deprecated",
        max_bins=None,
        store_leaf_values=False,
    ):
        super().__init__(
            estimator,
            n_estimators=n_estimators,
            estimator_params=estimator_params,
            bootstrap=bootstrap,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose,
            warm_start=warm_start,
            max_samples=max_samples,
            base_estimator=base_estimator,
            max_bins=max_bins,
            store_leaf_values=store_leaf_values,
        )

    def predict(self, X):
        """
        Predict regression target for X.

        The predicted regression target of an input sample is computed as the
        mean predicted regression targets of the trees in the forest.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The input samples. Internally, its dtype will be converted to
            ``dtype=np.float32``. If a sparse matrix is provided, it will be
            converted into a sparse ``csr_matrix``.

        Returns
        -------
        y : ndarray of shape (n_samples,) or (n_samples, n_outputs)
            The predicted values.
        """
        check_is_fitted(self)
        # Check data
        X = self._validate_X_predict(X)

        # if we trained a binning tree, then we should re-bin the data
        # XXX: this is inefficient and should be improved to be in line with what
        # the Histogram Gradient Boosting Tree does, where the binning thresholds
        # are passed into the tree itself, thus allowing us to set the node feature
        # value thresholds within the tree itself.
        if self.max_bins is not None:
            X = self._bin_data(X, is_training_data=False).astype(DTYPE)

        # Assign chunk of trees to jobs
        n_jobs, _, _ = _partition_estimators(self.n_estimators, self.n_jobs)

        # avoid storing the output of every estimator by summing them here
        if self.n_outputs_ > 1:
            y_hat = np.zeros((X.shape[0], self.n_outputs_), dtype=np.float64)
        else:
            y_hat = np.zeros((X.shape[0]), dtype=np.float64)

        # Parallel loop
        lock = threading.Lock()
        Parallel(n_jobs=n_jobs, verbose=self.verbose, require="sharedmem")(
            delayed(_accumulate_prediction)(e.predict, X, [y_hat], lock)
            for e in self.estimators_
        )

        y_hat /= len(self.estimators_)

        return y_hat

    @staticmethod
    def _get_oob_predictions(tree, X):
        """Compute the OOB predictions for an individual tree.

        Parameters
        ----------
        tree : DecisionTreeRegressor object
            A single decision tree regressor.
        X : ndarray of shape (n_samples, n_features)
            The OOB samples.

        Returns
        -------
        y_pred : ndarray of shape (n_samples, 1, n_outputs)
            The OOB associated predictions.
        """
        y_pred = tree.predict(X, check_input=False)
        if y_pred.ndim == 1:
            # single output regression
            y_pred = y_pred[:, np.newaxis, np.newaxis]
        else:
            # multioutput regression
            y_pred = y_pred[:, np.newaxis, :]
        return y_pred

    def _set_oob_score_and_attributes(self, X, y, scoring_function=None):
        """Compute and set the OOB score and attributes.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            The data matrix.
        y : ndarray of shape (n_samples, n_outputs)
            The target matrix.
        scoring_function : callable, default=None
            Scoring function for OOB score. Defaults to `r2_score`.
        """
        self.oob_prediction_ = super()._compute_oob_predictions(X, y).squeeze(axis=1)
        if self.oob_prediction_.shape[-1] == 1:
            # drop the n_outputs axis if there is a single output
            self.oob_prediction_ = self.oob_prediction_.squeeze(axis=-1)

        if scoring_function is None:
            scoring_function = r2_score

        self.oob_score_ = scoring_function(y, self.oob_prediction_)

    def _compute_partial_dependence_recursion(self, grid, target_features):
        """Fast partial dependence computation.

        Parameters
        ----------
        grid : ndarray of shape (n_samples, n_target_features)
            The grid points on which the partial dependence should be
            evaluated.
        target_features : ndarray of shape (n_target_features)
            The set of target features for which the partial dependence
            should be evaluated.

        Returns
        -------
        averaged_predictions : ndarray of shape (n_samples,)
            The value of the partial dependence function on each grid point.
        """
        grid = np.asarray(grid, dtype=DTYPE, order="C")
        averaged_predictions = np.zeros(
            shape=grid.shape[0], dtype=np.float64, order="C"
        )

        for tree in self.estimators_:
            # Note: we don't sum in parallel because the GIL isn't released in
            # the fast method.
            tree.tree_.compute_partial_dependence(
                grid, target_features, averaged_predictions
            )
        # Average over the forest
        averaged_predictions /= len(self.estimators_)

        return averaged_predictions

    def _more_tags(self):
        return {"multilabel": True}


class RandomForestClassifier(ForestClassifier):
    """
    A random forest classifier.

    A random forest is a meta estimator that fits a number of decision tree
    classifiers on various sub-samples of the dataset and uses averaging to
    improve the predictive accuracy and control over-fitting.
    The sub-sample size is controlled with the `max_samples` parameter if
    `bootstrap=True` (default), otherwise the whole dataset is used to build
    each tree.

    For a comparison between tree-based ensemble models see the example
    :ref:`sphx_glr_auto_examples_ensemble_plot_forest_hist_grad_boosting_comparison.py`.

    Read more in the :ref:`User Guide <forest>`.

    Parameters
    ----------
    n_estimators : int, default=100
        The number of trees in the forest.

        .. versionchanged:: 0.22
           The default value of ``n_estimators`` changed from 10 to 100
           in 0.22.

    criterion : {"gini", "entropy", "log_loss"}, default="gini"
        The function to measure the quality of a split. Supported criteria are
        "gini" for the Gini impurity and "log_loss" and "entropy" both for the
        Shannon information gain, see :ref:`tree_mathematical_formulation`.
        Note: This parameter is tree-specific.

    max_depth : int, default=None
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than
        min_samples_split samples.

    min_samples_split : int or float, default=2
        The minimum number of samples required to split an internal node:

        - If int, then consider `min_samples_split` as the minimum number.
        - If float, then `min_samples_split` is a fraction and
          `ceil(min_samples_split * n_samples)` are the minimum
          number of samples for each split.

        .. versionchanged:: 0.18
           Added float values for fractions.

    min_samples_leaf : int or float, default=1
        The minimum number of samples required to be at a leaf node.
        A split point at any depth will only be considered if it leaves at
        least ``min_samples_leaf`` training samples in each of the left and
        right branches.  This may have the effect of smoothing the model,
        especially in regression.

        - If int, then consider `min_samples_leaf` as the minimum number.
        - If float, then `min_samples_leaf` is a fraction and
          `ceil(min_samples_leaf * n_samples)` are the minimum
          number of samples for each node.

        .. versionchanged:: 0.18
           Added float values for fractions.

    min_weight_fraction_leaf : float, default=0.0
        The minimum weighted fraction of the sum total of weights (of all
        the input samples) required to be at a leaf node. Samples have
        equal weight when sample_weight is not provided.

    max_features : {"sqrt", "log2", None}, int or float, default="sqrt"
        The number of features to consider when looking for the best split:

        - If int, then consider `max_features` features at each split.
        - If float, then `max_features` is a fraction and
          `max(1, int(max_features * n_features_in_))` features are considered at each
          split.
        - If "sqrt", then `max_features=sqrt(n_features)`.
        - If "log2", then `max_features=log2(n_features)`.
        - If None, then `max_features=n_features`.

        .. versionchanged:: 1.1
            The default of `max_features` changed from `"auto"` to `"sqrt"`.

        Note: the search for a split does not stop until at least one
        valid partition of the node samples is found, even if it requires to
        effectively inspect more than ``max_features`` features.

    max_leaf_nodes : int, default=None
        Grow trees with ``max_leaf_nodes`` in best-first fashion.
        Best nodes are defined as relative reduction in impurity.
        If None then unlimited number of leaf nodes.

    min_impurity_decrease : float, default=0.0
        A node will be split if this split induces a decrease of the impurity
        greater than or equal to this value.

        The weighted impurity decrease equation is the following::

            N_t / N * (impurity - N_t_R / N_t * right_impurity
                                - N_t_L / N_t * left_impurity)

        where ``N`` is the total number of samples, ``N_t`` is the number of
        samples at the current node, ``N_t_L`` is the number of samples in the
        left child, and ``N_t_R`` is the number of samples in the right child.

        ``N``, ``N_t``, ``N_t_R`` and ``N_t_L`` all refer to the weighted sum,
        if ``sample_weight`` is passed.

        .. versionadded:: 0.19

    bootstrap : bool, default=True
        Whether bootstrap samples are used when building trees. If False, the
        whole dataset is used to build each tree.

    oob_score : bool or callable, default=False
        Whether to use out-of-bag samples to estimate the generalization score.
        By default, :func:`~sklearn.metrics.accuracy_score` is used.
        Provide a callable with signature `metric(y_true, y_pred)` to use a
        custom metric. Only available if `bootstrap=True`.

    n_jobs : int, default=None
        The number of jobs to run in parallel. :meth:`fit`, :meth:`predict`,
        :meth:`decision_path` and :meth:`apply` are all parallelized over the
        trees. ``None`` means 1 unless in a :obj:`joblib.parallel_backend`
        context. ``-1`` means using all processors. See :term:`Glossary
        <n_jobs>` for more details.

    random_state : int, RandomState instance or None, default=None
        Controls both the randomness of the bootstrapping of the samples used
        when building trees (if ``bootstrap=True``) and the sampling of the
        features to consider when looking for the best split at each node
        (if ``max_features < n_features``).
        See :term:`Glossary <random_state>` for details.

    verbose : int, default=0
        Controls the verbosity when fitting and predicting.

    warm_start : bool, default=False
        When set to ``True``, reuse the solution of the previous call to fit
        and add more estimators to the ensemble, otherwise, just fit a whole
        new forest. See :term:`Glossary <warm_start>` and
        :ref:`gradient_boosting_warm_start` for details.

    class_weight : {"balanced", "balanced_subsample"}, dict or list of dicts, \
            default=None
        Weights associated with classes in the form ``{class_label: weight}``.
        If not given, all classes are supposed to have weight one. For
        multi-output problems, a list of dicts can be provided in the same
        order as the columns of y.

        Note that for multioutput (including multilabel) weights should be
        defined for each class of every column in its own dict. For example,
        for four-class multilabel classification weights should be
        [{0: 1, 1: 1}, {0: 1, 1: 5}, {0: 1, 1: 1}, {0: 1, 1: 1}] instead of
        [{1:1}, {2:5}, {3:1}, {4:1}].

        The "balanced" mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies in the input data
        as ``n_samples / (n_classes * np.bincount(y))``

        The "balanced_subsample" mode is the same as "balanced" except that
        weights are computed based on the bootstrap sample for every tree
        grown.

        For multi-output, the weights of each column of y will be multiplied.

        Note that these weights will be multiplied with sample_weight (passed
        through the fit method) if sample_weight is specified.

    ccp_alpha : non-negative float, default=0.0
        Complexity parameter used for Minimal Cost-Complexity Pruning. The
        subtree with the largest cost complexity that is smaller than
        ``ccp_alpha`` will be chosen. By default, no pruning is performed. See
        :ref:`minimal_cost_complexity_pruning` for details.

        .. versionadded:: 0.22

    max_samples : int or float, default=None
        If bootstrap is True, the number of samples to draw from X
        to train each base estimator.

        - If None (default), then draw `X.shape[0]` samples.
        - If int, then draw `max_samples` samples.
        - If float, then draw `max(round(n_samples * max_samples), 1)` samples. Thus,
          `max_samples` should be in the interval `(0.0, 1.0]`.

        .. versionadded:: 0.22

    max_bins : int, default=255
        The maximum number of bins to use for non-missing values.

        **This is an experimental feature**.

    store_leaf_values : bool, default=False
        Whether to store the leaf values in the ``get_leaf_node_samples`` function.

        **This is an experimental feature**.

    monotonic_cst : array-like of int of shape (n_features), default=None
        Indicates the monotonicity constraint to enforce on each feature.
          - 1: monotonic increase
          - 0: no constraint
          - -1: monotonic decrease

        If monotonic_cst is None, no constraints are applied.

        Monotonicity constraints are not supported for:
          - multiclass classifications (i.e. when `n_classes > 2`),
          - multioutput classifications (i.e. when `n_outputs_ > 1`),
          - classifications trained on data with missing values.

        The constraints hold over the probability of the positive class.

        Read more in the :ref:`User Guide <monotonic_cst_gbdt>`.

        .. versionadded:: 1.4

    Attributes
    ----------
    estimator_ : :class:`~sklearn.tree.DecisionTreeClassifier`
        The child estimator template used to create the collection of fitted
        sub-estimators.

        .. versionadded:: 1.2
           `base_estimator_` was renamed to `estimator_`.

    base_estimator_ : DecisionTreeClassifier
        The child estimator template used to create the collection of fitted
        sub-estimators.

        .. deprecated:: 1.2
            `base_estimator_` is deprecated and will be removed in 1.4.
            Use `estimator_` instead.

    estimators_ : list of DecisionTreeClassifier
        The collection of fitted sub-estimators.

    classes_ : ndarray of shape (n_classes,) or a list of such arrays
        The classes labels (single output problem), or a list of arrays of
        class labels (multi-output problem).

    n_classes_ : int or list
        The number of classes (single output problem), or a list containing the
        number of classes for each output (multi-output problem).

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    n_outputs_ : int
        The number of outputs when ``fit`` is performed.

    feature_importances_ : ndarray of shape (n_features,)
        The impurity-based feature importances.
        The higher, the more important the feature.
        The importance of a feature is computed as the (normalized)
        total reduction of the criterion brought by that feature.  It is also
        known as the Gini importance.

        Warning: impurity-based feature importances can be misleading for
        high cardinality features (many unique values). See
        :func:`sklearn.inspection.permutation_importance` as an alternative.

    oob_score_ : float
        Score of the training dataset obtained using an out-of-bag estimate.
        This attribute exists only when ``oob_score`` is True.

    oob_decision_function_ : ndarray of shape (n_samples, n_classes) or \
            (n_samples, n_classes, n_outputs)
        Decision function computed with out-of-bag estimate on the training
        set. If n_estimators is small it might be possible that a data point
        was never left out during the bootstrap. In this case,
        `oob_decision_function_` might contain NaN. This attribute exists
        only when ``oob_score`` is True.

    See Also
    --------
    sklearn.tree.DecisionTreeClassifier : A decision tree classifier.
    sklearn.ensemble.ExtraTreesClassifier : Ensemble of extremely randomized
        tree classifiers.
    sklearn.ensemble.HistGradientBoostingClassifier : A Histogram-based Gradient
        Boosting Classification Tree, very fast for big datasets (n_samples >=
        10_000).

    Notes
    -----
    The default values for the parameters controlling the size of the trees
    (e.g. ``max_depth``, ``min_samples_leaf``, etc.) lead to fully grown and
    unpruned trees which can potentially be very large on some data sets. To
    reduce memory consumption, the complexity and size of the trees should be
    controlled by setting those parameter values.

    The features are always randomly permuted at each split. Therefore,
    the best found split may vary, even with the same training data,
    ``max_features=n_features`` and ``bootstrap=False``, if the improvement
    of the criterion is identical for several splits enumerated during the
    search of the best split. To obtain a deterministic behaviour during
    fitting, ``random_state`` has to be fixed.

    References
    ----------
    .. [1] L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32, 2001.

    Examples
    --------
    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from sklearn.datasets import make_classification
    >>> X, y = make_classification(n_samples=1000, n_features=4,
    ...                            n_informative=2, n_redundant=0,
    ...                            random_state=0, shuffle=False)
    >>> clf = RandomForestClassifier(max_depth=2, random_state=0)
    >>> clf.fit(X, y)
    RandomForestClassifier(...)
    >>> print(clf.predict([[0, 0, 0, 0]]))
    [1]
    """

    _parameter_constraints: dict = {
        **ForestClassifier._parameter_constraints,
        **DecisionTreeClassifier._parameter_constraints,
        "class_weight": [
            StrOptions({"balanced_subsample", "balanced"}),
            dict,
            list,
            None,
        ],
    }
    _parameter_constraints.pop("splitter")

    def __init__(
        self,
        n_estimators=100,
        *,
        criterion="gini",
        max_depth=None,
        min_samples_split=2,
        min_samples_leaf=1,
        min_weight_fraction_leaf=0.0,
        max_features="sqrt",
        max_leaf_nodes=None,
        min_impurity_decrease=0.0,
        bootstrap=True,
        oob_score=False,
        n_jobs=None,
        random_state=None,
        verbose=0,
        warm_start=False,
        class_weight=None,
        ccp_alpha=0.0,
        max_samples=None,
        max_bins=None,
        store_leaf_values=False,
        monotonic_cst=None,
    ):
        super().__init__(
            estimator=DecisionTreeClassifier(),
            n_estimators=n_estimators,
            estimator_params=(
                "criterion",
                "max_depth",
                "min_samples_split",
                "min_samples_leaf",
                "min_weight_fraction_leaf",
                "max_features",
                "max_leaf_nodes",
                "min_impurity_decrease",
                "random_state",
                "ccp_alpha",
                "store_leaf_values",
                "monotonic_cst",
            ),
            bootstrap=bootstrap,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose,
            warm_start=warm_start,
            class_weight=class_weight,
            max_samples=max_samples,
            max_bins=max_bins,
            store_leaf_values=store_leaf_values,
        )

        self.criterion = criterion
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_fraction_leaf = min_weight_fraction_leaf
        self.max_features = max_features
        self.max_leaf_nodes = max_leaf_nodes
        self.min_impurity_decrease = min_impurity_decrease
        self.monotonic_cst = monotonic_cst
        self.ccp_alpha = ccp_alpha


class RandomForestRegressor(ForestRegressor):
    """
    A random forest regressor.

    A random forest is a meta estimator that fits a number of classifying
    decision trees on various sub-samples of the dataset and uses averaging
    to improve the predictive accuracy and control over-fitting.
    The sub-sample size is controlled with the `max_samples` parameter if
    `bootstrap=True` (default), otherwise the whole dataset is used to build
    each tree.

    For a comparison between tree-based ensemble models see the example
    :ref:`sphx_glr_auto_examples_ensemble_plot_forest_hist_grad_boosting_comparison.py`.

    Read more in the :ref:`User Guide <forest>`.

    Parameters
    ----------
    n_estimators : int, default=100
        The number of trees in the forest.

        .. versionchanged:: 0.22
           The default value of ``n_estimators`` changed from 10 to 100
           in 0.22.

    criterion : {"squared_error", "absolute_error", "friedman_mse", "poisson"}, \
            default="squared_error"
        The function to measure the quality of a split. Supported criteria
        are "squared_error" for the mean squared error, which is equal to
        variance reduction as feature selection criterion and minimizes the L2
        loss using the mean of each terminal node, "friedman_mse", which uses
        mean squared error with Friedman's improvement score for potential
        splits, "absolute_error" for the mean absolute error, which minimizes
        the L1 loss using the median of each terminal node, and "poisson" which
        uses reduction in Poisson deviance to find splits.
        Training using "absolute_error" is significantly slower
        than when using "squared_error".

        .. versionadded:: 0.18
           Mean Absolute Error (MAE) criterion.

        .. versionadded:: 1.0
           Poisson criterion.

    max_depth : int, default=None
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than
        min_samples_split samples.

    min_samples_split : int or float, default=2
        The minimum number of samples required to split an internal node:

        - If int, then consider `min_samples_split` as the minimum number.
        - If float, then `min_samples_split` is a fraction and
          `ceil(min_samples_split * n_samples)` are the minimum
          number of samples for each split.

        .. versionchanged:: 0.18
           Added float values for fractions.

    min_samples_leaf : int or float, default=1
        The minimum number of samples required to be at a leaf node.
        A split point at any depth will only be considered if it leaves at
        least ``min_samples_leaf`` training samples in each of the left and
        right branches.  This may have the effect of smoothing the model,
        especially in regression.

        - If int, then consider `min_samples_leaf` as the minimum number.
        - If float, then `min_samples_leaf` is a fraction and
          `ceil(min_samples_leaf * n_samples)` are the minimum
          number of samples for each node.

        .. versionchanged:: 0.18
           Added float values for fractions.

    min_weight_fraction_leaf : float, default=0.0
        The minimum weighted fraction of the sum total of weights (of all
        the input samples) required to be at a leaf node. Samples have
        equal weight when sample_weight is not provided.

    max_features : {"sqrt", "log2", None}, int or float, default=1.0
        The number of features to consider when looking for the best split:

        - If int, then consider `max_features` features at each split.
        - If float, then `max_features` is a fraction and
          `max(1, int(max_features * n_features_in_))` features are considered at each
          split.
        - If "sqrt", then `max_features=sqrt(n_features)`.
        - If "log2", then `max_features=log2(n_features)`.
        - If None or 1.0, then `max_features=n_features`.

        .. note::
            The default of 1.0 is equivalent to bagged trees and more
            randomness can be achieved by setting smaller values, e.g. 0.3.

        .. versionchanged:: 1.1
            The default of `max_features` changed from `"auto"` to 1.0.

        Note: the search for a split does not stop until at least one
        valid partition of the node samples is found, even if it requires to
        effectively inspect more than ``max_features`` features.

    max_leaf_nodes : int, default=None
        Grow trees with ``max_leaf_nodes`` in best-first fashion.
        Best nodes are defined as relative reduction in impurity.
        If None then unlimited number of leaf nodes.

    min_impurity_decrease : float, default=0.0
        A node will be split if this split induces a decrease of the impurity
        greater than or equal to this value.

        The weighted impurity decrease equation is the following::

            N_t / N * (impurity - N_t_R / N_t * right_impurity
                                - N_t_L / N_t * left_impurity)

        where ``N`` is the total number of samples, ``N_t`` is the number of
        samples at the current node, ``N_t_L`` is the number of samples in the
        left child, and ``N_t_R`` is the number of samples in the right child.

        ``N``, ``N_t``, ``N_t_R`` and ``N_t_L`` all refer to the weighted sum,
        if ``sample_weight`` is passed.

        .. versionadded:: 0.19

    bootstrap : bool, default=True
        Whether bootstrap samples are used when building trees. If False, the
        whole dataset is used to build each tree.

    oob_score : bool or callable, default=False
        Whether to use out-of-bag samples to estimate the generalization score.
        By default, :func:`~sklearn.metrics.r2_score` is used.
        Provide a callable with signature `metric(y_true, y_pred)` to use a
        custom metric. Only available if `bootstrap=True`.

    n_jobs : int, default=None
        The number of jobs to run in parallel. :meth:`fit`, :meth:`predict`,
        :meth:`decision_path` and :meth:`apply` are all parallelized over the
        trees. ``None`` means 1 unless in a :obj:`joblib.parallel_backend`
        context. ``-1`` means using all processors. See :term:`Glossary
        <n_jobs>` for more details.

    random_state : int, RandomState instance or None, default=None
        Controls both the randomness of the bootstrapping of the samples used
        when building trees (if ``bootstrap=True``) and the sampling of the
        features to consider when looking for the best split at each node
        (if ``max_features < n_features``).
        See :term:`Glossary <random_state>` for details.

    verbose : int, default=0
        Controls the verbosity when fitting and predicting.

    warm_start : bool, default=False
        When set to ``True``, reuse the solution of the previous call to fit
        and add more estimators to the ensemble, otherwise, just fit a whole
        new forest. See :term:`Glossary <warm_start>` and
        :ref:`gradient_boosting_warm_start` for details.

    ccp_alpha : non-negative float, default=0.0
        Complexity parameter used for Minimal Cost-Complexity Pruning. The
        subtree with the largest cost complexity that is smaller than
        ``ccp_alpha`` will be chosen. By default, no pruning is performed. See
        :ref:`minimal_cost_complexity_pruning` for details.

        .. versionadded:: 0.22

    max_samples : int or float, default=None
        If bootstrap is True, the number of samples to draw from X
        to train each base estimator.

        - If None (default), then draw `X.shape[0]` samples.
        - If int, then draw `max_samples` samples.
        - If float, then draw `max(round(n_samples * max_samples), 1)` samples. Thus,
          `max_samples` should be in the interval `(0.0, 1.0]`.

        .. versionadded:: 0.22

    max_bins : int, default=255
        The maximum number of bins to use for non-missing values. Used for
        speeding up training time.

        **This is an experimental feature**.

    store_leaf_values : bool, default=False
        Whether to store the leaf values in the ``get_leaf_node_samples`` function.

        **This is an experimental feature**.

    monotonic_cst : array-like of int of shape (n_features), default=None
        Indicates the monotonicity constraint to enforce on each feature.
          - 1: monotonically increasing
          - 0: no constraint
          - -1: monotonically decreasing

        If monotonic_cst is None, no constraints are applied.

        Monotonicity constraints are not supported for:
          - multioutput regressions (i.e. when `n_outputs_ > 1`),
          - regressions trained on data with missing values.

        Read more in the :ref:`User Guide <monotonic_cst_gbdt>`.

        .. versionadded:: 1.4

    Attributes
    ----------
    estimator_ : :class:`~sklearn.tree.DecisionTreeRegressor`
        The child estimator template used to create the collection of fitted
        sub-estimators.

        .. versionadded:: 1.2
           `base_estimator_` was renamed to `estimator_`.

    base_estimator_ : DecisionTreeRegressor
        The child estimator template used to create the collection of fitted
        sub-estimators.

        .. deprecated:: 1.2
            `base_estimator_` is deprecated and will be removed in 1.4.
            Use `estimator_` instead.

    estimators_ : list of DecisionTreeRegressor
        The collection of fitted sub-estimators.

    feature_importances_ : ndarray of shape (n_features,)
        The impurity-based feature importances.
        The higher, the more important the feature.
        The importance of a feature is computed as the (normalized)
        total reduction of the criterion brought by that feature.  It is also
        known as the Gini importance.

        Warning: impurity-based feature importances can be misleading for
        high cardinality features (many unique values). See
        :func:`sklearn.inspection.permutation_importance` as an alternative.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    n_outputs_ : int
        The number of outputs when ``fit`` is performed.

    oob_score_ : float
        Score of the training dataset obtained using an out-of-bag estimate.
        This attribute exists only when ``oob_score`` is True.

    oob_prediction_ : ndarray of shape (n_samples,) or (n_samples, n_outputs)
        Prediction computed with out-of-bag estimate on the training set.
        This attribute exists only when ``oob_score`` is True.

    See Also
    --------
    sklearn.tree.DecisionTreeRegressor : A decision tree regressor.
    sklearn.ensemble.ExtraTreesRegressor : Ensemble of extremely randomized
        tree regressors.
    sklearn.ensemble.HistGradientBoostingRegressor : A Histogram-based Gradient
        Boosting Regression Tree, very fast for big datasets (n_samples >=
        10_000).

    Notes
    -----
    The default values for the parameters controlling the size of the trees
    (e.g. ``max_depth``, ``min_samples_leaf``, etc.) lead to fully grown and
    unpruned trees which can potentially be very large on some data sets. To
    reduce memory consumption, the complexity and size of the trees should be
    controlled by setting those parameter values.

    The features are always randomly permuted at each split. Therefore,
    the best found split may vary, even with the same training data,
    ``max_features=n_features`` and ``bootstrap=False``, if the improvement
    of the criterion is identical for several splits enumerated during the
    search of the best split. To obtain a deterministic behaviour during
    fitting, ``random_state`` has to be fixed.

    The default value ``max_features=1.0`` uses ``n_features``
    rather than ``n_features / 3``. The latter was originally suggested in
    [1], whereas the former was more recently justified empirically in [2].

    References
    ----------
    .. [1] L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32, 2001.

    .. [2] P. Geurts, D. Ernst., and L. Wehenkel, "Extremely randomized
           trees", Machine Learning, 63(1), 3-42, 2006.

    Examples
    --------
    >>> from sklearn.ensemble import RandomForestRegressor
    >>> from sklearn.datasets import make_regression
    >>> X, y = make_regression(n_features=4, n_informative=2,
    ...                        random_state=0, shuffle=False)
    >>> regr = RandomForestRegressor(max_depth=2, random_state=0)
    >>> regr.fit(X, y)
    RandomForestRegressor(...)
    >>> print(regr.predict([[0, 0, 0, 0]]))
    [-8.32987858]
    """

    _parameter_constraints: dict = {
        **ForestRegressor._parameter_constraints,
        **DecisionTreeRegressor._parameter_constraints,
    }
    _parameter_constraints.pop("splitter")

    def __init__(
        self,
        n_estimators=100,
        *,
        criterion="squared_error",
        max_depth=None,
        min_samples_split=2,
        min_samples_leaf=1,
        min_weight_fraction_leaf=0.0,
        max_features=1.0,
        max_leaf_nodes=None,
        min_impurity_decrease=0.0,
        bootstrap=True,
        oob_score=False,
        n_jobs=None,
        random_state=None,
        verbose=0,
        warm_start=False,
        ccp_alpha=0.0,
        max_samples=None,
        max_bins=None,
        store_leaf_values=False,
        monotonic_cst=None,
    ):
        super().__init__(
            estimator=DecisionTreeRegressor(),
            n_estimators=n_estimators,
            estimator_params=(
                "criterion",
                "max_depth",
                "min_samples_split",
                "min_samples_leaf",
                "min_weight_fraction_leaf",
                "max_features",
                "max_leaf_nodes",
                "min_impurity_decrease",
                "random_state",
                "ccp_alpha",
                "store_leaf_values",
                "monotonic_cst",
            ),
            bootstrap=bootstrap,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose,
            warm_start=warm_start,
            max_samples=max_samples,
            max_bins=max_bins,
            store_leaf_values=store_leaf_values,
        )

        self.criterion = criterion
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_fraction_leaf = min_weight_fraction_leaf
        self.max_features = max_features
        self.max_leaf_nodes = max_leaf_nodes
        self.min_impurity_decrease = min_impurity_decrease
        self.ccp_alpha = ccp_alpha
        self.monotonic_cst = monotonic_cst


class ExtraTreesClassifier(ForestClassifier):
    """
    An extra-trees classifier.

    This class implements a meta estimator that fits a number of
    randomized decision trees (a.k.a. extra-trees) on various sub-samples
    of the dataset and uses averaging to improve the predictive accuracy
    and control over-fitting.

    Read more in the :ref:`User Guide <forest>`.

    Parameters
    ----------
    n_estimators : int, default=100
        The number of trees in the forest.

        .. versionchanged:: 0.22
           The default value of ``n_estimators`` changed from 10 to 100
           in 0.22.

    criterion : {"gini", "entropy", "log_loss"}, default="gini"
        The function to measure the quality of a split. Supported criteria are
        "gini" for the Gini impurity and "log_loss" and "entropy" both for the
        Shannon information gain, see :ref:`tree_mathematical_formulation`.
        Note: This parameter is tree-specific.

    max_depth : int, default=None
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than
        min_samples_split samples.

    min_samples_split : int or float, default=2
        The minimum number of samples required to split an internal node:

        - If int, then consider `min_samples_split` as the minimum number.
        - If float, then `min_samples_split` is a fraction and
          `ceil(min_samples_split * n_samples)` are the minimum
          number of samples for each split.

        .. versionchanged:: 0.18
           Added float values for fractions.

    min_samples_leaf : int or float, default=1
        The minimum number of samples required to be at a leaf node.
        A split point at any depth will only be considered if it leaves at
        least ``min_samples_leaf`` training samples in each of the left and
        right branches.  This may have the effect of smoothing the model,
        especially in regression.

        - If int, then consider `min_samples_leaf` as the minimum number.
        - If float, then `min_samples_leaf` is a fraction and
          `ceil(min_samples_leaf * n_samples)` are the minimum
          number of samples for each node.

        .. versionchanged:: 0.18
           Added float values for fractions.

    min_weight_fraction_leaf : float, default=0.0
        The minimum weighted fraction of the sum total of weights (of all
        the input samples) required to be at a leaf node. Samples have
        equal weight when sample_weight is not provided.

    max_features : {"sqrt", "log2", None}, int or float, default="sqrt"
        The number of features to consider when looking for the best split:

        - If int, then consider `max_features` features at each split.
        - If float, then `max_features` is a fraction and
          `max(1, int(max_features * n_features_in_))` features are considered at each
          split.
        - If "sqrt", then `max_features=sqrt(n_features)`.
        - If "log2", then `max_features=log2(n_features)`.
        - If None, then `max_features=n_features`.

        .. versionchanged:: 1.1
            The default of `max_features` changed from `"auto"` to `"sqrt"`.

        Note: the search for a split does not stop until at least one
        valid partition of the node samples is found, even if it requires to
        effectively inspect more than ``max_features`` features.

    max_leaf_nodes : int, default=None
        Grow trees with ``max_leaf_nodes`` in best-first fashion.
        Best nodes are defined as relative reduction in impurity.
        If None then unlimited number of leaf nodes.

    min_impurity_decrease : float, default=0.0
        A node will be split if this split induces a decrease of the impurity
        greater than or equal to this value.

        The weighted impurity decrease equation is the following::

            N_t / N * (impurity - N_t_R / N_t * right_impurity
                                - N_t_L / N_t * left_impurity)

        where ``N`` is the total number of samples, ``N_t`` is the number of
        samples at the current node, ``N_t_L`` is the number of samples in the
        left child, and ``N_t_R`` is the number of samples in the right child.

        ``N``, ``N_t``, ``N_t_R`` and ``N_t_L`` all refer to the weighted sum,
        if ``sample_weight`` is passed.

        .. versionadded:: 0.19

    bootstrap : bool, default=False
        Whether bootstrap samples are used when building trees. If False, the
        whole dataset is used to build each tree.

    oob_score : bool or callable, default=False
        Whether to use out-of-bag samples to estimate the generalization score.
        By default, :func:`~sklearn.metrics.accuracy_score` is used.
        Provide a callable with signature `metric(y_true, y_pred)` to use a
        custom metric. Only available if `bootstrap=True`.

    n_jobs : int, default=None
        The number of jobs to run in parallel. :meth:`fit`, :meth:`predict`,
        :meth:`decision_path` and :meth:`apply` are all parallelized over the
        trees. ``None`` means 1 unless in a :obj:`joblib.parallel_backend`
        context. ``-1`` means using all processors. See :term:`Glossary
        <n_jobs>` for more details.

    random_state : int, RandomState instance or None, default=None
        Controls 3 sources of randomness:

        - the bootstrapping of the samples used when building trees
          (if ``bootstrap=True``)
        - the sampling of the features to consider when looking for the best
          split at each node (if ``max_features < n_features``)
        - the draw of the splits for each of the `max_features`

        See :term:`Glossary <random_state>` for details.

    verbose : int, default=0
        Controls the verbosity when fitting and predicting.

    warm_start : bool, default=False
        When set to ``True``, reuse the solution of the previous call to fit
        and add more estimators to the ensemble, otherwise, just fit a whole
        new forest. See :term:`Glossary <warm_start>` and
        :ref:`gradient_boosting_warm_start` for details.

    class_weight : {"balanced", "balanced_subsample"}, dict or list of dicts, \
            default=None
        Weights associated with classes in the form ``{class_label: weight}``.
        If not given, all classes are supposed to have weight one. For
        multi-output problems, a list of dicts can be provided in the same
        order as the columns of y.

        Note that for multioutput (including multilabel) weights should be
        defined for each class of every column in its own dict. For example,
        for four-class multilabel classification weights should be
        [{0: 1, 1: 1}, {0: 1, 1: 5}, {0: 1, 1: 1}, {0: 1, 1: 1}] instead of
        [{1:1}, {2:5}, {3:1}, {4:1}].

        The "balanced" mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies in the input data
        as ``n_samples / (n_classes * np.bincount(y))``

        The "balanced_subsample" mode is the same as "balanced" except that
        weights are computed based on the bootstrap sample for every tree
        grown.

        For multi-output, the weights of each column of y will be multiplied.

        Note that these weights will be multiplied with sample_weight (passed
        through the fit method) if sample_weight is specified.

    ccp_alpha : non-negative float, default=0.0
        Complexity parameter used for Minimal Cost-Complexity Pruning. The
        subtree with the largest cost complexity that is smaller than
        ``ccp_alpha`` will be chosen. By default, no pruning is performed. See
        :ref:`minimal_cost_complexity_pruning` for details.

        .. versionadded:: 0.22

    max_samples : int or float, default=None
        If bootstrap is True, the number of samples to draw from X
        to train each base estimator.

        - If None (default), then draw `X.shape[0]` samples.
        - If int, then draw `max_samples` samples.
        - If float, then draw `max_samples * X.shape[0]` samples. Thus,
          `max_samples` should be in the interval `(0.0, 1.0]`.

        .. versionadded:: 0.22

    max_bins : int, default=255
        The maximum number of bins to use for non-missing values.

        **This is an experimental feature**.

    store_leaf_values : bool, default=False
        Whether to store the leaf values in the ``get_leaf_node_samples`` function.

        **This is an experimental feature**.

    monotonic_cst : array-like of int of shape (n_features), default=None
        Indicates the monotonicity constraint to enforce on each feature.
          - 1: monotonically increasing
          - 0: no constraint
          - -1: monotonically decreasing

        If monotonic_cst is None, no constraints are applied.

        Monotonicity constraints are not supported for:
          - multiclass classifications (i.e. when `n_classes > 2`),
          - multioutput classifications (i.e. when `n_outputs_ > 1`),
          - classifications trained on data with missing values.

        The constraints hold over the probability of the positive class.

        Read more in the :ref:`User Guide <monotonic_cst_gbdt>`.

        .. versionadded:: 1.4

    Attributes
    ----------
    estimator_ : :class:`~sklearn.tree.ExtraTreeClassifier`
        The child estimator template used to create the collection of fitted
        sub-estimators.

        .. versionadded:: 1.2
           `base_estimator_` was renamed to `estimator_`.

    base_estimator_ : ExtraTreesClassifier
        The child estimator template used to create the collection of fitted
        sub-estimators.

        .. deprecated:: 1.2
            `base_estimator_` is deprecated and will be removed in 1.4.
            Use `estimator_` instead.

    estimators_ : list of DecisionTreeClassifier
        The collection of fitted sub-estimators.

    classes_ : ndarray of shape (n_classes,) or a list of such arrays
        The classes labels (single output problem), or a list of arrays of
        class labels (multi-output problem).

    n_classes_ : int or list
        The number of classes (single output problem), or a list containing the
        number of classes for each output (multi-output problem).

    feature_importances_ : ndarray of shape (n_features,)
        The impurity-based feature importances.
        The higher, the more important the feature.
        The importance of a feature is computed as the (normalized)
        total reduction of the criterion brought by that feature.  It is also
        known as the Gini importance.

        Warning: impurity-based feature importances can be misleading for
        high cardinality features (many unique values). See
        :func:`sklearn.inspection.permutation_importance` as an alternative.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    n_outputs_ : int
        The number of outputs when ``fit`` is performed.

    oob_score_ : float
        Score of the training dataset obtained using an out-of-bag estimate.
        This attribute exists only when ``oob_score`` is True.

    oob_decision_function_ : ndarray of shape (n_samples, n_classes) or \
            (n_samples, n_classes, n_outputs)
        Decision function computed with out-of-bag estimate on the training
        set. If n_estimators is small it might be possible that a data point
        was never left out during the bootstrap. In this case,
        `oob_decision_function_` might contain NaN. This attribute exists
        only when ``oob_score`` is True.

    See Also
    --------
    ExtraTreesRegressor : An extra-trees regressor with random splits.
    RandomForestClassifier : A random forest classifier with optimal splits.
    RandomForestRegressor : Ensemble regressor using trees with optimal splits.

    Notes
    -----
    The default values for the parameters controlling the size of the trees
    (e.g. ``max_depth``, ``min_samples_leaf``, etc.) lead to fully grown and
    unpruned trees which can potentially be very large on some data sets. To
    reduce memory consumption, the complexity and size of the trees should be
    controlled by setting those parameter values.

    References
    ----------
    .. [1] P. Geurts, D. Ernst., and L. Wehenkel, "Extremely randomized
           trees", Machine Learning, 63(1), 3-42, 2006.

    Examples
    --------
    >>> from sklearn.ensemble import ExtraTreesClassifier
    >>> from sklearn.datasets import make_classification
    >>> X, y = make_classification(n_features=4, random_state=0)
    >>> clf = ExtraTreesClassifier(n_estimators=100, random_state=0)
    >>> clf.fit(X, y)
    ExtraTreesClassifier(random_state=0)
    >>> clf.predict([[0, 0, 0, 0]])
    array([1])
    """

    _parameter_constraints: dict = {
        **ForestClassifier._parameter_constraints,
        **DecisionTreeClassifier._parameter_constraints,
        "class_weight": [
            StrOptions({"balanced_subsample", "balanced"}),
            dict,
            list,
            None,
        ],
    }
    _parameter_constraints.pop("splitter")

    def __init__(
        self,
        n_estimators=100,
        *,
        criterion="gini",
        max_depth=None,
        min_samples_split=2,
        min_samples_leaf=1,
        min_weight_fraction_leaf=0.0,
        max_features="sqrt",
        max_leaf_nodes=None,
        min_impurity_decrease=0.0,
        bootstrap=False,
        oob_score=False,
        n_jobs=None,
        random_state=None,
        verbose=0,
        warm_start=False,
        class_weight=None,
        ccp_alpha=0.0,
        max_samples=None,
        max_bins=None,
        store_leaf_values=False,
        monotonic_cst=None,
    ):
        super().__init__(
            estimator=ExtraTreeClassifier(),
            n_estimators=n_estimators,
            estimator_params=(
                "criterion",
                "max_depth",
                "min_samples_split",
                "min_samples_leaf",
                "min_weight_fraction_leaf",
                "max_features",
                "max_leaf_nodes",
                "min_impurity_decrease",
                "random_state",
                "ccp_alpha",
                "store_leaf_values",
                "monotonic_cst",
            ),
            bootstrap=bootstrap,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose,
            warm_start=warm_start,
            class_weight=class_weight,
            max_samples=max_samples,
            max_bins=max_bins,
            store_leaf_values=store_leaf_values,
        )

        self.criterion = criterion
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_fraction_leaf = min_weight_fraction_leaf
        self.max_features = max_features
        self.max_leaf_nodes = max_leaf_nodes
        self.min_impurity_decrease = min_impurity_decrease
        self.ccp_alpha = ccp_alpha
        self.monotonic_cst = monotonic_cst


class ExtraTreesRegressor(ForestRegressor):
    """
    An extra-trees regressor.

    This class implements a meta estimator that fits a number of
    randomized decision trees (a.k.a. extra-trees) on various sub-samples
    of the dataset and uses averaging to improve the predictive accuracy
    and control over-fitting.

    Read more in the :ref:`User Guide <forest>`.

    Parameters
    ----------
    n_estimators : int, default=100
        The number of trees in the forest.

        .. versionchanged:: 0.22
           The default value of ``n_estimators`` changed from 10 to 100
           in 0.22.

    criterion : {"squared_error", "absolute_error", "friedman_mse", "poisson"}, \
            default="squared_error"
        The function to measure the quality of a split. Supported criteria
        are "squared_error" for the mean squared error, which is equal to
        variance reduction as feature selection criterion and minimizes the L2
        loss using the mean of each terminal node, "friedman_mse", which uses
        mean squared error with Friedman's improvement score for potential
        splits, "absolute_error" for the mean absolute error, which minimizes
        the L1 loss using the median of each terminal node, and "poisson" which
        uses reduction in Poisson deviance to find splits.
        Training using "absolute_error" is significantly slower
        than when using "squared_error".

        .. versionadded:: 0.18
           Mean Absolute Error (MAE) criterion.

    max_depth : int, default=None
        The maximum depth of the tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than
        min_samples_split samples.

    min_samples_split : int or float, default=2
        The minimum number of samples required to split an internal node:

        - If int, then consider `min_samples_split` as the minimum number.
        - If float, then `min_samples_split` is a fraction and
          `ceil(min_samples_split * n_samples)` are the minimum
          number of samples for each split.

        .. versionchanged:: 0.18
           Added float values for fractions.

    min_samples_leaf : int or float, default=1
        The minimum number of samples required to be at a leaf node.
        A split point at any depth will only be considered if it leaves at
        least ``min_samples_leaf`` training samples in each of the left and
        right branches.  This may have the effect of smoothing the model,
        especially in regression.

        - If int, then consider `min_samples_leaf` as the minimum number.
        - If float, then `min_samples_leaf` is a fraction and
          `ceil(min_samples_leaf * n_samples)` are the minimum
          number of samples for each node.

        .. versionchanged:: 0.18
           Added float values for fractions.

    min_weight_fraction_leaf : float, default=0.0
        The minimum weighted fraction of the sum total of weights (of all
        the input samples) required to be at a leaf node. Samples have
        equal weight when sample_weight is not provided.

    max_features : {"sqrt", "log2", None}, int or float, default=1.0
        The number of features to consider when looking for the best split:

        - If int, then consider `max_features` features at each split.
        - If float, then `max_features` is a fraction and
          `max(1, int(max_features * n_features_in_))` features are considered at each
          split.
        - If "sqrt", then `max_features=sqrt(n_features)`.
        - If "log2", then `max_features=log2(n_features)`.
        - If None or 1.0, then `max_features=n_features`.

        .. note::
            The default of 1.0 is equivalent to bagged trees and more
            randomness can be achieved by setting smaller values, e.g. 0.3.

        .. versionchanged:: 1.1
            The default of `max_features` changed from `"auto"` to 1.0.

        Note: the search for a split does not stop until at least one
        valid partition of the node samples is found, even if it requires to
        effectively inspect more than ``max_features`` features.

    max_leaf_nodes : int, default=None
        Grow trees with ``max_leaf_nodes`` in best-first fashion.
        Best nodes are defined as relative reduction in impurity.
        If None then unlimited number of leaf nodes.

    min_impurity_decrease : float, default=0.0
        A node will be split if this split induces a decrease of the impurity
        greater than or equal to this value.

        The weighted impurity decrease equation is the following::

            N_t / N * (impurity - N_t_R / N_t * right_impurity
                                - N_t_L / N_t * left_impurity)

        where ``N`` is the total number of samples, ``N_t`` is the number of
        samples at the current node, ``N_t_L`` is the number of samples in the
        left child, and ``N_t_R`` is the number of samples in the right child.

        ``N``, ``N_t``, ``N_t_R`` and ``N_t_L`` all refer to the weighted sum,
        if ``sample_weight`` is passed.

        .. versionadded:: 0.19

    bootstrap : bool, default=False
        Whether bootstrap samples are used when building trees. If False, the
        whole dataset is used to build each tree.

    oob_score : bool or callable, default=False
        Whether to use out-of-bag samples to estimate the generalization score.
        By default, :func:`~sklearn.metrics.r2_score` is used.
        Provide a callable with signature `metric(y_true, y_pred)` to use a
        custom metric. Only available if `bootstrap=True`.

    n_jobs : int, default=None
        The number of jobs to run in parallel. :meth:`fit`, :meth:`predict`,
        :meth:`decision_path` and :meth:`apply` are all parallelized over the
        trees. ``None`` means 1 unless in a :obj:`joblib.parallel_backend`
        context. ``-1`` means using all processors. See :term:`Glossary
        <n_jobs>` for more details.

    random_state : int, RandomState instance or None, default=None
        Controls 3 sources of randomness:

        - the bootstrapping of the samples used when building trees
          (if ``bootstrap=True``)
        - the sampling of the features to consider when looking for the best
          split at each node (if ``max_features < n_features``)
        - the draw of the splits for each of the `max_features`

        See :term:`Glossary <random_state>` for details.

    verbose : int, default=0
        Controls the verbosity when fitting and predicting.

    warm_start : bool, default=False
        When set to ``True``, reuse the solution of the previous call to fit
        and add more estimators to the ensemble, otherwise, just fit a whole
        new forest. See :term:`Glossary <warm_start>` and
        :ref:`gradient_boosting_warm_start` for details.

    ccp_alpha : non-negative float, default=0.0
        Complexity parameter used for Minimal Cost-Complexity Pruning. The
        subtree with the largest cost complexity that is smaller than
        ``ccp_alpha`` will be chosen. By default, no pruning is performed. See
        :ref:`minimal_cost_complexity_pruning` for details.

        .. versionadded:: 0.22

    max_samples : int or float, default=None
        If bootstrap is True, the number of samples to draw from X
        to train each base estimator.

        - If None (default), then draw `X.shape[0]` samples.
        - If int, then draw `max_samples` samples.
        - If float, then draw `max_samples * X.shape[0]` samples. Thus,
          `max_samples` should be in the interval `(0.0, 1.0]`.

        .. versionadded:: 0.22

    max_bins : int, default=255
        The maximum number of bins to use for non-missing values.

        **This is an experimental feature**.

    store_leaf_values : bool, default=False
        Whether to store the leaf values in the ``get_leaf_node_samples`` function.

        **This is an experimental feature**.

    monotonic_cst : array-like of int of shape (n_features), default=None
        Indicates the monotonicity constraint to enforce on each feature.
          - 1: monotonically increasing
          - 0: no constraint
          - -1: monotonically decreasing

        If monotonic_cst is None, no constraints are applied.

        Monotonicity constraints are not supported for:
          - multioutput regressions (i.e. when `n_outputs_ > 1`),
          - regressions trained on data with missing values.

        Read more in the :ref:`User Guide <monotonic_cst_gbdt>`.

        .. versionadded:: 1.4

    Attributes
    ----------
    estimator_ : :class:`~sklearn.tree.ExtraTreeRegressor`
        The child estimator template used to create the collection of fitted
        sub-estimators.

        .. versionadded:: 1.2
           `base_estimator_` was renamed to `estimator_`.

    base_estimator_ : ExtraTreeRegressor
        The child estimator template used to create the collection of fitted
        sub-estimators.

        .. deprecated:: 1.2
            `base_estimator_` is deprecated and will be removed in 1.4.
            Use `estimator_` instead.

    estimators_ : list of DecisionTreeRegressor
        The collection of fitted sub-estimators.

    feature_importances_ : ndarray of shape (n_features,)
        The impurity-based feature importances.
        The higher, the more important the feature.
        The importance of a feature is computed as the (normalized)
        total reduction of the criterion brought by that feature.  It is also
        known as the Gini importance.

        Warning: impurity-based feature importances can be misleading for
        high cardinality features (many unique values). See
        :func:`sklearn.inspection.permutation_importance` as an alternative.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    n_outputs_ : int
        The number of outputs.

    oob_score_ : float
        Score of the training dataset obtained using an out-of-bag estimate.
        This attribute exists only when ``oob_score`` is True.

    oob_prediction_ : ndarray of shape (n_samples,) or (n_samples, n_outputs)
        Prediction computed with out-of-bag estimate on the training set.
        This attribute exists only when ``oob_score`` is True.

    See Also
    --------
    ExtraTreesClassifier : An extra-trees classifier with random splits.
    RandomForestClassifier : A random forest classifier with optimal splits.
    RandomForestRegressor : Ensemble regressor using trees with optimal splits.

    Notes
    -----
    The default values for the parameters controlling the size of the trees
    (e.g. ``max_depth``, ``min_samples_leaf``, etc.) lead to fully grown and
    unpruned trees which can potentially be very large on some data sets. To
    reduce memory consumption, the complexity and size of the trees should be
    controlled by setting those parameter values.

    References
    ----------
    .. [1] P. Geurts, D. Ernst., and L. Wehenkel, "Extremely randomized trees",
           Machine Learning, 63(1), 3-42, 2006.

    Examples
    --------
    >>> from sklearn.datasets import load_diabetes
    >>> from sklearn.model_selection import train_test_split
    >>> from sklearn.ensemble import ExtraTreesRegressor
    >>> X, y = load_diabetes(return_X_y=True)
    >>> X_train, X_test, y_train, y_test = train_test_split(
    ...     X, y, random_state=0)
    >>> reg = ExtraTreesRegressor(n_estimators=100, random_state=0).fit(
    ...    X_train, y_train)
    >>> reg.score(X_test, y_test)
    0.2727...
    """

    _parameter_constraints: dict = {
        **ForestRegressor._parameter_constraints,
        **DecisionTreeRegressor._parameter_constraints,
    }
    _parameter_constraints.pop("splitter")

    def __init__(
        self,
        n_estimators=100,
        *,
        criterion="squared_error",
        max_depth=None,
        min_samples_split=2,
        min_samples_leaf=1,
        min_weight_fraction_leaf=0.0,
        max_features=1.0,
        max_leaf_nodes=None,
        min_impurity_decrease=0.0,
        bootstrap=False,
        oob_score=False,
        n_jobs=None,
        random_state=None,
        verbose=0,
        warm_start=False,
        ccp_alpha=0.0,
        max_samples=None,
        max_bins=None,
        store_leaf_values=False,
        monotonic_cst=None,
    ):
        super().__init__(
            estimator=ExtraTreeRegressor(),
            n_estimators=n_estimators,
            estimator_params=(
                "criterion",
                "max_depth",
                "min_samples_split",
                "min_samples_leaf",
                "min_weight_fraction_leaf",
                "max_features",
                "max_leaf_nodes",
                "min_impurity_decrease",
                "random_state",
                "ccp_alpha",
                "store_leaf_values",
                "monotonic_cst",
            ),
            bootstrap=bootstrap,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose,
            warm_start=warm_start,
            max_samples=max_samples,
            max_bins=max_bins,
            store_leaf_values=store_leaf_values,
        )

        self.criterion = criterion
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_fraction_leaf = min_weight_fraction_leaf
        self.max_features = max_features
        self.max_leaf_nodes = max_leaf_nodes
        self.min_impurity_decrease = min_impurity_decrease
        self.ccp_alpha = ccp_alpha
        self.monotonic_cst = monotonic_cst


class RandomTreesEmbedding(TransformerMixin, BaseForest):
    """
    An ensemble of totally random trees.

    An unsupervised transformation of a dataset to a high-dimensional
    sparse representation. A datapoint is coded according to which leaf of
    each tree it is sorted into. Using a one-hot encoding of the leaves,
    this leads to a binary coding with as many ones as there are trees in
    the forest.

    The dimensionality of the resulting representation is
    ``n_out <= n_estimators * max_leaf_nodes``. If ``max_leaf_nodes == None``,
    the number of leaf nodes is at most ``n_estimators * 2 ** max_depth``.

    Read more in the :ref:`User Guide <random_trees_embedding>`.

    Parameters
    ----------
    n_estimators : int, default=100
        Number of trees in the forest.

        .. versionchanged:: 0.22
           The default value of ``n_estimators`` changed from 10 to 100
           in 0.22.

    max_depth : int, default=5
        The maximum depth of each tree. If None, then nodes are expanded until
        all leaves are pure or until all leaves contain less than
        min_samples_split samples.

    min_samples_split : int or float, default=2
        The minimum number of samples required to split an internal node:

        - If int, then consider `min_samples_split` as the minimum number.
        - If float, then `min_samples_split` is a fraction and
          `ceil(min_samples_split * n_samples)` is the minimum
          number of samples for each split.

        .. versionchanged:: 0.18
           Added float values for fractions.

    min_samples_leaf : int or float, default=1
        The minimum number of samples required to be at a leaf node.
        A split point at any depth will only be considered if it leaves at
        least ``min_samples_leaf`` training samples in each of the left and
        right branches.  This may have the effect of smoothing the model,
        especially in regression.

        - If int, then consider `min_samples_leaf` as the minimum number.
        - If float, then `min_samples_leaf` is a fraction and
          `ceil(min_samples_leaf * n_samples)` is the minimum
          number of samples for each node.

        .. versionchanged:: 0.18
           Added float values for fractions.

    min_weight_fraction_leaf : float, default=0.0
        The minimum weighted fraction of the sum total of weights (of all
        the input samples) required to be at a leaf node. Samples have
        equal weight when sample_weight is not provided.

    max_leaf_nodes : int, default=None
        Grow trees with ``max_leaf_nodes`` in best-first fashion.
        Best nodes are defined as relative reduction in impurity.
        If None then unlimited number of leaf nodes.

    min_impurity_decrease : float, default=0.0
        A node will be split if this split induces a decrease of the impurity
        greater than or equal to this value.

        The weighted impurity decrease equation is the following::

            N_t / N * (impurity - N_t_R / N_t * right_impurity
                                - N_t_L / N_t * left_impurity)

        where ``N`` is the total number of samples, ``N_t`` is the number of
        samples at the current node, ``N_t_L`` is the number of samples in the
        left child, and ``N_t_R`` is the number of samples in the right child.

        ``N``, ``N_t``, ``N_t_R`` and ``N_t_L`` all refer to the weighted sum,
        if ``sample_weight`` is passed.

        .. versionadded:: 0.19

    sparse_output : bool, default=True
        Whether or not to return a sparse CSR matrix, as default behavior,
        or to return a dense array compatible with dense pipeline operators.

    n_jobs : int, default=None
        The number of jobs to run in parallel. :meth:`fit`, :meth:`transform`,
        :meth:`decision_path` and :meth:`apply` are all parallelized over the
        trees. ``None`` means 1 unless in a :obj:`joblib.parallel_backend`
        context. ``-1`` means using all processors. See :term:`Glossary
        <n_jobs>` for more details.

    random_state : int, RandomState instance or None, default=None
        Controls the generation of the random `y` used to fit the trees
        and the draw of the splits for each feature at the trees' nodes.
        See :term:`Glossary <random_state>` for details.

    verbose : int, default=0
        Controls the verbosity when fitting and predicting.

    warm_start : bool, default=False
        When set to ``True``, reuse the solution of the previous call to fit
        and add more estimators to the ensemble, otherwise, just fit a whole
        new forest. See :term:`Glossary <warm_start>` and
        :ref:`gradient_boosting_warm_start` for details.

    store_leaf_values : bool, default=False
        Whether to store the leaf values in the ``get_leaf_node_samples`` function.

    Attributes
    ----------
    estimator_ : :class:`~sklearn.tree.ExtraTreeRegressor` instance
        The child estimator template used to create the collection of fitted
        sub-estimators.

        .. versionadded:: 1.2
           `base_estimator_` was renamed to `estimator_`.

    base_estimator_ : :class:`~sklearn.tree.ExtraTreeRegressor` instance
        The child estimator template used to create the collection of fitted
        sub-estimators.

        .. deprecated:: 1.2
            `base_estimator_` is deprecated and will be removed in 1.4.
            Use `estimator_` instead.

    estimators_ : list of :class:`~sklearn.tree.ExtraTreeRegressor` instances
        The collection of fitted sub-estimators.

    feature_importances_ : ndarray of shape (n_features,)
        The feature importances (the higher, the more important the feature).

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    n_outputs_ : int
        The number of outputs when ``fit`` is performed.

    one_hot_encoder_ : OneHotEncoder instance
        One-hot encoder used to create the sparse embedding.

    See Also
    --------
    ExtraTreesClassifier : An extra-trees classifier.
    ExtraTreesRegressor : An extra-trees regressor.
    RandomForestClassifier : A random forest classifier.
    RandomForestRegressor : A random forest regressor.
    sklearn.tree.ExtraTreeClassifier: An extremely randomized
        tree classifier.
    sklearn.tree.ExtraTreeRegressor : An extremely randomized
        tree regressor.

    References
    ----------
    .. [1] P. Geurts, D. Ernst., and L. Wehenkel, "Extremely randomized trees",
           Machine Learning, 63(1), 3-42, 2006.
    .. [2] Moosmann, F. and Triggs, B. and Jurie, F.  "Fast discriminative
           visual codebooks using randomized clustering forests"
           NIPS 2007

    Examples
    --------
    >>> from sklearn.ensemble import RandomTreesEmbedding
    >>> X = [[0,0], [1,0], [0,1], [-1,0], [0,-1]]
    >>> random_trees = RandomTreesEmbedding(
    ...    n_estimators=5, random_state=0, max_depth=1).fit(X)
    >>> X_sparse_embedding = random_trees.transform(X)
    >>> X_sparse_embedding.toarray()
    array([[0., 1., 1., 0., 1., 0., 0., 1., 1., 0.],
           [0., 1., 1., 0., 1., 0., 0., 1., 1., 0.],
           [0., 1., 0., 1., 0., 1., 0., 1., 0., 1.],
           [1., 0., 1., 0., 1., 0., 1., 0., 1., 0.],
           [0., 1., 1., 0., 1., 0., 0., 1., 1., 0.]])
    """

    _parameter_constraints: dict = {
        "n_estimators": [Interval(Integral, 1, None, closed="left")],
        "n_jobs": [Integral, None],
        "verbose": ["verbose"],
        "warm_start": ["boolean"],
        **BaseDecisionTree._parameter_constraints,
        "sparse_output": ["boolean"],
    }
    for param in ("max_features", "ccp_alpha", "splitter", "monotonic_cst"):
        _parameter_constraints.pop(param)

    criterion = "squared_error"
    max_features = 1

    def __init__(
        self,
        n_estimators=100,
        *,
        max_depth=5,
        min_samples_split=2,
        min_samples_leaf=1,
        min_weight_fraction_leaf=0.0,
        max_leaf_nodes=None,
        min_impurity_decrease=0.0,
        sparse_output=True,
        n_jobs=None,
        random_state=None,
        verbose=0,
        warm_start=False,
        store_leaf_values=False,
    ):
        super().__init__(
            estimator=ExtraTreeRegressor(),
            n_estimators=n_estimators,
            estimator_params=(
                "criterion",
                "max_depth",
                "min_samples_split",
                "min_samples_leaf",
                "min_weight_fraction_leaf",
                "max_features",
                "max_leaf_nodes",
                "min_impurity_decrease",
                "random_state",
                "store_leaf_values",
            ),
            bootstrap=False,
            oob_score=False,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose,
            warm_start=warm_start,
            max_samples=None,
            store_leaf_values=store_leaf_values,
        )

        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_fraction_leaf = min_weight_fraction_leaf
        self.max_leaf_nodes = max_leaf_nodes
        self.min_impurity_decrease = min_impurity_decrease
        self.sparse_output = sparse_output

    def _set_oob_score_and_attributes(self, X, y, scoring_function=None):
        raise NotImplementedError("OOB score not supported by tree embedding")

    def fit(self, X, y=None, sample_weight=None, classes=None):
        """
        Fit estimator.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The input samples. Use ``dtype=np.float32`` for maximum
            efficiency. Sparse matrices are also supported, use sparse
            ``csc_matrix`` for maximum efficiency.

        y : Ignored
            Not used, present for API consistency by convention.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights. If None, then samples are equally weighted. Splits
            that would create child nodes with net zero or negative weight are
            ignored while searching for a split in each node. In the case of
            classification, splits are also ignored if they would result in any
            single class carrying a negative weight in either child node.

        classes : array-like of shape (n_classes,), default=None
            List of all the classes that can possibly appear in the y vector.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        # Parameters are validated in fit_transform
        self.fit_transform(X, y, sample_weight=sample_weight, classes=classes)
        return self

    @_fit_context(prefer_skip_nested_validation=True)
    def fit_transform(self, X, y=None, sample_weight=None, classes=None):
        """
        Fit estimator and transform dataset.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Input data used to build forests. Use ``dtype=np.float32`` for
            maximum efficiency.

        y : Ignored
            Not used, present for API consistency by convention.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights. If None, then samples are equally weighted. Splits
            that would create child nodes with net zero or negative weight are
            ignored while searching for a split in each node. In the case of
            classification, splits are also ignored if they would result in any
            single class carrying a negative weight in either child node.

        classes : array-like of shape (n_classes,), default=None
            List of all the classes that can possibly appear in the y vector.

        Returns
        -------
        X_transformed : sparse matrix of shape (n_samples, n_out)
            Transformed dataset.
        """
        rnd = check_random_state(self.random_state)
        y = rnd.uniform(size=_num_samples(X))
        super().fit(X, y, sample_weight=sample_weight, classes=classes)

        self.one_hot_encoder_ = OneHotEncoder(sparse_output=self.sparse_output)
        output = self.one_hot_encoder_.fit_transform(self.apply(X))
        self._n_features_out = output.shape[1]
        return output

    def get_feature_names_out(self, input_features=None):
        """Get output feature names for transformation.

        Parameters
        ----------
        input_features : array-like of str or None, default=None
            Only used to validate feature names with the names seen in :meth:`fit`.

        Returns
        -------
        feature_names_out : ndarray of str objects
            Transformed feature names, in the format of
            `randomtreesembedding_{tree}_{leaf}`, where `tree` is the tree used
            to generate the leaf and `leaf` is the index of a leaf node
            in that tree. Note that the node indexing scheme is used to
            index both nodes with children (split nodes) and leaf nodes.
            Only the latter can be present as output features.
            As a consequence, there are missing indices in the output
            feature names.
        """
        check_is_fitted(self, "_n_features_out")
        _check_feature_names_in(
            self, input_features=input_features, generate_names=False
        )

        feature_names = [
            f"randomtreesembedding_{tree}_{leaf}"
            for tree in range(self.n_estimators)
            for leaf in self.one_hot_encoder_.categories_[tree]
        ]
        return np.asarray(feature_names, dtype=object)

    def transform(self, X):
        """
        Transform dataset.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Input data to be transformed. Use ``dtype=np.float32`` for maximum
            efficiency. Sparse matrices are also supported, use sparse
            ``csr_matrix`` for maximum efficiency.

        Returns
        -------
        X_transformed : sparse matrix of shape (n_samples, n_out)
            Transformed dataset.
        """
        check_is_fitted(self)
        return self.one_hot_encoder_.transform(self.apply(X))
