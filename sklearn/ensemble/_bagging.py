"""Bagging meta-estimator."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import itertools
import numbers
from abc import ABCMeta, abstractmethod
from functools import partial
from numbers import Integral
from warnings import warn

import numpy as np

from sklearn.base import ClassifierMixin, RegressorMixin, _fit_context
from sklearn.ensemble._base import BaseEnsemble, _partition_estimators
from sklearn.metrics import accuracy_score, r2_score
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.utils import Bunch, _safe_indexing, check_random_state, column_or_1d
from sklearn.utils._mask import indices_to_mask
from sklearn.utils._param_validation import HasMethods, Interval, RealNotInt
from sklearn.utils._tags import get_tags
from sklearn.utils.metadata_routing import (
    MetadataRouter,
    MethodMapping,
    _raise_for_params,
    _routing_enabled,
    get_routing_for_object,
    process_routing,
)
from sklearn.utils.metaestimators import available_if
from sklearn.utils.multiclass import check_classification_targets
from sklearn.utils.parallel import Parallel, delayed
from sklearn.utils.random import sample_without_replacement
from sklearn.utils.validation import (
    _check_method_params,
    _check_sample_weight,
    _estimator_has,
    check_is_fitted,
    has_fit_parameter,
    validate_data,
)

__all__ = ["BaggingClassifier", "BaggingRegressor"]

MAX_INT = np.iinfo(np.int32).max


def _generate_indices(random_state, bootstrap, n_population, n_samples):
    """Draw randomly sampled indices."""
    # Draw sample indices
    if bootstrap:
        indices = random_state.randint(0, n_population, n_samples)
    else:
        indices = sample_without_replacement(
            n_population, n_samples, random_state=random_state
        )

    return indices


def _generate_bagging_indices(
    random_state,
    bootstrap_features,
    bootstrap_samples,
    n_features,
    n_samples,
    max_features,
    max_samples,
    sample_weight,
):
    """Randomly draw feature and sample indices."""
    # Get valid random state
    random_state = check_random_state(random_state)

    # Draw indices
    feature_indices = _generate_indices(
        random_state, bootstrap_features, n_features, max_features
    )
    if sample_weight is None:
        sample_indices = _generate_indices(
            random_state, bootstrap_samples, n_samples, max_samples
        )
    else:
        normalized_sample_weight = sample_weight / np.sum(sample_weight)
        sample_indices = random_state.choice(
            n_samples,
            max_samples,
            replace=bootstrap_samples,
            p=normalized_sample_weight,
        )

    return feature_indices, sample_indices


def _consumes_sample_weight(estimator):
    if _routing_enabled():
        request_or_router = get_routing_for_object(estimator)
        consumes_sample_weight = request_or_router.consumes("fit", ("sample_weight",))
    else:
        consumes_sample_weight = has_fit_parameter(estimator, "sample_weight")
    return consumes_sample_weight


def _parallel_build_estimators(
    n_estimators,
    ensemble,
    X,
    y,
    sample_weight,
    seeds,
    total_n_estimators,
    verbose,
    check_input,
    fit_params,
):
    """Private function used to build a batch of estimators within a job."""
    # Retrieve settings
    n_samples, n_features = X.shape
    max_features = ensemble._max_features
    max_samples = ensemble._max_samples
    bootstrap = ensemble.bootstrap
    bootstrap_features = ensemble.bootstrap_features
    has_check_input = has_fit_parameter(ensemble.estimator_, "check_input")
    requires_feature_indexing = bootstrap_features or max_features != n_features
    consumes_sample_weight = _consumes_sample_weight(ensemble.estimator_)

    # Build estimators
    estimators = []
    estimators_features = []

    for i in range(n_estimators):
        if verbose > 1:
            print(
                "Building estimator %d of %d for this parallel run (total %d)..."
                % (i + 1, n_estimators, total_n_estimators)
            )

        random_state = seeds[i]
        estimator = ensemble._make_estimator(append=False, random_state=random_state)

        if has_check_input:
            estimator_fit = partial(estimator.fit, check_input=check_input)
        else:
            estimator_fit = estimator.fit

        # Draw random feature, sample indices (using normalized sample_weight
        # as probabilities if provided).
        features, indices = _generate_bagging_indices(
            random_state,
            bootstrap_features,
            bootstrap,
            n_features,
            n_samples,
            max_features,
            max_samples,
            sample_weight,
        )

        fit_params_ = fit_params.copy()

        # Note: Row sampling can be achieved either through setting sample_weight or
        # by indexing. The former is more memory efficient. Therefore, use this method
        # if possible, otherwise use indexing.
        if consumes_sample_weight:
            # Row sampling by setting sample_weight
            indices_as_sample_weight = np.bincount(indices, minlength=n_samples)
            fit_params_["sample_weight"] = indices_as_sample_weight
            X_ = X[:, features] if requires_feature_indexing else X
            estimator_fit(X_, y, **fit_params_)
        else:
            # Row sampling by indexing
            y_ = _safe_indexing(y, indices)
            X_ = _safe_indexing(X, indices)
            fit_params_ = _check_method_params(X, params=fit_params_, indices=indices)
            if requires_feature_indexing:
                X_ = X_[:, features]
            estimator_fit(X_, y_, **fit_params_)

        estimators.append(estimator)
        estimators_features.append(features)

    return estimators, estimators_features


def _parallel_predict_proba(
    estimators,
    estimators_features,
    X,
    n_classes,
    predict_params=None,
    predict_proba_params=None,
):
    """Private function used to compute (proba-)predictions within a job."""
    n_samples = X.shape[0]
    proba = np.zeros((n_samples, n_classes))

    for estimator, features in zip(estimators, estimators_features):
        if hasattr(estimator, "predict_proba"):
            proba_estimator = estimator.predict_proba(
                X[:, features], **(predict_params or {})
            )

            if n_classes == len(estimator.classes_):
                proba += proba_estimator

            else:
                proba[:, estimator.classes_] += proba_estimator[
                    :, range(len(estimator.classes_))
                ]

        else:
            # Resort to voting
            predictions = estimator.predict(
                X[:, features], **(predict_proba_params or {})
            )

            for i in range(n_samples):
                proba[i, predictions[i]] += 1

    return proba


def _parallel_predict_log_proba(estimators, estimators_features, X, n_classes, params):
    """Private function used to compute log probabilities within a job."""
    n_samples = X.shape[0]
    log_proba = np.empty((n_samples, n_classes))
    log_proba.fill(-np.inf)
    all_classes = np.arange(n_classes, dtype=int)

    for estimator, features in zip(estimators, estimators_features):
        log_proba_estimator = estimator.predict_log_proba(X[:, features], **params)

        if n_classes == len(estimator.classes_):
            log_proba = np.logaddexp(log_proba, log_proba_estimator)

        else:
            log_proba[:, estimator.classes_] = np.logaddexp(
                log_proba[:, estimator.classes_],
                log_proba_estimator[:, range(len(estimator.classes_))],
            )

            missing = np.setdiff1d(all_classes, estimator.classes_)
            log_proba[:, missing] = np.logaddexp(log_proba[:, missing], -np.inf)

    return log_proba


def _parallel_decision_function(estimators, estimators_features, X, params):
    """Private function used to compute decisions within a job."""
    return sum(
        estimator.decision_function(X[:, features], **params)
        for estimator, features in zip(estimators, estimators_features)
    )


def _parallel_predict_regression(estimators, estimators_features, X, params):
    """Private function used to compute predictions within a job."""
    return sum(
        estimator.predict(X[:, features], **params)
        for estimator, features in zip(estimators, estimators_features)
    )


class BaseBagging(BaseEnsemble, metaclass=ABCMeta):
    """Base class for Bagging meta-estimator.

    Warning: This class should not be used directly. Use derived classes
    instead.
    """

    _parameter_constraints: dict = {
        "estimator": [HasMethods(["fit", "predict"]), None],
        "n_estimators": [Interval(Integral, 1, None, closed="left")],
        "max_samples": [
            Interval(Integral, 1, None, closed="left"),
            Interval(RealNotInt, 0, 1, closed="right"),
        ],
        "max_features": [
            Interval(Integral, 1, None, closed="left"),
            Interval(RealNotInt, 0, 1, closed="right"),
        ],
        "bootstrap": ["boolean"],
        "bootstrap_features": ["boolean"],
        "oob_score": ["boolean"],
        "warm_start": ["boolean"],
        "n_jobs": [None, Integral],
        "random_state": ["random_state"],
        "verbose": ["verbose"],
    }

    @abstractmethod
    def __init__(
        self,
        estimator=None,
        n_estimators=10,
        *,
        max_samples=1.0,
        max_features=1.0,
        bootstrap=True,
        bootstrap_features=False,
        oob_score=False,
        warm_start=False,
        n_jobs=None,
        random_state=None,
        verbose=0,
    ):
        super().__init__(
            estimator=estimator,
            n_estimators=n_estimators,
        )
        self.max_samples = max_samples
        self.max_features = max_features
        self.bootstrap = bootstrap
        self.bootstrap_features = bootstrap_features
        self.oob_score = oob_score
        self.warm_start = warm_start
        self.n_jobs = n_jobs
        self.random_state = random_state
        self.verbose = verbose

    @_fit_context(
        # BaseBagging.estimator is not validated yet
        prefer_skip_nested_validation=False
    )
    def fit(self, X, y, sample_weight=None, **fit_params):
        """Build a Bagging ensemble of estimators from the training set (X, y).

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The training input samples. Sparse matrices are accepted only if
            they are supported by the base estimator.

        y : array-like of shape (n_samples,)
            The target values (class labels in classification, real numbers in
            regression).

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights. If None, then samples are equally weighted. Used as
            probabilities to sample the training set. Note that the expected
            frequency semantics for the `sample_weight` parameter are only
            fulfilled when sampling with replacement `bootstrap=True`.

        **fit_params : dict
            Parameters to pass to the underlying estimators.

            .. versionadded:: 1.5

                Only available if `enable_metadata_routing=True`,
                which can be set by using
                ``sklearn.set_config(enable_metadata_routing=True)``.
                See :ref:`Metadata Routing User Guide <metadata_routing>` for
                more details.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        _raise_for_params(fit_params, self, "fit")

        # Convert data (X is required to be 2d and indexable)
        X, y = validate_data(
            self,
            X,
            y,
            accept_sparse=["csr", "csc"],
            dtype=None,
            ensure_all_finite=False,
            multi_output=True,
        )

        if sample_weight is not None:
            sample_weight = _check_sample_weight(sample_weight, X, dtype=None)

            if not self.bootstrap:
                warn(
                    f"When fitting {self.__class__.__name__} with sample_weight "
                    f"it is recommended to use bootstrap=True, got {self.bootstrap}."
                )

        return self._fit(
            X,
            y,
            max_samples=self.max_samples,
            sample_weight=sample_weight,
            **fit_params,
        )

    def _parallel_args(self):
        return {}

    def _fit(
        self,
        X,
        y,
        max_samples=None,
        max_depth=None,
        check_input=True,
        sample_weight=None,
        **fit_params,
    ):
        """Build a Bagging ensemble of estimators from the training
           set (X, y).

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The training input samples. Sparse matrices are accepted only if
            they are supported by the base estimator.

        y : array-like of shape (n_samples,)
            The target values (class labels in classification, real numbers in
            regression).

        max_samples : int or float, default=None
            Argument to use instead of self.max_samples.

        max_depth : int, default=None
            Override value used when constructing base estimator. Only
            supported if the base estimator has a max_depth parameter.

        check_input : bool, default=True
            Override value used when fitting base estimator. Only supported
            if the base estimator has a check_input parameter for fit function.
            If the meta-estimator already checks the input, set this value to
            False to prevent redundant input validation.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights. If None, then samples are equally weighted.

        **fit_params : dict, default=None
            Parameters to pass to the :term:`fit` method of the underlying
            estimator.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        random_state = check_random_state(self.random_state)

        # Remap output
        n_samples = X.shape[0]
        self._n_samples = n_samples
        y = self._validate_y(y)

        # Check parameters
        self._validate_estimator(self._get_estimator())

        if _routing_enabled():
            routed_params = process_routing(self, "fit", **fit_params)
        else:
            routed_params = Bunch()
            routed_params.estimator = Bunch(fit=fit_params)

        if max_depth is not None:
            self.estimator_.max_depth = max_depth

        # Validate max_samples
        if max_samples is None:
            max_samples = self.max_samples

        if not isinstance(max_samples, numbers.Integral):
            if sample_weight is None:
                max_samples = max(int(max_samples * X.shape[0]), 1)
            else:
                sw_sum = np.sum(sample_weight)
                if sw_sum <= 1:
                    raise ValueError(
                        f"The total sum of sample weights is {sw_sum}, which prevents "
                        "resampling with a fractional value for max_samples="
                        f"{max_samples}. Either pass max_samples as an integer or "
                        "use a larger sample_weight."
                    )
                max_samples = max(int(max_samples * sw_sum), 1)

        if not self.bootstrap and max_samples > X.shape[0]:
            raise ValueError(
                f"Effective max_samples={max_samples} must be <= n_samples="
                f"{X.shape[0]} to be able to sample without replacement."
            )

        # Store validated integer row sampling value
        self._max_samples = max_samples

        # Validate max_features
        if isinstance(self.max_features, numbers.Integral):
            max_features = self.max_features
        elif isinstance(self.max_features, float):
            max_features = int(self.max_features * self.n_features_in_)

        if max_features > self.n_features_in_:
            raise ValueError("max_features must be <= n_features")

        max_features = max(1, int(max_features))

        # Store validated integer feature sampling value
        self._max_features = max_features

        # Store sample_weight (needed in _get_estimators_indices). Note that
        # we intentionally do not materialize `sample_weight=None` as an array
        # of ones to avoid unnecessarily cluttering trained estimator pickles.
        self._sample_weight = sample_weight

        # Other checks
        if not self.bootstrap and self.oob_score:
            raise ValueError("Out of bag estimation only available if bootstrap=True")

        if self.warm_start and self.oob_score:
            raise ValueError("Out of bag estimate only available if warm_start=False")

        if hasattr(self, "oob_score_") and self.warm_start:
            del self.oob_score_

        if not self.warm_start or not hasattr(self, "estimators_"):
            # Free allocated memory, if any
            self.estimators_ = []
            self.estimators_features_ = []

        n_more_estimators = self.n_estimators - len(self.estimators_)

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
            return self

        # Parallel loop
        n_jobs, n_estimators, starts = _partition_estimators(
            n_more_estimators, self.n_jobs
        )
        total_n_estimators = sum(n_estimators)

        # Advance random state to state after training
        # the first n_estimators
        if self.warm_start and len(self.estimators_) > 0:
            random_state.randint(MAX_INT, size=len(self.estimators_))

        seeds = random_state.randint(MAX_INT, size=n_more_estimators)
        self._seeds = seeds

        all_results = Parallel(
            n_jobs=n_jobs, verbose=self.verbose, **self._parallel_args()
        )(
            delayed(_parallel_build_estimators)(
                n_estimators[i],
                self,
                X,
                y,
                sample_weight,
                seeds[starts[i] : starts[i + 1]],
                total_n_estimators,
                verbose=self.verbose,
                check_input=check_input,
                fit_params=routed_params.estimator.fit,
            )
            for i in range(n_jobs)
        )

        # Reduce
        self.estimators_ += list(
            itertools.chain.from_iterable(t[0] for t in all_results)
        )
        self.estimators_features_ += list(
            itertools.chain.from_iterable(t[1] for t in all_results)
        )

        if self.oob_score:
            self._set_oob_score(X, y)

        return self

    @abstractmethod
    def _set_oob_score(self, X, y):
        """Calculate out of bag predictions and score."""

    def _validate_y(self, y):
        if len(y.shape) == 1 or y.shape[1] == 1:
            return column_or_1d(y, warn=True)
        return y

    def _get_estimators_indices(self):
        # Get drawn indices along both sample and feature axes
        for seed in self._seeds:
            # Operations accessing random_state must be performed identically
            # to those in `_parallel_build_estimators()`
            feature_indices, sample_indices = _generate_bagging_indices(
                seed,
                self.bootstrap_features,
                self.bootstrap,
                self.n_features_in_,
                self._n_samples,
                self._max_features,
                self._max_samples,
                self._sample_weight,
            )

            yield feature_indices, sample_indices

    @property
    def estimators_samples_(self):
        """
        The subset of drawn samples for each base estimator.

        Returns a dynamically generated list of indices identifying
        the samples used for fitting each member of the ensemble, i.e.,
        the in-bag samples.

        Note: the list is re-created at each call to the property in order
        to reduce the object memory footprint by not storing the sampling
        data. Thus fetching the property may be slower than expected.
        """
        return [sample_indices for _, sample_indices in self._get_estimators_indices()]

    def get_metadata_routing(self):
        """Get metadata routing of this object.

        Please check :ref:`User Guide <metadata_routing>` on how the routing
        mechanism works.

        .. versionadded:: 1.5

        Returns
        -------
        routing : MetadataRouter
            A :class:`~sklearn.utils.metadata_routing.MetadataRouter` encapsulating
            routing information.
        """
        router = MetadataRouter(owner=self.__class__.__name__)

        method_mapping = MethodMapping()
        method_mapping.add(caller="fit", callee="fit").add(
            caller="decision_function", callee="decision_function"
        )

        # the router needs to be built depending on whether the sub-estimator has a
        # `predict_proba` method (as BaggingClassifier decides dynamically at runtime):
        if hasattr(self._get_estimator(), "predict_proba"):
            (
                method_mapping.add(caller="predict", callee="predict_proba").add(
                    caller="predict_proba", callee="predict_proba"
                )
            )

        else:
            (
                method_mapping.add(caller="predict", callee="predict").add(
                    caller="predict_proba", callee="predict"
                )
            )

        # the router needs to be built depending on whether the sub-estimator has a
        # `predict_log_proba` method (as BaggingClassifier decides dynamically at
        # runtime):
        if hasattr(self._get_estimator(), "predict_log_proba"):
            method_mapping.add(caller="predict_log_proba", callee="predict_log_proba")

        else:
            # if `predict_log_proba` is not available in BaggingClassifier's
            # sub-estimator, the routing should go to its `predict_proba` if it is
            # available or else to its `predict` method; according to how
            # `sample_weight` is passed to the respective methods dynamically at
            # runtime:
            if hasattr(self._get_estimator(), "predict_proba"):
                method_mapping.add(caller="predict_log_proba", callee="predict_proba")

            else:
                method_mapping.add(caller="predict_log_proba", callee="predict")

        router.add(estimator=self._get_estimator(), method_mapping=method_mapping)
        return router

    @abstractmethod
    def _get_estimator(self):
        """Resolve which estimator to return."""

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.input_tags.sparse = get_tags(self._get_estimator()).input_tags.sparse
        tags.input_tags.allow_nan = get_tags(self._get_estimator()).input_tags.allow_nan
        return tags


class BaggingClassifier(ClassifierMixin, BaseBagging):
    """A Bagging classifier.

    A Bagging classifier is an ensemble meta-estimator that fits base
    classifiers each on random subsets of the original dataset and then
    aggregate their individual predictions (either by voting or by averaging)
    to form a final prediction. Such a meta-estimator can typically be used as
    a way to reduce the variance of a black-box estimator (e.g., a decision
    tree), by introducing randomization into its construction procedure and
    then making an ensemble out of it.

    This algorithm encompasses several works from the literature. When random
    subsets of the dataset are drawn as random subsets of the samples, then
    this algorithm is known as Pasting [1]_. If samples are drawn with
    replacement, then the method is known as Bagging [2]_. When random subsets
    of the dataset are drawn as random subsets of the features, then the method
    is known as Random Subspaces [3]_. Finally, when base estimators are built
    on subsets of both samples and features, then the method is known as
    Random Patches [4]_.

    Read more in the :ref:`User Guide <bagging>`.

    .. versionadded:: 0.15

    Parameters
    ----------
    estimator : object, default=None
        The base estimator to fit on random subsets of the dataset.
        If None, then the base estimator is a
        :class:`~sklearn.tree.DecisionTreeClassifier`.

        .. versionadded:: 1.2
           `base_estimator` was renamed to `estimator`.

    n_estimators : int, default=10
        The number of base estimators in the ensemble.

    max_samples : int or float, default=1.0
        The number of samples to draw from X to train each base estimator (with
        replacement by default, see `bootstrap` for more details).

        - If int, then draw `max_samples` samples.
        - If float, then draw `max_samples * X.shape[0]` unweighted samples
          or `max_samples * sample_weight.sum()` weighted samples.

    max_features : int or float, default=1.0
        The number of features to draw from X to train each base estimator (
        without replacement by default, see `bootstrap_features` for more
        details).

        - If int, then draw `max_features` features.
        - If float, then draw `max(1, int(max_features * n_features_in_))` features.

    bootstrap : bool, default=True
        Whether samples are drawn with replacement. If False, sampling without
        replacement is performed. If fitting with `sample_weight`, it is
        strongly recommended to choose True, as only drawing with replacement
        will ensure the expected frequency semantics of `sample_weight`.

    bootstrap_features : bool, default=False
        Whether features are drawn with replacement.

    oob_score : bool, default=False
        Whether to use out-of-bag samples to estimate
        the generalization error. Only available if bootstrap=True.

    warm_start : bool, default=False
        When set to True, reuse the solution of the previous call to fit
        and add more estimators to the ensemble, otherwise, just fit
        a whole new ensemble. See :term:`the Glossary <warm_start>`.

        .. versionadded:: 0.17
           *warm_start* constructor parameter.

    n_jobs : int, default=None
        The number of jobs to run in parallel for both :meth:`fit` and
        :meth:`predict`. ``None`` means 1 unless in a
        :obj:`joblib.parallel_backend` context. ``-1`` means using all
        processors. See :term:`Glossary <n_jobs>` for more details.

    random_state : int, RandomState instance or None, default=None
        Controls the random resampling of the original dataset
        (sample wise and feature wise).
        If the base estimator accepts a `random_state` attribute, a different
        seed is generated for each instance in the ensemble.
        Pass an int for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    verbose : int, default=0
        Controls the verbosity when fitting and predicting.

    Attributes
    ----------
    estimator_ : estimator
        The base estimator from which the ensemble is grown.

        .. versionadded:: 1.2
           `base_estimator_` was renamed to `estimator_`.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    estimators_ : list of estimators
        The collection of fitted base estimators.

    estimators_samples_ : list of arrays
        The subset of drawn samples (i.e., the in-bag samples) for each base
        estimator. Each subset is defined by an array of the indices selected.

    estimators_features_ : list of arrays
        The subset of drawn features for each base estimator.

    classes_ : ndarray of shape (n_classes,)
        The classes labels.

    n_classes_ : int or list
        The number of classes.

    oob_score_ : float
        Score of the training dataset obtained using an out-of-bag estimate.
        This attribute exists only when ``oob_score`` is True.

    oob_decision_function_ : ndarray of shape (n_samples, n_classes)
        Decision function computed with out-of-bag estimate on the training
        set. If n_estimators is small it might be possible that a data point
        was never left out during the bootstrap. In this case,
        `oob_decision_function_` might contain NaN. This attribute exists
        only when ``oob_score`` is True.

    See Also
    --------
    BaggingRegressor : A Bagging regressor.

    References
    ----------

    .. [1] L. Breiman, "Pasting small votes for classification in large
           databases and on-line", Machine Learning, 36(1), 85-103, 1999.

    .. [2] L. Breiman, "Bagging predictors", Machine Learning, 24(2), 123-140,
           1996.

    .. [3] T. Ho, "The random subspace method for constructing decision
           forests", Pattern Analysis and Machine Intelligence, 20(8), 832-844,
           1998.

    .. [4] G. Louppe and P. Geurts, "Ensembles on Random Patches", Machine
           Learning and Knowledge Discovery in Databases, 346-361, 2012.

    Examples
    --------
    >>> from sklearn.svm import SVC
    >>> from sklearn.ensemble import BaggingClassifier
    >>> from sklearn.datasets import make_classification
    >>> X, y = make_classification(n_samples=100, n_features=4,
    ...                            n_informative=2, n_redundant=0,
    ...                            random_state=0, shuffle=False)
    >>> clf = BaggingClassifier(estimator=SVC(),
    ...                         n_estimators=10, random_state=0).fit(X, y)
    >>> clf.predict([[0, 0, 0, 0]])
    array([1])
    """

    def __init__(
        self,
        estimator=None,
        n_estimators=10,
        *,
        max_samples=1.0,
        max_features=1.0,
        bootstrap=True,
        bootstrap_features=False,
        oob_score=False,
        warm_start=False,
        n_jobs=None,
        random_state=None,
        verbose=0,
    ):
        super().__init__(
            estimator=estimator,
            n_estimators=n_estimators,
            max_samples=max_samples,
            max_features=max_features,
            bootstrap=bootstrap,
            bootstrap_features=bootstrap_features,
            oob_score=oob_score,
            warm_start=warm_start,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose,
        )

    def _get_estimator(self):
        """Resolve which estimator to return (default is DecisionTreeClassifier)"""
        if self.estimator is None:
            return DecisionTreeClassifier()
        return self.estimator

    def _set_oob_score(self, X, y):
        n_samples = y.shape[0]
        n_classes_ = self.n_classes_

        predictions = np.zeros((n_samples, n_classes_))

        for estimator, samples, features in zip(
            self.estimators_, self.estimators_samples_, self.estimators_features_
        ):
            # Create mask for OOB samples
            mask = ~indices_to_mask(samples, n_samples)

            if hasattr(estimator, "predict_proba"):
                predictions[mask, :] += estimator.predict_proba(
                    (X[mask, :])[:, features]
                )

            else:
                p = estimator.predict((X[mask, :])[:, features])
                j = 0

                for i in range(n_samples):
                    if mask[i]:
                        predictions[i, p[j]] += 1
                        j += 1

        if (predictions.sum(axis=1) == 0).any():
            warn(
                "Some inputs do not have OOB scores. "
                "This probably means too few estimators were used "
                "to compute any reliable oob estimates."
            )

        oob_decision_function = predictions / predictions.sum(axis=1)[:, np.newaxis]
        oob_score = accuracy_score(y, np.argmax(predictions, axis=1))

        self.oob_decision_function_ = oob_decision_function
        self.oob_score_ = oob_score

    def _validate_y(self, y):
        y = column_or_1d(y, warn=True)
        check_classification_targets(y)
        self.classes_, y = np.unique(y, return_inverse=True)
        self.n_classes_ = len(self.classes_)

        return y

    def predict(self, X, **params):
        """Predict class for X.

        The predicted class of an input sample is computed as the class with
        the highest mean predicted probability. If base estimators do not
        implement a ``predict_proba`` method, then it resorts to voting.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The training input samples. Sparse matrices are accepted only if
            they are supported by the base estimator.

        **params : dict
            Parameters routed to the `predict_proba` (if available) or the `predict`
            method (otherwise) of the sub-estimators via the metadata routing API.

            .. versionadded:: 1.7

                Only available if
                `sklearn.set_config(enable_metadata_routing=True)` is set. See
                :ref:`Metadata Routing User Guide <metadata_routing>` for more
                details.

        Returns
        -------
        y : ndarray of shape (n_samples,)
            The predicted classes.
        """
        _raise_for_params(params, self, "predict")

        predicted_probabilitiy = self.predict_proba(X, **params)
        return self.classes_.take((np.argmax(predicted_probabilitiy, axis=1)), axis=0)

    def predict_proba(self, X, **params):
        """Predict class probabilities for X.

        The predicted class probabilities of an input sample is computed as
        the mean predicted class probabilities of the base estimators in the
        ensemble. If base estimators do not implement a ``predict_proba``
        method, then it resorts to voting and the predicted class probabilities
        of an input sample represents the proportion of estimators predicting
        each class.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The training input samples. Sparse matrices are accepted only if
            they are supported by the base estimator.

        **params : dict
            Parameters routed to the `predict_proba` (if available) or the `predict`
            method (otherwise) of the sub-estimators via the metadata routing API.

            .. versionadded:: 1.7

                Only available if
                `sklearn.set_config(enable_metadata_routing=True)` is set. See
                :ref:`Metadata Routing User Guide <metadata_routing>` for more
                details.

        Returns
        -------
        p : ndarray of shape (n_samples, n_classes)
            The class probabilities of the input samples. The order of the
            classes corresponds to that in the attribute :term:`classes_`.
        """
        _raise_for_params(params, self, "predict_proba")

        check_is_fitted(self)
        # Check data
        X = validate_data(
            self,
            X,
            accept_sparse=["csr", "csc"],
            dtype=None,
            ensure_all_finite=False,
            reset=False,
        )

        if _routing_enabled():
            routed_params = process_routing(self, "predict_proba", **params)
        else:
            routed_params = Bunch()
            routed_params.estimator = Bunch(predict_proba=Bunch())

        # Parallel loop
        n_jobs, _, starts = _partition_estimators(self.n_estimators, self.n_jobs)

        all_proba = Parallel(
            n_jobs=n_jobs, verbose=self.verbose, **self._parallel_args()
        )(
            delayed(_parallel_predict_proba)(
                self.estimators_[starts[i] : starts[i + 1]],
                self.estimators_features_[starts[i] : starts[i + 1]],
                X,
                self.n_classes_,
                predict_params=routed_params.estimator.get("predict", None),
                predict_proba_params=routed_params.estimator.get("predict_proba", None),
            )
            for i in range(n_jobs)
        )

        # Reduce
        proba = sum(all_proba) / self.n_estimators

        return proba

    def predict_log_proba(self, X, **params):
        """Predict class log-probabilities for X.

        The predicted class log-probabilities of an input sample is computed as
        the log of the mean predicted class probabilities of the base
        estimators in the ensemble.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The training input samples. Sparse matrices are accepted only if
            they are supported by the base estimator.

        **params : dict
            Parameters routed to the `predict_log_proba`, the `predict_proba` or the
            `proba` method of the sub-estimators via the metadata routing API. The
            routing is tried in the mentioned order depending on whether this method is
            available on the sub-estimator.

            .. versionadded:: 1.7

                Only available if
                `sklearn.set_config(enable_metadata_routing=True)` is set. See
                :ref:`Metadata Routing User Guide <metadata_routing>` for more
                details.

        Returns
        -------
        p : ndarray of shape (n_samples, n_classes)
            The class log-probabilities of the input samples. The order of the
            classes corresponds to that in the attribute :term:`classes_`.
        """
        _raise_for_params(params, self, "predict_log_proba")

        check_is_fitted(self)

        if hasattr(self.estimator_, "predict_log_proba"):
            # Check data
            X = validate_data(
                self,
                X,
                accept_sparse=["csr", "csc"],
                dtype=None,
                ensure_all_finite=False,
                reset=False,
            )

            if _routing_enabled():
                routed_params = process_routing(self, "predict_log_proba", **params)
            else:
                routed_params = Bunch()
                routed_params.estimator = Bunch(predict_log_proba=Bunch())

            # Parallel loop
            n_jobs, _, starts = _partition_estimators(self.n_estimators, self.n_jobs)

            all_log_proba = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
                delayed(_parallel_predict_log_proba)(
                    self.estimators_[starts[i] : starts[i + 1]],
                    self.estimators_features_[starts[i] : starts[i + 1]],
                    X,
                    self.n_classes_,
                    params=routed_params.estimator.predict_log_proba,
                )
                for i in range(n_jobs)
            )

            # Reduce
            log_proba = all_log_proba[0]

            for j in range(1, len(all_log_proba)):
                log_proba = np.logaddexp(log_proba, all_log_proba[j])

            log_proba -= np.log(self.n_estimators)

        else:
            log_proba = np.log(self.predict_proba(X, **params))

        return log_proba

    @available_if(
        _estimator_has("decision_function", delegates=("estimators_", "estimator"))
    )
    def decision_function(self, X, **params):
        """Average of the decision functions of the base classifiers.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The training input samples. Sparse matrices are accepted only if
            they are supported by the base estimator.

        **params : dict
            Parameters routed to the `decision_function` method of the sub-estimators
            via the metadata routing API.

            .. versionadded:: 1.7

                Only available if
                `sklearn.set_config(enable_metadata_routing=True)` is set. See
                :ref:`Metadata Routing User Guide <metadata_routing>` for more
                details.

        Returns
        -------
        score : ndarray of shape (n_samples, k)
            The decision function of the input samples. The columns correspond
            to the classes in sorted order, as they appear in the attribute
            ``classes_``. Regression and binary classification are special
            cases with ``k == 1``, otherwise ``k==n_classes``.
        """
        _raise_for_params(params, self, "decision_function")

        check_is_fitted(self)

        # Check data
        X = validate_data(
            self,
            X,
            accept_sparse=["csr", "csc"],
            dtype=None,
            ensure_all_finite=False,
            reset=False,
        )

        if _routing_enabled():
            routed_params = process_routing(self, "decision_function", **params)
        else:
            routed_params = Bunch()
            routed_params.estimator = Bunch(decision_function=Bunch())

        # Parallel loop
        n_jobs, _, starts = _partition_estimators(self.n_estimators, self.n_jobs)

        all_decisions = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
            delayed(_parallel_decision_function)(
                self.estimators_[starts[i] : starts[i + 1]],
                self.estimators_features_[starts[i] : starts[i + 1]],
                X,
                params=routed_params.estimator.decision_function,
            )
            for i in range(n_jobs)
        )

        # Reduce
        decisions = sum(all_decisions) / self.n_estimators

        return decisions


class BaggingRegressor(RegressorMixin, BaseBagging):
    """A Bagging regressor.

    A Bagging regressor is an ensemble meta-estimator that fits base
    regressors each on random subsets of the original dataset and then
    aggregate their individual predictions (either by voting or by averaging)
    to form a final prediction. Such a meta-estimator can typically be used as
    a way to reduce the variance of a black-box estimator (e.g., a decision
    tree), by introducing randomization into its construction procedure and
    then making an ensemble out of it.

    This algorithm encompasses several works from the literature. When random
    subsets of the dataset are drawn as random subsets of the samples, then
    this algorithm is known as Pasting [1]_. If samples are drawn with
    replacement, then the method is known as Bagging [2]_. When random subsets
    of the dataset are drawn as random subsets of the features, then the method
    is known as Random Subspaces [3]_. Finally, when base estimators are built
    on subsets of both samples and features, then the method is known as
    Random Patches [4]_.

    Read more in the :ref:`User Guide <bagging>`.

    .. versionadded:: 0.15

    Parameters
    ----------
    estimator : object, default=None
        The base estimator to fit on random subsets of the dataset.
        If None, then the base estimator is a
        :class:`~sklearn.tree.DecisionTreeRegressor`.

        .. versionadded:: 1.2
           `base_estimator` was renamed to `estimator`.

    n_estimators : int, default=10
        The number of base estimators in the ensemble.

    max_samples : int or float, default=1.0
        The number of samples to draw from X to train each base estimator (with
        replacement by default, see `bootstrap` for more details).

        - If int, then draw `max_samples` samples.
        - If float, then draw `max_samples * X.shape[0]` samples.

    max_features : int or float, default=1.0
        The number of features to draw from X to train each base estimator (
        without replacement by default, see `bootstrap_features` for more
        details).

        - If int, then draw `max_features` features.
        - If float, then draw `max(1, int(max_features * n_features_in_))` features.

    bootstrap : bool, default=True
        Whether samples are drawn with replacement. If False, sampling without
        replacement is performed. If fitting with `sample_weight`, it is
        strongly recommended to choose True, as only drawing with replacement
        will ensure the expected frequency semantics of `sample_weight`.

    bootstrap_features : bool, default=False
        Whether features are drawn with replacement.

    oob_score : bool, default=False
        Whether to use out-of-bag samples to estimate
        the generalization error. Only available if bootstrap=True.

    warm_start : bool, default=False
        When set to True, reuse the solution of the previous call to fit
        and add more estimators to the ensemble, otherwise, just fit
        a whole new ensemble. See :term:`the Glossary <warm_start>`.

    n_jobs : int, default=None
        The number of jobs to run in parallel for both :meth:`fit` and
        :meth:`predict`. ``None`` means 1 unless in a
        :obj:`joblib.parallel_backend` context. ``-1`` means using all
        processors. See :term:`Glossary <n_jobs>` for more details.

    random_state : int, RandomState instance or None, default=None
        Controls the random resampling of the original dataset
        (sample wise and feature wise).
        If the base estimator accepts a `random_state` attribute, a different
        seed is generated for each instance in the ensemble.
        Pass an int for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    verbose : int, default=0
        Controls the verbosity when fitting and predicting.

    Attributes
    ----------
    estimator_ : estimator
        The base estimator from which the ensemble is grown.

        .. versionadded:: 1.2
           `base_estimator_` was renamed to `estimator_`.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    estimators_ : list of estimators
        The collection of fitted sub-estimators.

    estimators_samples_ : list of arrays
        The subset of drawn samples (i.e., the in-bag samples) for each base
        estimator. Each subset is defined by an array of the indices selected.

    estimators_features_ : list of arrays
        The subset of drawn features for each base estimator.

    oob_score_ : float
        Score of the training dataset obtained using an out-of-bag estimate.
        This attribute exists only when ``oob_score`` is True.

    oob_prediction_ : ndarray of shape (n_samples,)
        Prediction computed with out-of-bag estimate on the training
        set. If n_estimators is small it might be possible that a data point
        was never left out during the bootstrap. In this case,
        `oob_prediction_` might contain NaN. This attribute exists only
        when ``oob_score`` is True.

    See Also
    --------
    BaggingClassifier : A Bagging classifier.

    References
    ----------

    .. [1] L. Breiman, "Pasting small votes for classification in large
           databases and on-line", Machine Learning, 36(1), 85-103, 1999.

    .. [2] L. Breiman, "Bagging predictors", Machine Learning, 24(2), 123-140,
           1996.

    .. [3] T. Ho, "The random subspace method for constructing decision
           forests", Pattern Analysis and Machine Intelligence, 20(8), 832-844,
           1998.

    .. [4] G. Louppe and P. Geurts, "Ensembles on Random Patches", Machine
           Learning and Knowledge Discovery in Databases, 346-361, 2012.

    Examples
    --------
    >>> from sklearn.svm import SVR
    >>> from sklearn.ensemble import BaggingRegressor
    >>> from sklearn.datasets import make_regression
    >>> X, y = make_regression(n_samples=100, n_features=4,
    ...                        n_informative=2, n_targets=1,
    ...                        random_state=0, shuffle=False)
    >>> regr = BaggingRegressor(estimator=SVR(),
    ...                         n_estimators=10, random_state=0).fit(X, y)
    >>> regr.predict([[0, 0, 0, 0]])
    array([-2.8720])
    """

    def __init__(
        self,
        estimator=None,
        n_estimators=10,
        *,
        max_samples=1.0,
        max_features=1.0,
        bootstrap=True,
        bootstrap_features=False,
        oob_score=False,
        warm_start=False,
        n_jobs=None,
        random_state=None,
        verbose=0,
    ):
        super().__init__(
            estimator=estimator,
            n_estimators=n_estimators,
            max_samples=max_samples,
            max_features=max_features,
            bootstrap=bootstrap,
            bootstrap_features=bootstrap_features,
            oob_score=oob_score,
            warm_start=warm_start,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose,
        )

    def predict(self, X, **params):
        """Predict regression target for X.

        The predicted regression target of an input sample is computed as the
        mean predicted regression targets of the estimators in the ensemble.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The training input samples. Sparse matrices are accepted only if
            they are supported by the base estimator.

        **params : dict
            Parameters routed to the `predict` method of the sub-estimators via the
            metadata routing API.

            .. versionadded:: 1.7

                Only available if
                `sklearn.set_config(enable_metadata_routing=True)` is set. See
                :ref:`Metadata Routing User Guide <metadata_routing>` for more
                details.

        Returns
        -------
        y : ndarray of shape (n_samples,)
            The predicted values.
        """
        _raise_for_params(params, self, "predict")

        check_is_fitted(self)
        # Check data
        X = validate_data(
            self,
            X,
            accept_sparse=["csr", "csc"],
            dtype=None,
            ensure_all_finite=False,
            reset=False,
        )

        if _routing_enabled():
            routed_params = process_routing(self, "predict", **params)
        else:
            routed_params = Bunch()
            routed_params.estimator = Bunch(predict=Bunch())

        # Parallel loop
        n_jobs, _, starts = _partition_estimators(self.n_estimators, self.n_jobs)

        all_y_hat = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
            delayed(_parallel_predict_regression)(
                self.estimators_[starts[i] : starts[i + 1]],
                self.estimators_features_[starts[i] : starts[i + 1]],
                X,
                params=routed_params.estimator.predict,
            )
            for i in range(n_jobs)
        )

        # Reduce
        y_hat = sum(all_y_hat) / self.n_estimators

        return y_hat

    def _set_oob_score(self, X, y):
        n_samples = y.shape[0]

        predictions = np.zeros((n_samples,))
        n_predictions = np.zeros((n_samples,))

        for estimator, samples, features in zip(
            self.estimators_, self.estimators_samples_, self.estimators_features_
        ):
            # Create mask for OOB samples
            mask = ~indices_to_mask(samples, n_samples)

            predictions[mask] += estimator.predict((X[mask, :])[:, features])
            n_predictions[mask] += 1

        if (n_predictions == 0).any():
            warn(
                "Some inputs do not have OOB scores. "
                "This probably means too few estimators were used "
                "to compute any reliable oob estimates."
            )
            n_predictions[n_predictions == 0] = 1

        predictions /= n_predictions

        self.oob_prediction_ = predictions
        self.oob_score_ = r2_score(y, predictions)

    def _get_estimator(self):
        """Resolve which estimator to return (default is DecisionTreeClassifier)"""
        if self.estimator is None:
            return DecisionTreeRegressor()
        return self.estimator
