# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from collections.abc import MutableMapping
from numbers import Integral, Real

import numpy as np

from ..base import (
    BaseEstimator,
    ClassifierMixin,
    MetaEstimatorMixin,
    _fit_context,
    clone,
)
from ..exceptions import NotFittedError
from ..metrics import (
    check_scoring,
    get_scorer_names,
)
from ..metrics._scorer import (
    _CurveScorer,
    _threshold_scores_to_class_labels,
)
from ..utils import _safe_indexing
from ..utils._param_validation import HasMethods, Interval, RealNotInt, StrOptions
from ..utils._response import _get_response_values_binary
from ..utils.metadata_routing import (
    MetadataRouter,
    MethodMapping,
    _raise_for_params,
    process_routing,
)
from ..utils.metaestimators import available_if
from ..utils.multiclass import type_of_target
from ..utils.parallel import Parallel, delayed
from ..utils.validation import (
    _check_method_params,
    _estimator_has,
    _num_samples,
    check_is_fitted,
    indexable,
)
from ._split import StratifiedShuffleSplit, check_cv


def _check_is_fitted(estimator):
    try:
        check_is_fitted(estimator.estimator)
    except NotFittedError:
        check_is_fitted(estimator, "estimator_")


class BaseThresholdClassifier(ClassifierMixin, MetaEstimatorMixin, BaseEstimator):
    """Base class for binary classifiers that set a non-default decision threshold.

    In this base class, we define the following interface:

    - the validation of common parameters in `fit`;
    - the different prediction methods that can be used with the classifier.

    .. versionadded:: 1.5

    Parameters
    ----------
    estimator : estimator instance
        The binary classifier, fitted or not, for which we want to optimize
        the decision threshold used during `predict`.

    response_method : {"auto", "decision_function", "predict_proba"}, default="auto"
        Methods by the classifier `estimator` corresponding to the
        decision function for which we want to find a threshold. It can be:

        * if `"auto"`, it will try to invoke, for each classifier,
          `"predict_proba"` or `"decision_function"` in that order.
        * otherwise, one of `"predict_proba"` or `"decision_function"`.
          If the method is not implemented by the classifier, it will raise an
          error.
    """

    _parameter_constraints: dict = {
        "estimator": [
            HasMethods(["fit", "predict_proba"]),
            HasMethods(["fit", "decision_function"]),
        ],
        "response_method": [StrOptions({"auto", "predict_proba", "decision_function"})],
    }

    def __init__(self, estimator, *, response_method="auto"):
        self.estimator = estimator
        self.response_method = response_method

    def _get_response_method(self):
        """Define the response method."""
        if self.response_method == "auto":
            response_method = ["predict_proba", "decision_function"]
        else:
            response_method = self.response_method
        return response_method

    @_fit_context(
        # *ThresholdClassifier*.estimator is not validated yet
        prefer_skip_nested_validation=False
    )
    def fit(self, X, y, **params):
        """Fit the classifier.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.

        y : array-like of shape (n_samples,)
            Target values.

        **params : dict
            Parameters to pass to the `fit` method of the underlying
            classifier.

        Returns
        -------
        self : object
            Returns an instance of self.
        """
        _raise_for_params(params, self, None)

        X, y = indexable(X, y)

        y_type = type_of_target(y, input_name="y")
        if y_type != "binary":
            raise ValueError(
                f"Only binary classification is supported. Unknown label type: {y_type}"
            )

        self._fit(X, y, **params)

        if hasattr(self.estimator_, "n_features_in_"):
            self.n_features_in_ = self.estimator_.n_features_in_
        if hasattr(self.estimator_, "feature_names_in_"):
            self.feature_names_in_ = self.estimator_.feature_names_in_

        return self

    @property
    def classes_(self):
        """Classes labels."""
        return self.estimator_.classes_

    @available_if(_estimator_has("predict_proba"))
    def predict_proba(self, X):
        """Predict class probabilities for `X` using the fitted estimator.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        Returns
        -------
        probabilities : ndarray of shape (n_samples, n_classes)
            The class probabilities of the input samples.
        """
        _check_is_fitted(self)
        estimator = getattr(self, "estimator_", self.estimator)
        return estimator.predict_proba(X)

    @available_if(_estimator_has("predict_log_proba"))
    def predict_log_proba(self, X):
        """Predict logarithm class probabilities for `X` using the fitted estimator.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        Returns
        -------
        log_probabilities : ndarray of shape (n_samples, n_classes)
            The logarithm class probabilities of the input samples.
        """
        _check_is_fitted(self)
        estimator = getattr(self, "estimator_", self.estimator)
        return estimator.predict_log_proba(X)

    @available_if(_estimator_has("decision_function"))
    def decision_function(self, X):
        """Decision function for samples in `X` using the fitted estimator.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vectors, where `n_samples` is the number of samples and
            `n_features` is the number of features.

        Returns
        -------
        decisions : ndarray of shape (n_samples,)
            The decision function computed the fitted estimator.
        """
        _check_is_fitted(self)
        estimator = getattr(self, "estimator_", self.estimator)
        return estimator.decision_function(X)

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.classifier_tags.multi_class = False
        return tags


class FixedThresholdClassifier(BaseThresholdClassifier):
    """Binary classifier that manually sets the decision threshold.

    This classifier allows to change the default decision threshold used for
    converting posterior probability estimates (i.e. output of `predict_proba`) or
    decision scores (i.e. output of `decision_function`) into a class label.

    Here, the threshold is not optimized and is set to a constant value.

    Read more in the :ref:`User Guide <FixedThresholdClassifier>`.

    .. versionadded:: 1.5

    Parameters
    ----------
    estimator : estimator instance
        The binary classifier, fitted or not, for which we want to optimize
        the decision threshold used during `predict`.

    threshold : {"auto"} or float, default="auto"
        The decision threshold to use when converting posterior probability estimates
        (i.e. output of `predict_proba`) or decision scores (i.e. output of
        `decision_function`) into a class label. When `"auto"`, the threshold is set
        to 0.5 if `predict_proba` is used as `response_method`, otherwise it is set to
        0 (i.e. the default threshold for `decision_function`).

    pos_label : int, float, bool or str, default=None
        The label of the positive class. Used to process the output of the
        `response_method` method. When `pos_label=None`, if `y_true` is in `{-1, 1}` or
        `{0, 1}`, `pos_label` is set to 1, otherwise an error will be raised.

    response_method : {"auto", "decision_function", "predict_proba"}, default="auto"
        Methods by the classifier `estimator` corresponding to the
        decision function for which we want to find a threshold. It can be:

        * if `"auto"`, it will try to invoke `"predict_proba"` or `"decision_function"`
          in that order.
        * otherwise, one of `"predict_proba"` or `"decision_function"`.
          If the method is not implemented by the classifier, it will raise an
          error.

    Attributes
    ----------
    estimator_ : estimator instance
        The fitted classifier used when predicting.

    classes_ : ndarray of shape (n_classes,)
        The class labels.

    n_features_in_ : int
        Number of features seen during :term:`fit`. Only defined if the
        underlying estimator exposes such an attribute when fit.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Only defined if the
        underlying estimator exposes such an attribute when fit.

    See Also
    --------
    sklearn.model_selection.TunedThresholdClassifierCV : Classifier that post-tunes
        the decision threshold based on some metrics and using cross-validation.
    sklearn.calibration.CalibratedClassifierCV : Estimator that calibrates
        probabilities.

    Examples
    --------
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.linear_model import LogisticRegression
    >>> from sklearn.metrics import confusion_matrix
    >>> from sklearn.model_selection import FixedThresholdClassifier, train_test_split
    >>> X, y = make_classification(
    ...     n_samples=1_000, weights=[0.9, 0.1], class_sep=0.8, random_state=42
    ... )
    >>> X_train, X_test, y_train, y_test = train_test_split(
    ...     X, y, stratify=y, random_state=42
    ... )
    >>> classifier = LogisticRegression(random_state=0).fit(X_train, y_train)
    >>> print(confusion_matrix(y_test, classifier.predict(X_test)))
    [[217   7]
     [ 19   7]]
    >>> classifier_other_threshold = FixedThresholdClassifier(
    ...     classifier, threshold=0.1, response_method="predict_proba"
    ... ).fit(X_train, y_train)
    >>> print(confusion_matrix(y_test, classifier_other_threshold.predict(X_test)))
    [[184  40]
     [  6  20]]
    """

    _parameter_constraints: dict = {
        **BaseThresholdClassifier._parameter_constraints,
        "threshold": [StrOptions({"auto"}), Real],
        "pos_label": [Real, str, "boolean", None],
    }

    def __init__(
        self,
        estimator,
        *,
        threshold="auto",
        pos_label=None,
        response_method="auto",
    ):
        super().__init__(estimator=estimator, response_method=response_method)
        self.pos_label = pos_label
        self.threshold = threshold

    @property
    def classes_(self):
        if estimator := getattr(self, "estimator_", None):
            return estimator.classes_
        try:
            check_is_fitted(self.estimator)
            return self.estimator.classes_
        except NotFittedError:
            raise AttributeError(
                "The underlying estimator is not fitted yet."
            ) from NotFittedError

    def _fit(self, X, y, **params):
        """Fit the classifier.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.

        y : array-like of shape (n_samples,)
            Target values.

        **params : dict
            Parameters to pass to the `fit` method of the underlying
            classifier.

        Returns
        -------
        self : object
            Returns an instance of self.
        """
        routed_params = process_routing(self, "fit", **params)
        self.estimator_ = clone(self.estimator).fit(X, y, **routed_params.estimator.fit)
        return self

    def predict(self, X):
        """Predict the target of new samples.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The samples, as accepted by `estimator.predict`.

        Returns
        -------
        class_labels : ndarray of shape (n_samples,)
            The predicted class.
        """
        _check_is_fitted(self)

        estimator = getattr(self, "estimator_", self.estimator)

        y_score, _, response_method_used = _get_response_values_binary(
            estimator,
            X,
            self._get_response_method(),
            pos_label=self.pos_label,
            return_response_method_used=True,
        )

        if self.threshold == "auto":
            decision_threshold = 0.5 if response_method_used == "predict_proba" else 0.0
        else:
            decision_threshold = self.threshold

        return _threshold_scores_to_class_labels(
            y_score, decision_threshold, self.classes_, self.pos_label
        )

    def get_metadata_routing(self):
        """Get metadata routing of this object.

        Please check :ref:`User Guide <metadata_routing>` on how the routing
        mechanism works.

        Returns
        -------
        routing : MetadataRouter
            A :class:`~sklearn.utils.metadata_routing.MetadataRouter` encapsulating
            routing information.
        """
        router = MetadataRouter(owner=self.__class__.__name__).add(
            estimator=self.estimator,
            method_mapping=MethodMapping().add(callee="fit", caller="fit"),
        )
        return router


def _fit_and_score_over_thresholds(
    classifier,
    X,
    y,
    *,
    fit_params,
    train_idx,
    val_idx,
    curve_scorer,
    score_params,
):
    """Fit a classifier and compute the scores for different decision thresholds.

    Parameters
    ----------
    classifier : estimator instance
        The classifier to fit and use for scoring. If `classifier` is already fitted,
        it will be used as is.

    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        The entire dataset.

    y : array-like of shape (n_samples,)
        The entire target vector.

    fit_params : dict
        Parameters to pass to the `fit` method of the underlying classifier.

    train_idx : ndarray of shape (n_train_samples,) or None
        The indices of the training set. If `None`, `classifier` is expected to be
        already fitted.

    val_idx : ndarray of shape (n_val_samples,)
        The indices of the validation set used to score `classifier`. If `train_idx`,
        the entire set will be used.

    curve_scorer : scorer instance
        The scorer taking `classifier` and the validation set as input and outputting
        decision thresholds and scores as a curve. Note that this is different from
        the usual scorer that output a single score value:

        * when `score_method` is one of the four constraint metrics, the curve scorer
          will output a curve of two scores parametrized by the decision threshold, e.g.
          TPR/TNR or precision/recall curves for each threshold;
        * otherwise, the curve scorer will output a single score value for each
          threshold.

    score_params : dict
        Parameters to pass to the `score` method of the underlying scorer.

    Returns
    -------
    scores : ndarray of shape (thresholds,) or tuple of such arrays
        The scores computed for each decision threshold. When TPR/TNR or precision/
        recall are computed, `scores` is a tuple of two arrays.

    potential_thresholds : ndarray of shape (thresholds,)
        The decision thresholds used to compute the scores. They are returned in
        ascending order.
    """

    if train_idx is not None:
        X_train, X_val = _safe_indexing(X, train_idx), _safe_indexing(X, val_idx)
        y_train, y_val = _safe_indexing(y, train_idx), _safe_indexing(y, val_idx)
        fit_params_train = _check_method_params(X, fit_params, indices=train_idx)
        score_params_val = _check_method_params(X, score_params, indices=val_idx)
        classifier.fit(X_train, y_train, **fit_params_train)
    else:  # prefit estimator, only a validation set is provided
        X_val, y_val, score_params_val = X, y, score_params

    return curve_scorer(classifier, X_val, y_val, **score_params_val)


def _mean_interpolated_score(target_thresholds, cv_thresholds, cv_scores):
    """Compute the mean interpolated score across folds by defining common thresholds.

    Parameters
    ----------
    target_thresholds : ndarray of shape (thresholds,)
        The thresholds to use to compute the mean score.

    cv_thresholds : ndarray of shape (n_folds, thresholds_fold)
        The thresholds used to compute the scores for each fold.

    cv_scores : ndarray of shape (n_folds, thresholds_fold)
        The scores computed for each threshold for each fold.

    Returns
    -------
    mean_score : ndarray of shape (thresholds,)
        The mean score across all folds for each target threshold.
    """
    return np.mean(
        [
            np.interp(target_thresholds, split_thresholds, split_score)
            for split_thresholds, split_score in zip(cv_thresholds, cv_scores)
        ],
        axis=0,
    )


class TunedThresholdClassifierCV(BaseThresholdClassifier):
    """Classifier that post-tunes the decision threshold using cross-validation.

    This estimator post-tunes the decision threshold (cut-off point) that is
    used for converting posterior probability estimates (i.e. output of
    `predict_proba`) or decision scores (i.e. output of `decision_function`)
    into a class label. The tuning is done by optimizing a binary metric,
    potentially constrained by a another metric.

    Read more in the :ref:`User Guide <TunedThresholdClassifierCV>`.

    .. versionadded:: 1.5

    Parameters
    ----------
    estimator : estimator instance
        The classifier, fitted or not, for which we want to optimize
        the decision threshold used during `predict`.

    scoring : str or callable, default="balanced_accuracy"
        The objective metric to be optimized. Can be one of:

        * a string associated to a scoring function for binary classification
          (see :ref:`scoring_parameter`);
        * a scorer callable object created with :func:`~sklearn.metrics.make_scorer`;

    response_method : {"auto", "decision_function", "predict_proba"}, default="auto"
        Methods by the classifier `estimator` corresponding to the
        decision function for which we want to find a threshold. It can be:

        * if `"auto"`, it will try to invoke, for each classifier,
          `"predict_proba"` or `"decision_function"` in that order.
        * otherwise, one of `"predict_proba"` or `"decision_function"`.
          If the method is not implemented by the classifier, it will raise an
          error.

    thresholds : int or array-like, default=100
        The number of decision threshold to use when discretizing the output of the
        classifier `method`. Pass an array-like to manually specify the thresholds
        to use.

    cv : int, float, cross-validation generator, iterable or "prefit", default=None
        Determines the cross-validation splitting strategy to train classifier.
        Possible inputs for cv are:

        * `None`, to use the default 5-fold stratified K-fold cross validation;
        * An integer number, to specify the number of folds in a stratified k-fold;
        * A float number, to specify a single shuffle split. The floating number should
          be in (0, 1) and represent the size of the validation set;
        * An object to be used as a cross-validation generator;
        * An iterable yielding train, test splits;
        * `"prefit"`, to bypass the cross-validation.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

        .. warning::
            Using `cv="prefit"` and passing the same dataset for fitting `estimator`
            and tuning the cut-off point is subject to undesired overfitting. You can
            refer to :ref:`TunedThresholdClassifierCV_no_cv` for an example.

            This option should only be used when the set used to fit `estimator` is
            different from the one used to tune the cut-off point (by calling
            :meth:`TunedThresholdClassifierCV.fit`).

    refit : bool, default=True
        Whether or not to refit the classifier on the entire training set once
        the decision threshold has been found.
        Note that forcing `refit=False` on cross-validation having more
        than a single split will raise an error. Similarly, `refit=True` in
        conjunction with `cv="prefit"` will raise an error.

    n_jobs : int, default=None
        The number of jobs to run in parallel. When `cv` represents a
        cross-validation strategy, the fitting and scoring on each data split
        is done in parallel. ``None`` means 1 unless in a
        :obj:`joblib.parallel_backend` context. ``-1`` means using all
        processors. See :term:`Glossary <n_jobs>` for more details.

    random_state : int, RandomState instance or None, default=None
        Controls the randomness of cross-validation when `cv` is a float.
        See :term:`Glossary <random_state>`.

    store_cv_results : bool, default=False
        Whether to store all scores and thresholds computed during the cross-validation
        process.

    Attributes
    ----------
    estimator_ : estimator instance
        The fitted classifier used when predicting.

    best_threshold_ : float
        The new decision threshold.

    best_score_ : float or None
        The optimal score of the objective metric, evaluated at `best_threshold_`.

    cv_results_ : dict or None
        A dictionary containing the scores and thresholds computed during the
        cross-validation process. Only exist if `store_cv_results=True`. The
        keys are `"thresholds"` and `"scores"`.

    classes_ : ndarray of shape (n_classes,)
        The class labels.

    n_features_in_ : int
        Number of features seen during :term:`fit`. Only defined if the
        underlying estimator exposes such an attribute when fit.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Only defined if the
        underlying estimator exposes such an attribute when fit.

    See Also
    --------
    sklearn.model_selection.FixedThresholdClassifier : Classifier that uses a
        constant threshold.
    sklearn.calibration.CalibratedClassifierCV : Estimator that calibrates
        probabilities.

    Examples
    --------
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.ensemble import RandomForestClassifier
    >>> from sklearn.metrics import classification_report
    >>> from sklearn.model_selection import TunedThresholdClassifierCV, train_test_split
    >>> X, y = make_classification(
    ...     n_samples=1_000, weights=[0.9, 0.1], class_sep=0.8, random_state=42
    ... )
    >>> X_train, X_test, y_train, y_test = train_test_split(
    ...     X, y, stratify=y, random_state=42
    ... )
    >>> classifier = RandomForestClassifier(random_state=0).fit(X_train, y_train)
    >>> print(classification_report(y_test, classifier.predict(X_test)))
                  precision    recall  f1-score   support
    <BLANKLINE>
               0       0.94      0.99      0.96       224
               1       0.80      0.46      0.59        26
    <BLANKLINE>
        accuracy                           0.93       250
       macro avg       0.87      0.72      0.77       250
    weighted avg       0.93      0.93      0.92       250
    <BLANKLINE>
    >>> classifier_tuned = TunedThresholdClassifierCV(
    ...     classifier, scoring="balanced_accuracy"
    ... ).fit(X_train, y_train)
    >>> print(
    ...     f"Cut-off point found at {classifier_tuned.best_threshold_:.3f}"
    ... )
    Cut-off point found at 0.342
    >>> print(classification_report(y_test, classifier_tuned.predict(X_test)))
                  precision    recall  f1-score   support
    <BLANKLINE>
               0       0.96      0.95      0.96       224
               1       0.61      0.65      0.63        26
    <BLANKLINE>
        accuracy                           0.92       250
       macro avg       0.78      0.80      0.79       250
    weighted avg       0.92      0.92      0.92       250
    <BLANKLINE>
    """

    _parameter_constraints: dict = {
        **BaseThresholdClassifier._parameter_constraints,
        "scoring": [
            StrOptions(set(get_scorer_names())),
            callable,
            MutableMapping,
        ],
        "thresholds": [Interval(Integral, 1, None, closed="left"), "array-like"],
        "cv": [
            "cv_object",
            StrOptions({"prefit"}),
            Interval(RealNotInt, 0.0, 1.0, closed="neither"),
        ],
        "refit": ["boolean"],
        "n_jobs": [Integral, None],
        "random_state": ["random_state"],
        "store_cv_results": ["boolean"],
    }

    def __init__(
        self,
        estimator,
        *,
        scoring="balanced_accuracy",
        response_method="auto",
        thresholds=100,
        cv=None,
        refit=True,
        n_jobs=None,
        random_state=None,
        store_cv_results=False,
    ):
        super().__init__(estimator=estimator, response_method=response_method)
        self.scoring = scoring
        self.thresholds = thresholds
        self.cv = cv
        self.refit = refit
        self.n_jobs = n_jobs
        self.random_state = random_state
        self.store_cv_results = store_cv_results

    def _fit(self, X, y, **params):
        """Fit the classifier and post-tune the decision threshold.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.

        y : array-like of shape (n_samples,)
            Target values.

        **params : dict
            Parameters to pass to the `fit` method of the underlying
            classifier and to the `scoring` scorer.

        Returns
        -------
        self : object
            Returns an instance of self.
        """
        if isinstance(self.cv, Real) and 0 < self.cv < 1:
            cv = StratifiedShuffleSplit(
                n_splits=1, test_size=self.cv, random_state=self.random_state
            )
        elif self.cv == "prefit":
            if self.refit is True:
                raise ValueError("When cv='prefit', refit cannot be True.")
            try:
                check_is_fitted(self.estimator, "classes_")
            except NotFittedError as exc:
                raise NotFittedError(
                    """When cv='prefit', `estimator` must be fitted."""
                ) from exc
            cv = self.cv
        else:
            cv = check_cv(self.cv, y=y, classifier=True)
            if self.refit is False and cv.get_n_splits() > 1:
                raise ValueError("When cv has several folds, refit cannot be False.")

        routed_params = process_routing(self, "fit", **params)
        self._curve_scorer = self._get_curve_scorer()

        # in the following block, we:
        # - define the final classifier `self.estimator_` and train it if necessary
        # - define `classifier` to be used to post-tune the decision threshold
        # - define `split` to be used to fit/score `classifier`
        if cv == "prefit":
            self.estimator_ = self.estimator
            classifier = self.estimator_
            splits = [(None, range(_num_samples(X)))]
        else:
            self.estimator_ = clone(self.estimator)
            classifier = clone(self.estimator)
            splits = cv.split(X, y, **routed_params.splitter.split)

            if self.refit:
                # train on the whole dataset
                X_train, y_train, fit_params_train = X, y, routed_params.estimator.fit
            else:
                # single split cross-validation
                train_idx, _ = next(cv.split(X, y, **routed_params.splitter.split))
                X_train = _safe_indexing(X, train_idx)
                y_train = _safe_indexing(y, train_idx)
                fit_params_train = _check_method_params(
                    X, routed_params.estimator.fit, indices=train_idx
                )

            self.estimator_.fit(X_train, y_train, **fit_params_train)

        cv_scores, cv_thresholds = zip(
            *Parallel(n_jobs=self.n_jobs)(
                delayed(_fit_and_score_over_thresholds)(
                    clone(classifier) if cv != "prefit" else classifier,
                    X,
                    y,
                    fit_params=routed_params.estimator.fit,
                    train_idx=train_idx,
                    val_idx=val_idx,
                    curve_scorer=self._curve_scorer,
                    score_params=routed_params.scorer.score,
                )
                for train_idx, val_idx in splits
            )
        )

        if any(np.isclose(th[0], th[-1]) for th in cv_thresholds):
            raise ValueError(
                "The provided estimator makes constant predictions. Therefore, it is "
                "impossible to optimize the decision threshold."
            )

        # find the global min and max thresholds across all folds
        min_threshold = min(
            split_thresholds.min() for split_thresholds in cv_thresholds
        )
        max_threshold = max(
            split_thresholds.max() for split_thresholds in cv_thresholds
        )
        if isinstance(self.thresholds, Integral):
            decision_thresholds = np.linspace(
                min_threshold, max_threshold, num=self.thresholds
            )
        else:
            decision_thresholds = np.asarray(self.thresholds)

        objective_scores = _mean_interpolated_score(
            decision_thresholds, cv_thresholds, cv_scores
        )
        best_idx = objective_scores.argmax()
        self.best_score_ = objective_scores[best_idx]
        self.best_threshold_ = decision_thresholds[best_idx]
        if self.store_cv_results:
            self.cv_results_ = {
                "thresholds": decision_thresholds,
                "scores": objective_scores,
            }

        return self

    def predict(self, X):
        """Predict the target of new samples.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The samples, as accepted by `estimator.predict`.

        Returns
        -------
        class_labels : ndarray of shape (n_samples,)
            The predicted class.
        """
        check_is_fitted(self, "estimator_")
        pos_label = self._curve_scorer._get_pos_label()
        y_score, _ = _get_response_values_binary(
            self.estimator_,
            X,
            self._get_response_method(),
            pos_label=pos_label,
        )

        return _threshold_scores_to_class_labels(
            y_score, self.best_threshold_, self.classes_, pos_label
        )

    def get_metadata_routing(self):
        """Get metadata routing of this object.

        Please check :ref:`User Guide <metadata_routing>` on how the routing
        mechanism works.

        Returns
        -------
        routing : MetadataRouter
            A :class:`~sklearn.utils.metadata_routing.MetadataRouter` encapsulating
            routing information.
        """
        router = (
            MetadataRouter(owner=self.__class__.__name__)
            .add(
                estimator=self.estimator,
                method_mapping=MethodMapping().add(callee="fit", caller="fit"),
            )
            .add(
                splitter=self.cv,
                method_mapping=MethodMapping().add(callee="split", caller="fit"),
            )
            .add(
                scorer=self._get_curve_scorer(),
                method_mapping=MethodMapping().add(callee="score", caller="fit"),
            )
        )
        return router

    def _get_curve_scorer(self):
        """Get the curve scorer based on the objective metric used."""
        scoring = check_scoring(self.estimator, scoring=self.scoring)
        curve_scorer = _CurveScorer.from_scorer(
            scoring, self._get_response_method(), self.thresholds
        )
        return curve_scorer
