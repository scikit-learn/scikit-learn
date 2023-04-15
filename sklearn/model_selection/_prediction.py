from inspect import signature
from numbers import Integral, Real

import numpy as np

from ..base import BaseEstimator, ClassifierMixin, MetaEstimatorMixin, clone
from ..metrics import check_scoring, get_scorer_names, make_scorer, roc_curve
from ..metrics._scorer import _ContinuousScorer
from ..utils import _safe_indexing
from ..utils._param_validation import HasMethods, Interval, StrOptions
from ..utils._response import _get_response_values_binary
from ..utils.metaestimators import available_if
from ..utils.multiclass import type_of_target
from ..utils.parallel import Parallel, delayed
from ..utils.validation import (
    _check_sample_weight,
    _num_samples,
    check_is_fitted,
    indexable,
)

from ._split import check_cv


def _estimator_has(attr):
    """Check if we can delegate a method to the underlying estimator.

    First, we check the first fitted estimator if available, otherwise we
    check the unfitted estimator.
    """
    return lambda self: (
        hasattr(self.estimator_, attr)
        if hasattr(self, "estimator_")
        else hasattr(self.estimator, attr)
    )


def _fit_and_score(
    classifier, X, y, sample_weight, train_idx, val_idx, scorer, score_method
):
    fit_parameters = signature(classifier.fit).parameters
    supports_sw = "sample_weight" in fit_parameters

    if train_idx is not None:
        X_train, X_val = _safe_indexing(X, train_idx), _safe_indexing(X, val_idx)
        y_train, y_val = _safe_indexing(y, train_idx), _safe_indexing(y, val_idx)
        if sample_weight is not None:
            sw_train, sw_val = (
                _safe_indexing(sample_weight, train_idx),
                _safe_indexing(sample_weight, val_idx),
            )
        else:
            sw_train, sw_val = None, None
        if supports_sw:
            classifier.fit(X_train, y_train, sample_weight=sw_train)
        else:
            classifier.fit(X_train, y_train)
    else:  # prefit estimator, only a validation set is provided
        X_val, y_val, sw_val = X, y, sample_weight

    if score_method == {"tnr", "tpr"}:
        fpr, tpr, potential_thresholds = scorer(
            classifier, X_val, y_val, sample_weight=sw_val
        )
        if score_method == "tnr":
            return potential_thresholds[::-1], (1 - fpr)[::-1]
        return potential_thresholds, tpr
    return scorer(classifier, X_val, y_val, sample_weight=sw_val)


class CutOffClassifier(ClassifierMixin, MetaEstimatorMixin, BaseEstimator):
    """Decision threshold calibration for binary classification.

    Parameters
    ----------
    estimator : estimator instance
        The classifier, fitted or not fitted, for which we want to optimize
        the decision threshold used during `predict`.

    objective_metric : {"tpr", "tnr"}, str or callable,  \
            default="balanced_accuracy"
        The objective metric to be optimized. Can be one of:

        * a string associated to a scoring function (see model evaluation
          documentation);
        * a scorer callable object / function with the signature
          `metric(estimator, X, y)`;
        * `"tpr"`: find the decision threshold for a true positive ratio (TPR)
          of `objective_value`;
        * `"tnr"`: find the decision threshold for a true negative ratio (TNR)
          of `objective_value`.

    objective_value : float, default=None
        The value associated with the `objective_metric` metric for which we
        want to find the decision threshold when `objective_metric` is equal to
        `"tpr"` or `"tnr"`.

    response_method : {"auto", "decision_function", "predict_proba"}, \
            default="auto"
        Methods by the classifier `base_estimator` corresponding to the
        decision function for which we want to find a threshold. It can be:

        * if `"auto"`, it will try to invoke, for each classifier,
          `"predict_proba"` or `"decision_function"` in that order.
        * otherwise, one of `"predict_proba"` or `"decision_function"`.
          If the method is not implemented by the classifier, it will raise an
          error.

    n_thresholds : int, default=1000
        The number of decision threshold to use when discretizing the output
        of the classifier `method`.

    cv : int, float, cross-validation generator, iterable or "prefit", \
            default=None
        Determines the cross-validation splitting strategy to train classifier.
        Possible inputs for cv are:

        * None, to use the default 5-fold stratified K-fold cross validation;
        * An integer number, to specify the number of folds in a stratified
          k-fold;
        * A float number, to specify a single shuffle split. The floating
          number should be in (0, 1) and represent the size of the validation
          set;
        * An object to be used as a cross-validation generator;
        * An iterable yielding train, test splits;
        * "prefit", to bypass the cross-validation.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    n_jobs : int, default=None
        The number of jobs to run in parallel. When `cv` represents a
        cross-validation strategy, the fitting and scoring on each data split
        is done in parallel. ``None`` means 1 unless in a
        :obj:`joblib.parallel_backend` context. ``-1`` means using all
        processors. See :term:`Glossary <n_jobs>` for more details.

    Attributes
    ----------
    decision_threshold_ : float
        The new decision threshold.

    classes_ : ndarray of shape (n_classes,)
        The class labels.

    n_features_in_ : int
        Number of features seen during :term:`fit`. Only defined if the
        underlying estimator exposes such an attribute when fit.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Only defined if the
        underlying estimator exposes such an attribute when fit.
    """

    def __init__(
        self,
        estimator,
        *,
        objective_metric="balanced_accuracy",
        objective_value=None,
        response_method="auto",
        n_thresholds=1_000,
        cv=None,
        n_jobs=None,
    ):
        self.estimator = estimator
        self.objective_metric = objective_metric
        self.objective_value = objective_value
        self.response_method = response_method
        self.n_thresholds = n_thresholds
        self.cv = cv
        self.n_jobs = n_jobs

    _parameter_constraints: dict = {
        "estimator": [
            HasMethods(["fit", "predict_proba"]),
            HasMethods(["fit", "decision_function"]),
        ],
        "objective_metric": [
            StrOptions(set(get_scorer_names()) | {"tpr", "fpr"}),
            callable,
        ],
        "objective_value": [Real, None],
        "response_method": [StrOptions({"auto", "predict_proba", "decision_function"})],
        "n_thresholds": [Interval(Integral, 1, None, closed="left")],
        "cv": ["cv_object", StrOptions({"prefit"})],
        "n_jobs": [Integral, None],
    }

    def fit(self, X, y, sample_weight=None, **fit_params):
        """Fit the calibrated model.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.

        y : array-like of shape (n_samples,)
            Target values.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights. If `None`, then samples are equally weighted.

        **fit_params : dict
            Parameters to pass to the `fit` method of the underlying
            classifier.

        Returns
        -------
        self : object
            Returns an instance of self.
        """
        self._validate_params()
        X, y = indexable(X, y)

        y_type = type_of_target(y, input_name="y")
        if y_type != "binary":
            raise ValueError(
                f"Only binary classification is supported. Unknown label type: {y_type}"
            )

        if sample_weight is not None:
            sample_weight = _check_sample_weight(sample_weight, X)

        cv = self.cv if self.cv == "prefit" else check_cv(self.cv, y, classifier=True)

        if self.response_method == "auto":
            self._response_method = ["predict_proba", "decision_function"]
        else:
            self._response_method = self.response_method

        if self.objective_metric in {"tpr", "fpr"}:
            if self.objective_value is None:
                raise ValueError(
                    "When `objective_metric` is 'tpr' or 'fpr', `objective_value` must "
                    "be provided. Got None instead."
                )
            objective_value = self.objective_value
        else:
            objective_value = "highest"

        fit_parameters = signature(self.estimator.fit).parameters
        supports_sw = "sample_weight" in fit_parameters
        if sample_weight is not None and supports_sw:
            self.estimator_ = clone(self.estimator).fit(
                X, y, sample_weight, **fit_params
            )
        else:
            self.estimator_ = clone(self.estimator).fit(X, y, **fit_params)

        if cv == "prefit":
            classifier = self.estimator
            split = ([None, range(_num_samples(X))],)
        else:
            classifier = clone(self.estimator)
            split = cv.split(X, y)

        if self.objective_metric in {"tpr", "fpr"}:
            self._scorer = make_scorer(roc_curve, needs_threshold=True)
        else:
            scoring = check_scoring(classifier, scoring=self.objective_metric)
            # transform a binary metric into a curve metric for all possible decision
            # thresholds
            self._scorer = _ContinuousScorer(
                score_func=scoring._score_func,
                sign=scoring._sign,
                response_method=self._response_method,
                kwargs=scoring._kwargs,
            )

        thresholds, scores = zip(
            *Parallel(n_jobs=self.n_jobs)(
                delayed(_fit_and_score)(
                    classifier,
                    X,
                    y,
                    sample_weight,
                    train_idx,
                    val_idx,
                    self._scorer,
                    self.objective_metric,
                )
                for train_idx, val_idx in split
            )
        )

        min_threshold = np.min([th.min() for th in thresholds])
        max_threshold = np.max([th.max() for th in thresholds])
        ascending = thresholds[0].argmin() == 0
        start = min_threshold if ascending else max_threshold
        stop = max_threshold if ascending else min_threshold
        thresholds_interpolated = np.linspace(start, stop, num=self.n_thresholds)
        mean_score = np.mean(
            [
                np.interp(thresholds_interpolated, th, sc)
                for th, sc in zip(thresholds, scores)
            ],
            axis=0,
        )
        if objective_value == "highest":
            best_idx = mean_score.argmax()
        else:
            best_idx = np.searchsorted(mean_score, objective_value)
        self.decision_threshold_ = thresholds_interpolated[best_idx]

        if hasattr(self.estimator_, "n_features_in_"):
            self.n_features_in_ = self.estimator_.n_features_in_
        if hasattr(self.estimator_, "feature_names_in_"):
            self.feature_names_in_ = self.estimator_.feature_names_in_

        return self

    @property
    def classes_(self):
        """Classes labels."""
        return self.estimator_.classes_

    def predict(self, X):
        """Predict the target of new samples.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The samples, as accepted by `estimator.predict`.

        Returns
        -------
        C : ndarray of shape (n_samples,)
            The predicted class.
        """
        check_is_fitted(self, "estimator_")
        pos_label = self._scorer._get_pos_label()
        y_score, _ = _get_response_values_binary(
            self.estimator_, X, self._response_method, pos_label=pos_label
        )
        y_pred = (y_score >= self.decision_threshold_).astype(int)
        return self.classes_[y_pred]

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
        check_is_fitted(self, "estimator_")
        return self.estimator_.predict_proba(X)

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
        check_is_fitted(self, "estimator_")
        return self.estimator_.predict_log_proba(X)

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
        check_is_fitted(self, "estimator_")
        return self.estimator_.decision_function(X)

    def _more_tags(self):
        return {"binary_only": True}
