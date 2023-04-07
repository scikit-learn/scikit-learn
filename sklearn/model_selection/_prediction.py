from numbers import Integral, Real

import numpy as np

from ..base import BaseEstimator, ClassifierMixin, MetaEstimatorMixin, clone
from ..metrics import check_scoring, get_scorer_names, make_scorer, roc_curve
from ..metrics._scorer import _ContinuousScorer
from ..utils import _safe_indexing
from ..utils._param_validation import HasMethods, Interval, StrOptions
from ..utils._response import _get_response_values_binary
from ..utils.multiclass import type_of_target
from ..utils.parallel import Parallel, delayed
from ..utils.validation import _check_sample_weight, _num_samples

from ._split import check_cv


def _fit_and_score(classifier, X, y, train_idx, val_idx, scorer, score_method):
    if train_idx is not None:
        X_train, X_val = _safe_indexing(X, train_idx), _safe_indexing(X, val_idx)
        y_train, y_val = _safe_indexing(y, train_idx), _safe_indexing(y, val_idx)
        classifier.fit(X_train, y_train)
    else:  # prefit estimator, only a validation set is provided
        X_val, y_val = X, y

    if score_method == {"tnr", "tpr"}:
        fpr, tpr, potential_thresholds = scorer(classifier, X_val, y_val)
        if score_method == "tnr":
            return potential_thresholds[::-1], (1 - fpr)[::-1]
        return potential_thresholds, tpr
    return scorer(classifier, X_val, y_val)


class CutOffClassifier(ClassifierMixin, MetaEstimatorMixin, BaseEstimator):
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
        self._validate_params()

        y_type = type_of_target(y, input_name="y")
        if y_type != "binary":
            raise ValueError(
                f"Only binary classification is supported. Got {y_type} instead."
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

        self.estimator_ = clone(self.estimator).fit(X, y, sample_weight, **fit_params)

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

        return self

    @property
    def classes_(self):
        """Classes labels."""
        return self.estimator_.classes_

    def predict(self, X):
        pos_label = self._scorer._get_pos_label()
        y_score, _ = _get_response_values_binary(
            self.estimator_, X, self._response_method, pos_label=pos_label
        )
        y_pred = (y_score >= self.decision_threshold_).astype(int)
        return self.classes_[y_pred]
