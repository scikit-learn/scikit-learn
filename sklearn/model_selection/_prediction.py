from collections.abc import MutableMapping
from inspect import signature
from numbers import Integral, Real

import numpy as np

from ..base import BaseEstimator, ClassifierMixin, MetaEstimatorMixin, clone
from ..exceptions import NotFittedError
from ..metrics import (
    check_scoring,
    confusion_matrix,
    get_scorer_names,
    make_scorer,
    precision_recall_curve,
    roc_curve,
)
from ..metrics._scorer import _ContinuousScorer
from ..utils import _safe_indexing
from ..utils._param_validation import HasMethods, Interval, RealNotInt, StrOptions
from ..utils._response import _get_response_values_binary
from ..utils.metaestimators import available_if
from ..utils.multiclass import type_of_target
from ..utils.parallel import Parallel, delayed
from ..utils.validation import (
    _check_fit_params,
    _check_sample_weight,
    _check_pos_label_consistency,
    _num_samples,
    check_consistent_length,
    check_is_fitted,
    indexable,
)

from ._split import check_cv, StratifiedShuffleSplit


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
    classifier,
    X,
    y,
    sample_weight,
    fit_params,
    train_idx,
    val_idx,
    scorer,
    score_method,
):
    """Fit a classifier and compute the scores for different decision thresholds.

    Parameters
    ----------
    classifier : estimator instance
        The classifier to fit and used for scoring. If `classifier` is already fitted,
        it will be used as is.

    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        The entire dataset.

    y : array-like of shape (n_samples,)
        The entire target vector.

    sample_weight : array-like of shape (n_samples,)
        Some optional associated sample weights.

    fit_params : dict
        Parameters to pass to the `fit` method of the underlying classifier.

    train_idx : ndarray of shape (n_train_samples,) or None
        The indices of the training set. If `None`, `classifier` is expected to be
        already fitted.

    val_idx : ndarray of shape (n_val_samples,)
        The indices of the validation set used to score `classifier`. If `train_idx`,
        the entire set will be used.

    scorer : scorer instance
        The scorer taking `classifier` and the validation set as input and outputting
        decision thresholds and scores.

    score_method : str or callable
        The scoring method to use. Used to detect if we compute TPR/TNR or precision/
        recall.

    Returns
    -------
    thresholds : ndarray of shape (n_thresholds,)
        The decision thresholds used to compute the scores. They are returned in
        ascending order.

    scores : ndarray of shape (n_thresholds,) or tuple os such arrays
        The scores computed for each decision threshold. When TPR/TNR or precision/
        recall are computed, `scores` is a tuple of two arrays.
    """
    arrays = (X, y) if sample_weight is None else (X, y, sample_weight)
    check_consistent_length(*arrays)

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
        fit_params_train = _check_fit_params(X, fit_params, indices=train_idx)
        if supports_sw:
            classifier.fit(X_train, y_train, sample_weight=sw_train, **fit_params_train)
        else:
            classifier.fit(X_train, y_train, **fit_params_train)
    else:  # prefit estimator, only a validation set is provided
        X_val, y_val, sw_val = X, y, sample_weight
        check_is_fitted(classifier, "classes_")

    if isinstance(score_method, str):
        if score_method in {"max_tpr_at_tnr_constraint", "max_tnr_at_tpr_constraint"}:
            fpr, tpr, potential_thresholds = scorer(
                classifier, X_val, y_val, sample_weight=sw_val
            )
            # thresholds are in decreasing order
            return potential_thresholds[::-1], ((1 - fpr)[::-1], tpr[::-1])
        elif score_method in {
            "max_precision_at_recall_constraint",
            "max_recall_at_precision_constraint",
        }:
            precision, recall, potential_thresholds = scorer(
                classifier, X_val, y_val, sample_weight=sw_val
            )
            # thresholds are in increasing order, we also have one missing threshold
            # TODO: check what to do with the missing threshold or additional scores.
            return potential_thresholds, (precision[:-1], recall[:-1])
    return scorer(classifier, X_val, y_val, sample_weight=sw_val)


class CutOffClassifier(ClassifierMixin, MetaEstimatorMixin, BaseEstimator):
    """Decision threshold calibration for binary classification.

    Parameters
    ----------
    estimator : estimator instance
        The classifier, fitted or not fitted, for which we want to optimize
        the decision threshold used during `predict`.

    objective_metric : {"max_tpr_at_tnr_constraint", "max_tnr_at_tpr_constraint", \
            "max_precision_at_recall_constraint, "max_recall_at_precision_constraint"} \
            , str, dict or callable, default="balanced_accuracy"
        The objective metric to be optimized. Can be one of:

        * a string associated to a scoring function (see model evaluation
          documentation);
        * a scorer callable object created with :func:`~sklearn.metrics.make_scorer`;
        * `"max_tnr_at_tpr_constraint"`: find the decision threshold for a true
          positive ratio (TPR) of `constraint_value`;
        * `"max_tpr_at_tnr_constraint"`: find the decision threshold for a true
          negative ratio (TNR) of `constraint_value`.
        * `"max_precision_at_recall_constraint"`: find the decision threshold for a
          recall of `constraint_value`;
        * `"max_recall_at_precision_constraint"`: find the decision threshold for a
          precision of `constraint_value`.
        * a dictionary to be used as cost-sensitive matrix. The keys of the
          dictionary should be: `("tp", "fp", "tn", "fn")`. The values of the
          dictionary corresponds costs (negative values) and gains (positive
          values).

    constraint_value : float, default=None
        The value associated with the `objective_metric` metric for which we
        want to find the decision threshold when `objective_metric` is equal one of
        `"max_tnr_at_tpr_constraint"`, `"max_tpr_at_tnr_constraint"`,
        `"max_precision_at_recall_constraint"`, or
        `"max_recall_at_precision_constraint".

    pos_label : int, float, bool or str, default=None
        The label of the positive class. Used when `objective_metric` is
        `"max_tnr_at_tpr_constraint"`"`, `"max_tpr_at_tnr_constraint"`, or a dictionary.
        When `pos_label=None`, if `y_true` is in `{-1, 1}` or `{0, 1}`,
        `pos_label` is set to 1, otherwise an error will be raised. When using a
        scorer, `pos_label` can be passed as a keyword argument to
        :func:`~sklearn.metrics.make_scorer`.

    response_method : {"auto", "decision_function", "predict_proba"}, default="auto"
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

    cv : int, float, cross-validation generator, iterable or "prefit", default=None
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

    refit : "auto" or bool, default="auto"
        Whether or not to refit the classifier on the entire training set once
        the decision threshold has been found. By default, `refit="auto"` is
        equivalent to `refit=False` when `cv` is a float number using a single
        shuffle split or `cv="prefit"` otherwise `refit=True` in all other
        cases. Note that forcing `refit=False` on cross-validation having more
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

    Attributes
    ----------
    decision_threshold_ : float
        The new decision threshold.

    objective_score_ : float
        The score of the objective metric associated with the decision threshold found.

    classes_ : ndarray of shape (n_classes,)
        The class labels.

    n_features_in_ : int
        Number of features seen during :term:`fit`. Only defined if the
        underlying estimator exposes such an attribute when fit.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Only defined if the
        underlying estimator exposes such an attribute when fit.
    """

    _parameter_constraints: dict = {
        "estimator": [
            HasMethods(["fit", "predict_proba"]),
            HasMethods(["fit", "decision_function"]),
        ],
        "objective_metric": [
            StrOptions(
                set(get_scorer_names())
                | {
                    "max_tnr_at_tpr_constraint",
                    "max_tpr_at_tnr_constraint",
                    "max_precision_at_recall_constraint",
                    "max_recall_at_precision_constraint",
                }
            ),
            callable,
            MutableMapping,
        ],
        "constraint_value": [Real, None],
        "pos_label": [Real, str, "boolean", None],
        "response_method": [StrOptions({"auto", "predict_proba", "decision_function"})],
        "n_thresholds": [Interval(Integral, 1, None, closed="left")],
        "cv": [
            "cv_object",
            StrOptions({"prefit"}),
            Interval(RealNotInt, 0.0, 1.0, closed="right"),
        ],
        "refit": ["boolean", StrOptions({"auto"})],
        "n_jobs": [Integral, None],
        "random_state": ["random_state"],
    }

    def __init__(
        self,
        estimator,
        *,
        objective_metric="balanced_accuracy",
        constraint_value=None,
        pos_label=None,
        response_method="auto",
        n_thresholds=1_000,
        cv=None,
        refit="auto",
        n_jobs=None,
        random_state=None,
    ):
        self.estimator = estimator
        self.objective_metric = objective_metric
        self.constraint_value = constraint_value
        self.pos_label = pos_label
        self.response_method = response_method
        self.n_thresholds = n_thresholds
        self.cv = cv
        self.refit = refit
        self.n_jobs = n_jobs
        self.random_state = random_state

    def fit(self, X, y, sample_weight=None, **fit_params):
        """Fit the classifier and post-tune the decision threshold.

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

        if isinstance(self.cv, Real) and 0 < self.cv <= 1:
            cv = StratifiedShuffleSplit(
                n_splits=1, test_size=self.cv, random_state=self.random_state
            )
            refit = False if self.refit == "auto" else self.refit
        elif self.cv == "prefit":
            if self.refit is True:
                raise ValueError("When cv='prefit', refit cannot be True.")
            try:
                check_is_fitted(self.estimator, "classes_")
            except NotFittedError as exc:
                raise NotFittedError(
                    """When cv='prefit', `estimator` must be fitted."""
                ) from exc
            cv, refit = self.cv, False
        else:
            cv = check_cv(self.cv, y=y, classifier=True)
            if self.refit is False and cv.get_n_splits() > 1:
                raise ValueError("When cv has several folds, refit cannot be False.")
            if self.refit == "auto":
                refit = True if cv.get_n_splits() == 1 else False
            else:
                refit = self.refit

        if self.response_method == "auto":
            self._response_method = ["predict_proba", "decision_function"]
        else:
            self._response_method = self.response_method

        if isinstance(self.objective_metric, str) and self.objective_metric in {
            "max_tpr_at_tnr_constraint",
            "max_tnr_at_tpr_constraint",
            "max_precision_at_recall_constraint",
            "max_recall_at_precision_constraint",
        }:
            if self.constraint_value is None:
                raise ValueError(
                    "When `objective_metric` is 'max_tpr_at_tnr_constraint', "
                    "'max_tnr_at_tpr_constraint', 'max_precision_at_recall_constraint',"
                    " or 'max_recall_at_precision_constraint', `constraint_value` must "
                    "be provided. Got None instead."
                )
            constraint_value = self.constraint_value
        else:
            constraint_value = "highest"

        fit_parameters = signature(self.estimator.fit).parameters
        supports_sw = "sample_weight" in fit_parameters

        # in the following block, we:
        # - define the final classifier `self.estimator_` and train it if necessary
        # - define `classifier` to be used to post-tune the decision threshold
        # - define `split` to be used to fit/score `classifier`
        if cv == "prefit":
            self.estimator_ = self.estimator
            classifier = self.estimator_
            splits = ([None, range(_num_samples(X))],)
        else:
            self.estimator_ = clone(self.estimator)
            classifier = clone(self.estimator)
            splits = cv.split(X, y)

            if refit:
                # train on the whole dataset
                X_train, y_train, sw_train = X, y, sample_weight
                fit_params_train = _check_fit_params(X, fit_params, indices=None)
            else:
                # single split cross-validation
                train_idx, _ = next(cv.split(X, y))
                X_train = _safe_indexing(X, train_idx)
                y_train = _safe_indexing(y, train_idx)
                if sample_weight is not None:
                    sw_train = _safe_indexing(sample_weight, train_idx)
                else:
                    sw_train = None
                fit_params_train = _check_fit_params(X, fit_params, indices=train_idx)

            if sw_train is not None and supports_sw:
                self.estimator_.fit(
                    X_train, y_train, sample_weight=sw_train, **fit_params_train
                )
            else:
                self.estimator_.fit(X_train, y_train, **fit_params_train)

        if isinstance(self.objective_metric, MutableMapping):
            keys = set(self.objective_metric.keys())
            if not keys == {"tp", "tn", "fp", "fn"}:
                raise ValueError(
                    "Invalid keys in `objective_metric`. Valid keys are "
                    f"'tp', 'tn', 'fp', and 'fn'. Got {keys} instead."
                )
            pos_label = _check_pos_label_consistency(self.pos_label, y)

            def cost_sensitive_score_func(y_true, y_pred, **kwargs):
                costs_and_gain = np.array(
                    [
                        [kwargs["tn"], kwargs["fp"]],
                        [kwargs["fn"], kwargs["tp"]],
                    ]
                )

                sample_weight = kwargs.get("sample_weight", None)
                cm = confusion_matrix(y_true, y_pred, sample_weight=sample_weight)

                pos_label, classes = kwargs["pos_label"], np.unique(y_true)
                pos_label_idx = np.searchsorted(classes, pos_label)
                if pos_label_idx == 0:
                    # reorder the confusion matrix to be aligned with the cost-matrix
                    cm = cm[::-1, ::-1]

                return (costs_and_gain * cm).sum()

            self._scorer = _ContinuousScorer(
                score_func=cost_sensitive_score_func,
                sign=1,
                response_method=self._response_method,
                kwargs={
                    **self.objective_metric,
                    "pos_label": pos_label,
                },
            )
        elif self.objective_metric in {
            "max_tnr_at_tpr_constraint",
            "max_tpr_at_tnr_constraint",
            "max_precision_at_recall_constraint",
            "max_recall_at_precision_constraint",
        }:
            if self._response_method == "predict_proba":
                params_scorer = {"needs_proba": True, "pos_label": self.pos_label}
            elif (
                isinstance(self._response_method, list)
                and self._response_method[0] == "predict_proba"
                and hasattr(classifier, "predict_proba")
            ):
                # TODO: this is due to a limitation in `make_scorer`: ideally, we should
                # be able to pass a list of response methods to `make_scorer` and give
                # priority to `predict_proba` other `decision_function`.
                # Here, we manually check if the classifier provide `predict_proba` to
                # use `needs_proba` instead and ensure that no error will be raised.
                params_scorer = {"needs_proba": True, "pos_label": self.pos_label}
            else:
                params_scorer = {"needs_threshold": True, "pos_label": self.pos_label}

            if "tpr" in self.objective_metric:  # tpr/tnr
                score_func = roc_curve
            else:  # precision/recall
                score_func = precision_recall_curve
            self._scorer = make_scorer(score_func, **params_scorer)
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
                    fit_params,
                    train_idx,
                    val_idx,
                    self._scorer,
                    self.objective_metric,
                )
                for train_idx, val_idx in splits
            )
        )

        # we add/subtract an arbitrary value to the min/max thresholds to ensure that
        # we get the case where `y_pred` will be all zeros and all ones.
        if hasattr(classifier, "predict_proba") and (
            self._response_method == "predict_proba"
            or (
                isinstance(self._response_method, list)
                and self._response_method[0] == "predict_proba"
            )
        ):
            # `predict_proba` was used to compute scores
            min_threshold = 0.0 - np.finfo(np.float64).eps
            max_threshold = 1.0 + np.finfo(np.float64).eps
        else:
            # `decision_function` was used to compute scores
            min_threshold = np.min([th.min() for th in thresholds]) - 1.0
            max_threshold = np.max([th.max() for th in thresholds]) + 1.0

        # thresholds are sorted in ascending order which is necessary for the
        # interpolation of the score below
        thresholds_interpolated = np.linspace(
            min_threshold, max_threshold, num=self.n_thresholds
        )

        def _mean_interpolated_score(thresholds, scores):
            return np.mean(
                [
                    np.interp(thresholds_interpolated, th, sc)
                    for th, sc in zip(thresholds, scores)
                ],
                axis=0,
            )

        if constraint_value == "highest":  # find best score
            # we don't need to sort the scores and directly take the maximum
            mean_score = _mean_interpolated_score(thresholds, scores)
            best_idx = mean_score.argmax()
            self.objective_score_ = mean_score[best_idx]
            self.decision_threshold_ = thresholds_interpolated[best_idx]
        else:
            if "tpr" in self.objective_metric:  # tpr/tnr
                mean_tnr, mean_tpr = [
                    _mean_interpolated_score(thresholds, sc) for sc in zip(*scores)
                ]
            else:  # precision/recall
                mean_precision, mean_recall = [
                    _mean_interpolated_score(thresholds, sc) for sc in zip(*scores)
                ]

            def _get_best_idx(constrained_score, maximized_score):
                indices = np.arange(len(constrained_score))
                mask = constrained_score >= constraint_value
                mask_idx = maximized_score[mask].argmax()
                return indices[mask][mask_idx]

            if self.objective_metric == "max_tpr_at_tnr_constraint":
                constrained_score, maximized_score = mean_tnr, mean_tpr
            elif self.objective_metric == "max_tnr_at_tpr_constraint":
                constrained_score, maximized_score = mean_tpr, mean_tnr
            elif self.objective_metric == "max_precision_at_recall_constraint":
                constrained_score, maximized_score = mean_recall, mean_precision
            else:  # max_recall_at_precision_constraint
                constrained_score, maximized_score = mean_precision, mean_recall

            best_idx = _get_best_idx(constrained_score, maximized_score)
            self.objective_score_ = (
                constrained_score[best_idx],
                maximized_score[best_idx],
            )
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
        return {
            "binary_only": True,
            "_xfail_checks": {
                "check_classifiers_train": "Threshold at probability 0.5 does not hold"
            },
        }
