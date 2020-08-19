from inspect import signature
import numbers

import numpy as np
from joblib import Parallel, delayed

from ._split import check_cv
from ._split import StratifiedShuffleSplit

from ..base import clone
from ..base import BaseEstimator
from ..base import ClassifierMixin
from ..base import MetaEstimatorMixin
from ..metrics import balanced_accuracy_score
from ..metrics import check_scoring
from ..metrics import make_scorer
from ..metrics import roc_curve
from ..metrics._plot.base import _check_classifier_response_method
from ..metrics._scorer import _BaseScorer
from ..preprocessing import LabelEncoder
from ..utils import check_array
from ..utils import _safe_indexing
from ..utils.multiclass import check_classification_targets
from ..utils.multiclass import type_of_target
from ..utils.validation import check_is_fitted


class _ContinuousScorer(_BaseScorer):
    def __init__(self, score_func, sign, response_method, kwargs):
        super().__init__(score_func=score_func, sign=sign, kwargs=kwargs)
        self.response_method = response_method

    def _score(self, method_caller, estimator, X, y_true, sample_weight=None):
        """Evaluate predicted target values for X relative to y_true.

        Parameters
        ----------
        method_caller : callable
            Returns predictions given an estimator, method name, and other
            arguments, potentially caching results.

        estimator : object
            Trained estimator to use for scoring. Must have a predict_proba
            method; the output of that is used to compute the score.

        X : {array-like, sparse matrix}
            Test data that will be fed to estimator.predict.

        y_true : array-like
            Gold standard target values for X.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        Returns
        -------
        score : float
            Score function applied to prediction of estimator on X.
        """
        response_method = _check_classifier_response_method(
            estimator=estimator, response_method=self.response_method
        )
        y_score = response_method(X)
        if response_method.__name__ == "decision_function":
            y_score = self._check_decision_function(
                y_score, estimator.classes_
            )
        else:
            y_score = self._select_proba(
                y_score, estimator.classes_, support_multi_class=False
            )

        # `np.unique` returned sorted array, thus no need to sort values
        potential_thresholds = np.unique(y_score)
        score_thresholds = []
        for th in potential_thresholds:
            y_score_thresholded = estimator.classes_[
                (y_score >= th).astype(int)
            ]
            if sample_weight is not None:
                score_thresholds.append(
                    self._sign
                    * self._score_func(
                        y_true,
                        y_score_thresholded,
                        sample_weight=sample_weight,
                        **self._kwargs,
                    )
                )
            else:
                score_thresholds.append(
                    self._sign
                    * self._score_func(
                        y_true, y_score_thresholded, **self._kwargs
                    )
                )
        return np.array(potential_thresholds), np.array(score_thresholds)


class CutoffClassifier(MetaEstimatorMixin, ClassifierMixin, BaseEstimator):
    """Decision threshold calibration for binary classification.

    Estimator that calibrates the decision threshold (cutoff point) that is
    used for prediction. The methods for picking cutoff points make use of
    traditional binary classification evaluation statistics such as the true
    positive and true negative rates or any metrics accepting true labels and
    the output of a scoring function from a scikit-learn estimator.

    Parameters
    ----------
    base_estimator : estimator object
        The classifier, fitted or not fitted, from which we want to optimize
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

    n_threshold : int, default=1000
        The number of decision threshold to use when discretizing the output
        of the classifier `method`.

    cv : int, float, cross-validation generator, iterable or "prefit", \
            default=None
        Determines the cross-validation splitting strategy used in
        `cross_val_predict` to train classifier. Possible inputs for cv are:

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

    random_state : int or RandomState, default=None
        Controls the randomness of the training and testing indices produced
        when `cv` is a single shuffle split (i.e., giving a float number).
        See :term:`Glossary <random_state>`.

    n_jobs : int, default=None
        The number of jobs to run in parallel all `estimators` `fit`.
        `None` means 1 unless in a `joblib.parallel_backend` context. -1 means
        using all processors. See Glossary for more details.

    Attributes
    ----------
    decision_threshold_ : float
        The new decision threshold.

    classes_ : array of shape (n_classes,)
        The class labels.

    Examples
    --------
    First, we will load the breast cancer databases and make it highly
    imbalanced.

    >>> import numpy as np
    >>> from sklearn.datasets import load_breast_cancer
    >>> X, y = load_breast_cancer(return_X_y=True)
    >>> pos_idx = np.flatnonzero(y == 1)[:10].tolist()
    >>> neg_idx = np.flatnonzero(y == 0).tolist()
    >>> X, y = X[pos_idx + neg_idx, :], y[pos_idx + neg_idx]

    Then, we can split into a training and testing set and keep the
    same imbalance level in both sets.

    >>> from sklearn.model_selection import train_test_split
    >>> X_train, X_test, y_train, y_test = train_test_split(
    ...     X, y, stratify=y, random_state=0
    ... )

    We can check the performance of a logistic regression model.

    >>> from sklearn.preprocessing import StandardScaler
    >>> from sklearn.linear_model import LogisticRegression
    >>> from sklearn.pipeline import make_pipeline
    >>> model = make_pipeline(StandardScaler(), LogisticRegression())
    >>> model.fit(X_train, y_train)
    Pipeline(steps=[('standardscaler', StandardScaler()),
                    ('logisticregression', LogisticRegression())])
    >>> from sklearn.metrics import balanced_accuracy_score
    >>> y_pred = model.predict(X_test)
    >>> print(f"Score: {balanced_accuracy_score(y_test, y_pred):.3f}")
    Score: 0.833

    We will try to correct the decision threshold which is impacted by the
    class imbalanced.

    >>> from sklearn.model_selection import CutoffClassifier
    >>> model_optimized = CutoffClassifier(
    ...     base_estimator=model, objective_metric=balanced_accuracy_score
    ... )
    >>> model_optimized.fit(X, y)
    CutoffClassifier(...)
    >>> y_pred = model_optimized.predict(X_test)
    >>> print(f"Score: {balanced_accuracy_score(y_test, y_pred):.3f}")
    Score: 0.972
    """

    def __init__(
        self,
        base_estimator,
        objective_metric="balanced_accuracy",
        objective_value=None,
        response_method="auto",
        n_threshold=1000,
        cv=None,
        refit="auto",
        random_state=None,
        n_jobs=None,
    ):
        self.base_estimator = base_estimator
        self.objective_metric = objective_metric
        self.objective_value = objective_value
        self.response_method = response_method
        self.n_threshold = n_threshold
        self.cv = cv
        self.refit = refit
        self.random_state = random_state
        self.n_jobs = n_jobs

    def _validate_parameters(self):
        """Validate the input parameters."""
        if (
            self.objective_metric not in ("tpr", "tnr")
            and self.objective_value is not None
        ):
            raise ValueError(
                f"When 'objective_metric' is a scoring function, "
                f"'objective_value' should be None. Got "
                f"{self.objective_value} instead."
            )

        if (
            not isinstance(self.n_threshold, numbers.Integral)
            or self.n_threshold < 0
        ):
            raise ValueError(
                f"'n_threshold' should be a strictly positive integer. "
                f"Got {self.n_threshold} instead."
            )

    @staticmethod
    def _validate_data(X, y):
        y = check_array(y, ensure_2d=False, dtype=None)
        check_classification_targets(y)
        y_type = type_of_target(y)
        if y_type != 'binary':
            raise ValueError(f'Expected target of binary type. Got {y_type}.')
        return X, y

    def _check_cv_refit(self, cv, refit, y, random_state):
        if isinstance(cv, numbers.Real) and 0 < cv < 1:
            cv = StratifiedShuffleSplit(
                n_splits=1, test_size=cv, random_state=random_state
            )
            refit = False if refit == "auto" else refit
        elif cv == "prefit":
            if refit is True:
                raise ValueError("When cv='prefit', refit cannot be True.")
            refit = False
        else:
            cv = check_cv(cv, y=y, classifier=True)
            if refit is False:
                raise ValueError(
                    "When cv has several folds, refit cannot be False"
                )
            refit = True
        return cv, refit

    @staticmethod
    def _fit_and_score(
        estimator, X, y, train_idx, val_idx, scorer, score_method
    ):
        if train_idx is not None:
            X_train = _safe_indexing(X, train_idx)
            X_val = _safe_indexing(X, val_idx)
            y_train = _safe_indexing(y, train_idx)
            y_val = _safe_indexing(y, val_idx)

            estimator.fit(X_train, y_train)
        else:
            X_val, y_val = X, y

        if score_method in ("tnr", "tpr"):
            fpr, tpr, potential_thresholds = scorer(estimator, X_val, y_val)
            score_thresholds = tpr
            if score_method == "tnr":
                score_thresholds = (1 - fpr)[::-1]
                potential_thresholds = potential_thresholds[::-1]
        else:
            potential_thresholds, score_thresholds = scorer(
                estimator, X_val, y_val
            )

        return potential_thresholds, score_thresholds

    @staticmethod
    def _find_decision_threshold(
        thresholds, scores, n_thresholds, objective_score
    ):
        min_threshold = np.min([th.min() for th in thresholds])
        max_threshold = np.max([th.max() for th in thresholds])
        ascending = thresholds[0].argmin() == 0
        start = min_threshold if ascending else max_threshold
        stop = max_threshold if ascending else min_threshold
        thresholds_interpolated = np.linspace(start, stop, num=n_thresholds)
        mean_score = np.mean(
            [np.interp(thresholds_interpolated,
                       thresholds[fold_idx], scores[fold_idx])
             for fold_idx in range(len(scores))],
            axis=0
        )
        if objective_score == "highest":
            threshold_idx = mean_score.argmax()
        else:
            threshold_idx = np.searchsorted(mean_score, objective_score)
        return thresholds_interpolated[threshold_idx]

    def fit(self, X, y):
        """Find the decision threshold.

        Parameters
        ----------
        X : {array-like, sparse matrix, dataframe} of shape \
                (n_samples, n_features)
            The training data.

        y : array-like of shape (n_samples,)
            Target values. It should be a binary target.

        Returns
        -------
        self : object
            Returns an instance of self.
        """
        X, y = self._validate_data(X, y)

        cv, refit = self._check_cv_refit(
            self.cv, self.refit, y, self.random_state
        )

        # Start by fitting the final estimator
        if refit:
            self._estimator = clone(self.base_estimator).fit(X, y)
        elif cv == "prefit":
            check_is_fitted(self.base_estimator, attributes=["classes_"])
            self._estimator = self.base_estimator
        else:  # single shuffle split CV
            train_idx, _ = next(cv.split(X, y))
            X_train = _safe_indexing(X, train_idx)
            y_train = _safe_indexing(y, train_idx)
            self._estimator = clone(self.base_estimator).fit(X_train, y_train)

        self.classes_ = self._estimator.classes_
        if len(self.classes_) == 1:
            raise ValueError(
                f"This classifier needs samples from 2 classes in the data "
                f"to be trained but the data contains only the class: "
                f"{self.classes_.item(0)}"
            )

        # delayed the parameters check until we have a fitted base estimator
        # with known classes
        self._validate_parameters()

        if cv == "prefit" or not refit:
            model = self._estimator
            splits = ([None, range(len(X))],)
        else:
            model = clone(self.base_estimator)
            splits = cv.split(X, y)

        if self.objective_metric in ("tpr", "tnr"):
            scoring = make_scorer(roc_curve, needs_threshold=True)
        else:
            scoring = check_scoring(
                estimator=model, scoring=self.objective_metric
            )
            if isinstance(scoring, _BaseScorer):
                scoring = _ContinuousScorer(
                    score_func=scoring._score_func,
                    sign=scoring._sign,
                    response_method=self.response_method,
                    kwargs=scoring._kwargs,
                )
        self._scorer = check_scoring(estimator=model, scoring=scoring)

        thresholds, scores = zip(
            *Parallel(n_jobs=self.n_jobs)(
                delayed(self._fit_and_score)(
                    model,
                    X,
                    y,
                    train_idx,
                    val_idx,
                    self._scorer,
                    self.objective_metric,
                )
                for train_idx, val_idx in splits
            )
        )

        if self.objective_metric in ("tnr", "tpr"):
            objective_value = self.objective_value
        else:
            objective_value = "highest"
        self.decision_threshold_ = self._find_decision_threshold(
            thresholds, scores, self.n_threshold, objective_value
        )

        return self

    def predict(self, X):
        """Predict using the calibrated decision threshold

        Parameters
        ----------
        X : {array-like, sparse matrix, dataframe} of shape \
                (n_samples, n_features)
            The data matrix.

        Returns
        -------
        C : ndarray of shape (n_samples,)
            The predicted class.
        """
        check_is_fitted(self)

        response_method = _check_classifier_response_method(
            estimator=self._estimator, response_method=self.response_method
        )

        y_score = response_method(X)
        if response_method.__name__ == "decision_function":
            y_score = self._scorer._check_decision_function(
                y_score, self.classes_
            )
        else:
            y_score = self._scorer._select_proba(
                y_score, self.classes_, support_multi_class=False
            )
        y_class_indices = (y_score >= self.decision_threshold_).astype(int)

        return self.classes_[y_class_indices]

    def _more_tags(self):
        return {
            "binary_only": True,
            "_xfail_test": {
                "check_classifiers_classes":
                "requires non default 'pos_label='two'' parameter",
                "check_fit2d_1feature":
                "requires non default 'pos_label=2' parameter",
            }
        }
