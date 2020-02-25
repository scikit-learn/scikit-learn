from inspect import signature

import numpy as np

from ..base import clone
from ..base import BaseEstimator
from ..base import ClassifierMixin
from ..base import MetaEstimatorMixin
from ..exceptions import NotFittedError
from ..metrics import balanced_accuracy_score
from ..metrics import roc_curve
from ..preprocessing import LabelEncoder
from ..utils import check_array
from ..utils.multiclass import check_classification_targets
from ..utils.multiclass import type_of_target
from ..utils.validation import check_is_fitted


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

    objective_metric : {"tpr", "tnr"} or callable, \
            default=balanced_accuracy_score
        The objective metric to be optimized. Can be on of:

        * `"tpr"`: Find the decision threshold for a true positive ratio (TPR)
          of `objective_value`.
        * `"tnr"`: Find the decision threshold for a true negative ratio (TNR)
          of `objective_value`.
        * a callable with the signature `metric(y_true, y_score, **kwargs)`.

    objective_metric_params : dict, default=None
        Some extra parameters to pass to `objective_metric`.

    objective_value : float, default=None
        The value associated with the `objective_metric` metric for which we
        want to find the decision threshold. Only apply when `objective_metric`
        is `"tpr"` or `"tnr"`

    method : {"auto", "decision_function", "predict_proba"}, default="auto"
        Methods by the classifier `base_estimator` corresponding to the
        decision function for which we want to find a threshold. It can be:

        * if `"auto"`, it will try to invoke, for each estimator,
          `"decision_function` or `"predict_proba"` in that order.
        * otherwise, one of `"predict_proba"` or `"decision_function"`.
          If the method is not implemented by the estimator, it will raise an
          error.

    pos_label : int or str, default=None
        The label of the positive class. When `pos_label=None`, if `y_true` is
        in `{-1, 1}` or `{0, 1}`, `pos_label` is set to 1, otherwise an error
        will be raised.

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
    Score: 0.962
    """

    def __init__(
        self,
        base_estimator,
        objective_metric=balanced_accuracy_score,
        objective_metric_params=None,
        objective_value=None,
        method="auto",
        pos_label=None
    ):
        self.base_estimator = base_estimator
        self.objective_metric = objective_metric
        self.objective_metric_params = objective_metric_params
        self.objective_value = objective_value
        self.method = method
        self.pos_label = pos_label

    def _validate_parameters(self):
        """Validate the input parameters."""
        supported_methods = ("decision_function", "predict_proba")
        if self.method == "auto":
            has_methods = [
                hasattr(self.base_estimator, m) for m in supported_methods
            ]
            if not any(has_methods):
                raise TypeError(
                    f"'base_estimator' must implement one of the "
                    f"{', '.join(supported_methods)} methods."
                )
            self._method = next(
                (m for m, i in zip(supported_methods, has_methods) if i), None
            )
        else:
            if self.method not in supported_methods:
                raise ValueError(
                    f"'method' should be one of {', '.join(supported_methods)}"
                    f". Got {self.method} instead."
                )
            elif not hasattr(self.base_estimator, self.method):
                raise TypeError(
                    f"'base_estimator' does not implement {self.method}."
                )
            self._method = self.method
        if (self.objective_metric not in ("tnr", "tpr") and
                self.objective_value is not None):
            raise ValueError(
                f"When 'objective_metric' is a scoring function, "
                f"'objective_value' should be None. Got {self.objective_value}"
                f" instead."
            )

        # ensure binary classification if `pos_label` is not specified
        # `classes.dtype.kind` in ('O', 'U', 'S') is required to avoid
        # triggering a FutureWarning by calling np.array_equal(a, b)
        # when elements in the two arrays are not comparable.
        if (self.pos_label is None and (
                self.classes_.dtype.kind in ('O', 'U', 'S') or
                not (np.array_equal(self.classes_, [0, 1]) or
                     np.array_equal(self.classes_, [-1, 1]) or
                     np.array_equal(self.classes_, [0]) or
                     np.array_equal(self.classes_, [-1]) or
                     np.array_equal(self.classes_, [1])))):
            classes_repr = ", ".join(repr(c) for c in self.classes_)
            raise ValueError(
                f"'y_true' takes value in {classes_repr} and 'pos_label' is "
                f"not specified: either make 'y_true' take value in "
                "{{0, 1}} or {{-1, 1}} or pass pos_label explicitly."
            )
        elif self.pos_label is None:
            self._pos_label = 1
        else:
            self._pos_label = self.pos_label

    @staticmethod
    def _validate_data(X, y):
        y = check_array(y, ensure_2d=False, dtype=None)
        check_classification_targets(y)
        y_type = type_of_target(y)
        if y_type != 'binary':
            raise ValueError(f'Expected target of binary type. Got {y_type}.')
        return X, y

    @staticmethod
    def _get_pos_label_score(y_score, classes, pos_label):
        """Get score of the positive class."""
        if y_score.ndim == 2:
            pos_label_encoded = np.flatnonzero(classes == pos_label).item(0)
            y_score = y_score[:, pos_label_encoded]
        return y_score

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

        try:
            check_is_fitted(self.base_estimator, attributes=["classes_"])
            self._estimator = self.base_estimator
        except NotFittedError:
            self._estimator = clone(self.base_estimator).fit(X, y)
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

        # warm start a label encoder using the fitted estimator
        label_encoder = LabelEncoder()
        label_encoder.classes_ = self.classes_

        y_encoded = label_encoder.transform(y)
        self._pos_label_encoded = np.flatnonzero(
            self.classes_ == self._pos_label
        ).item(0)

        y_score = getattr(self._estimator, self._method)(X)
        if self.objective_metric in ("tpr", "tnr"):
            fpr, tpr, thresholds = roc_curve(
                y_encoded, y_score, pos_label=self._pos_label_encoded
            )
            metric = tpr
            if self.objective_metric == "tnr":
                tnr, thresholds = (1 - fpr)[::-1], thresholds[::-1]
                metric = tnr

            threshold_idx = np.searchsorted(
                metric, self.objective_value
            )
            self.decision_threshold_ = thresholds[threshold_idx]
        else:
            # `np.unique` is already sorting the value, no need to call
            #  `thresholds.sort()`
            thresholds = np.unique(
                self._get_pos_label_score(
                    y_score, self.classes_, self.pos_label
                )
            )
            params = ({} if self.objective_metric_params is None
                      else self.objective_metric_params)
            metric_signature = signature(self.objective_metric)
            if "pos_label" in metric_signature.parameters:
                params["pos_label"] = self._pos_label_encoded
            scores = [
                self.objective_metric(
                    y_encoded, (y_score >= th).astype(int), **params
                )
                for th in thresholds
            ]
            self.decision_threshold_ = thresholds[np.argmax(scores)]

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

        decision_function = getattr(self._estimator, self._method)
        y_score = self._get_pos_label_score(
            decision_function(X), self.classes_, self._pos_label
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
