"""Utilities to get the response values of a classifier or a regressor.

It allows to make uniform checks and validation.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from sklearn.base import is_classifier, is_outlier_detector, is_regressor
from sklearn.utils.multiclass import type_of_target
from sklearn.utils.validation import _check_response_method, check_is_fitted


def _process_predict_proba(*, y_pred, target_type, classes, pos_label):
    """Get the response values when the response method is `predict_proba`.

    This function process the `y_pred` array in the binary and multi-label cases.
    In the binary case, it selects the column corresponding to the positive
    class. In the multi-label case, it stacks the predictions if they are not
    in the "compressed" format `(n_samples, n_outputs)`.

    Parameters
    ----------
    y_pred : ndarray
        Output of `estimator.predict_proba`. The shape depends on the target type:

        - for binary classification, it is a 2d array of shape `(n_samples, 2)`;
        - for multiclass classification, it is a 2d array of shape
          `(n_samples, n_classes)`;
        - for multilabel classification, it is either a list of 2d arrays of shape
          `(n_samples, 2)` (e.g. `RandomForestClassifier` or `KNeighborsClassifier`) or
          an array of shape `(n_samples, n_outputs)` (e.g. `MLPClassifier` or
          `RidgeClassifier`).

    target_type : {"binary", "multiclass", "multilabel-indicator"}
        Type of the target.

    classes : ndarray of shape (n_classes,) or list of such arrays
        Class labels as reported by `estimator.classes_`.

    pos_label : int, float, bool or str
        Only used with binary and multiclass targets.

    Returns
    -------
    y_pred : ndarray of shape (n_samples,), (n_samples, n_classes) or \
            (n_samples, n_output)
        Compressed predictions format as requested by the metrics.
    """
    if target_type == "binary" and y_pred.shape[1] < 2:
        # We don't handle classifiers trained on a single class.
        raise ValueError(
            f"Got predict_proba of shape {y_pred.shape}, but need "
            "classifier with two classes."
        )

    if target_type == "binary":
        col_idx = np.flatnonzero(classes == pos_label)[0]
        return y_pred[:, col_idx]
    elif target_type == "multilabel-indicator":
        # Use a compress format of shape `(n_samples, n_output)`.
        if isinstance(y_pred, list):
            # list of arrays of shape `(n_samples, 2)`
            return np.vstack([p[:, -1] for p in y_pred]).T
        else:
            # array of shape `(n_samples, n_outputs)`
            return y_pred

    return y_pred


def _process_decision_function(*, y_pred, target_type, classes, pos_label):
    """Get the response values when the response method is `decision_function`."""

    if target_type == "binary" and pos_label == classes[0]:
        return -1 * y_pred
    return y_pred


def _get_response_values(
    estimator,
    X,
    response_method,
    pos_label=None,
    return_response_method_used=False,
):
    """Compute the response values of a classifier, an outlier detector, or a regressor."""

    if is_classifier(estimator):
        prediction_method = _check_response_method(estimator, response_method)
        classes = estimator.classes_
        target_type = type_of_target(classes)

        if target_type in ("binary", "multiclass"):
            if pos_label is not None and pos_label not in classes.tolist():
                raise ValueError(
                    f"pos_label={pos_label} is not a valid label: It should be "
                    f"one of {classes}"
                )
            elif pos_label is None and target_type == "binary":
                pos_label = classes[-1]

        y_pred = prediction_method(X)

        if prediction_method.__name__ in ("predict_proba", "predict_log_proba"):
            y_pred = _process_predict_proba(
                y_pred=y_pred,
                target_type=target_type,
                classes=classes,
                pos_label=pos_label,
            )
        elif prediction_method.__name__ == "decision_function":
            y_pred = _process_decision_function(
                y_pred=y_pred,
                target_type=target_type,
                classes=classes,
                pos_label=pos_label,
            )

    elif is_outlier_detector(estimator):
        prediction_method = _check_response_method(estimator, response_method)
        y_pred, pos_label = prediction_method(X), None

    elif is_regressor(estimator):
        if response_method != "predict":
            raise ValueError(
                f"{estimator.__class__.__name__} should either be a classifier to be "
                f"used with response_method={response_method} or the response_method "
                "should be 'predict'. Got a regressor with response_method="
                f"{response_method} instead."
            )
        prediction_method = estimator.predict
        y_pred, pos_label = prediction_method(X), None

    else:
        raise ValueError(
            f"{estimator.__class__.__name__} is not recognized as a classifier, "
            "an outlier detector, or a regressor."
        )

    if return_response_method_used:
        return y_pred, pos_label, prediction_method.__name__
    return y_pred, pos_label


def _get_response_values_binary(
    estimator, X, response_method, pos_label=None, return_response_method_used=False
):
    """Compute the response values of a binary classifier."""

    classification_error = "Expected 'estimator' to be a binary classifier."

    check_is_fitted(estimator)
    if not is_classifier(estimator):
        raise ValueError(
            classification_error + f" Got {estimator.__class__.__name__} instead."
        )
    elif len(estimator.classes_) != 2:
        raise ValueError(
            classification_error + f" Got {len(estimator.classes_)} classes instead."
        )

    if response_method == "auto":
        response_method = ["predict_proba", "decision_function"]

    return _get_response_values(
        estimator,
        X,
        response_method,
        pos_label=pos_label,
        return_response_method_used=return_response_method_used,
    )
