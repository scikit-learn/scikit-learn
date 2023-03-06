import numpy as np

from .validation import _check_response_method


def _get_response_values(
    estimator,
    X,
    response_method,
    pos_label=None,
):
    """Compute the response values of a classifier or a regressor.

    The response values are predictions, one scalar value for each sample in X
    that depends on the specific choice of `response_method`.

    This helper only accepts multiclass classifiers with the `predict` response
    method.

    If `estimator` is a binary classifier, also return the label for the
    effective positive class.

    .. versionadded:: 1.3

    Parameters
    ----------
    estimator : estimator instance
        Fitted classifier or regressor or a fitted :class:`~sklearn.pipeline.Pipeline`
        in which the last estimator is a classifier or a regressor.

    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        Input values.

    response_method : {"predict_proba", "decision_function", "predict"} or \
            list of such str
        Specifies the response method to use get prediction from an estimator
        (i.e. :term:`predict_proba`, :term:`decision_function` or
        :term:`predict`). Possible choices are:

        - if `str`, it corresponds to the name to the method to return;
        - if a list of `str`, it provides the method names in order of
          preference. The method returned corresponds to the first method in
          the list and which is implemented by `estimator`.

    pos_label : str or int, default=None
        The class considered as the positive class when computing
        the metrics. By default, `estimators.classes_[1]` is
        considered as the positive class.

    Returns
    -------
    y_pred : ndarray of shape (n_samples,)
        Target scores calculated from the provided response_method
        and `pos_label`.

    pos_label : str, int or None
        The class considered as the positive class when computing
        the metrics. Returns `None` if `estimator` is a regressor.

    Raises
    ------
    ValueError
        If `pos_label` is not a valid label.
        If the shape of `y_pred` is not consistent for binary classifier.
        If the response method can be applied to a classifier only and
        `estimator` is a regressor.
    """
    from sklearn.base import is_classifier  # noqa

    if is_classifier(estimator):
        prediction_method = _check_response_method(estimator, response_method)
        y_pred = prediction_method(X)
        classes = estimator.classes_

        target_type = "binary" if len(classes) <= 2 else "multiclass"

        if target_type == "multiclass" and prediction_method.__name__ != "predict":
            raise ValueError(
                "With a multiclass estimator, the response method should be "
                f"predict, got {prediction_method.__name__} instead."
            )

        if pos_label is not None and pos_label not in classes.tolist():
            raise ValueError(
                f"pos_label={pos_label} is not a valid label: It should be "
                f"one of {classes}"
            )
        elif pos_label is None and target_type == "binary":
            pos_label = pos_label if pos_label is not None else classes[-1]

        if prediction_method.__name__ == "predict_proba":
            if target_type == "binary" and y_pred.shape[1] <= 2:
                if y_pred.shape[1] == 2:
                    col_idx = np.flatnonzero(classes == pos_label)[0]
                    y_pred = y_pred[:, col_idx]
                else:
                    err_msg = (
                        f"Got predict_proba of shape {y_pred.shape}, but need "
                        "classifier with two classes."
                    )
                    raise ValueError(err_msg)
        elif prediction_method.__name__ == "decision_function":
            if target_type == "binary":
                if pos_label == classes[0]:
                    y_pred *= -1
    else:
        if response_method != "predict":
            raise ValueError(f"{estimator.__class__.__name__} should be a classifier")
        y_pred, pos_label = estimator.predict(X), None

    return y_pred, pos_label
