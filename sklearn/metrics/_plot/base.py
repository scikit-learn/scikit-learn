import numpy as np

from sklearn.base import is_classifier


def _check_classifier_response_method(estimator, response_method):
    """Return prediction method from the response_method

    Parameters
    ----------
    estimator: object
        Classifier to check

    response_method: {'auto', 'predict_proba', 'decision_function'}
        Specifies whether to use :term:`predict_proba` or
        :term:`decision_function` as the target response. If set to 'auto',
        :term:`predict_proba` is tried first and if it does not exist
        :term:`decision_function` is tried next.

    Returns
    -------
    prediction_method: callable
        prediction method of estimator
    """

    if response_method not in ("predict_proba", "decision_function", "auto"):
        raise ValueError("response_method must be 'predict_proba', "
                         "'decision_function' or 'auto'")

    error_msg = "response method {} is not defined in {}"
    if response_method != "auto":
        prediction_method = getattr(estimator, response_method, None)
        if prediction_method is None:
            raise ValueError(error_msg.format(response_method,
                                              estimator.__class__.__name__))
    else:
        predict_proba = getattr(estimator, 'predict_proba', None)
        decision_function = getattr(estimator, 'decision_function', None)
        prediction_method = predict_proba or decision_function
        if prediction_method is None:
            raise ValueError(error_msg.format(
                "decision_function or predict_proba",
                estimator.__class__.__name__))

    return prediction_method


def _get_response(X, estimator, response_method, pos_label=None):
    """Return response and positive label.

    Parameters
    ----------
    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        Input values.

    estimator : estimator instance
        Fitted classifier or a fitted :class:`~sklearn.pipeline.Pipeline`
        in which the last estimator is a classifier.

    response_method: {'auto', 'predict_proba', 'decision_function'}
        Specifies whether to use :term:`predict_proba` or
        :term:`decision_function` as the target response. If set to 'auto',
        :term:`predict_proba` is tried first and if it does not exist
        :term:`decision_function` is tried next.

    pos_label : str or int, default=None
        The class considered as the positive class when computing
        the metrics. By default, `estimators.classes_[1]` is
        considered as the positive class.

    Returns
    -------
    y_pred: ndarray of shape (n_samples,)
        Target scores calculated from the provided response_method
        and pos_label.

    pos_label: str or int
        The class considered as the positive class when computing
        the metrics.
    """
    classification_error = (
        "{} should be a binary classifier".format(estimator.__class__.__name__)
    )

    if not is_classifier(estimator):
        raise ValueError(classification_error)

    prediction_method = _check_classifier_response_method(
        estimator, response_method)

    y_pred = prediction_method(X)

    if pos_label is not None and pos_label not in estimator.classes_:
        raise ValueError(
            f"The class provided by 'pos_label' is unknown. Got "
            f"{pos_label} instead of one of {estimator.classes_}"
        )

    if y_pred.ndim != 1:  # `predict_proba`
        if y_pred.shape[1] != 2:
            raise ValueError(classification_error)
        if pos_label is None:
            pos_label = estimator.classes_[1]
            y_pred = y_pred[:, 1]
        else:
            class_idx = np.flatnonzero(estimator.classes_ == pos_label)
            y_pred = y_pred[:, class_idx]
    else:
        if pos_label is None:
            pos_label = estimator.classes_[1]
        elif pos_label == estimator.classes_[0]:
            y_pred *= -1

    return y_pred, pos_label
