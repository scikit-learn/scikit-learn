from ...base import is_classifier
from ...utils._response import _get_response_values
from ...utils.validation import check_is_fitted


def _get_response_binary(estimator, X, response_method, pos_label=None):
    """Return response and positive label.

    Parameters
    ----------
    estimator : estimator instance
        Fitted classifier or a fitted :class:`~sklearn.pipeline.Pipeline`
        in which the last estimator is a classifier.

    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        Input values.

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
        "Expected 'estimator' to be a binary classifier, but got"
        f" {estimator.__class__.__name__}"
    )

    check_is_fitted(estimator)
    if not is_classifier(estimator) or len(estimator.classes_) != 2:
        raise ValueError(classification_error)

    if response_method == "auto":
        response_method = ["predict_proba", "decision_function"]

    return _get_response_values(
        estimator,
        X,
        response_method,
        pos_label=pos_label,
    )
