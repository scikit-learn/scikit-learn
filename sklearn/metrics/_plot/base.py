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
