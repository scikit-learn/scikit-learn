import numpy as np

from ..base import is_classifier, is_regressor
from ..exceptions import NotFittedError


def _check_feature_names(X, feature_names=None):
    """Check feature names.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Input data.

    feature_names : None or array-like of shape (n_names,), dtype=str
        Feature names to check or `None`.

    Returns
    -------
    feature_names : list of str
        Feature names validated. If `feature_names` is `None`, then a list of
        feature names is provided, i.e. the column names of a pandas dataframe
        or a generic list of feature names (e.g. `["x0", "x1", ...]`) for a
        NumPy array.
    """
    if feature_names is None:
        if hasattr(X, "columns") and hasattr(X.columns, "tolist"):
            # get the column names for a pandas dataframe
            feature_names = X.columns.tolist()
        else:
            # define a list of numbered indices for a numpy array
            feature_names = [f"x{i}" for i in range(X.shape[1])]
    elif hasattr(feature_names, "tolist"):
        # convert numpy array or pandas index to a list
        feature_names = feature_names.tolist()
    if len(set(feature_names)) != len(feature_names):
        raise ValueError("feature_names should not contain duplicates.")

    return feature_names


def _get_feature_index(fx, feature_names=None):
    """Get feature index.

    Parameters
    ----------
    fx : int or str
        Feature index or name.

    feature_names : list of str, default=None
        All feature names from which to search the indices.

    Returns
    -------
    idx : int
        Feature index.
    """
    if isinstance(fx, str):
        if feature_names is None:
            raise ValueError(
                f"Cannot plot partial dependence for feature {fx!r} since "
                "the list of feature names was not provided, neither as "
                "column names of a pandas data-frame nor via the feature_names "
                "parameter."
            )
        try:
            return feature_names.index(fx)
        except ValueError as e:
            raise ValueError(f"Feature {fx!r} not in feature_names") from e
    return fx


def _robust_predict_for_scatter(X, est, response_method):
    """Get predictions from the estimator for one-way PD scatter plots.

    Parameters
    ----------
    X : ndarray of shape (n_samples, n_features)
        The data on which to compute the predictions.

    est : estimator
        The estimator from which the predictions are computed.

    response_method : str
        The method used to compute the response. Must be one of
        {'auto', 'predict_proba', 'decision_function'}.

    """
    predictions = []

    # define the prediction_method (predict, predict_proba, decision_function).
    if is_regressor(est):
        prediction_method = est.predict
    else:
        predict_proba = getattr(est, "predict_proba", None)
        decision_function = getattr(est, "decision_function", None)
        if response_method == "auto":
            # try predict_proba, then decision_function if it doesn't exist
            prediction_method = predict_proba or decision_function
        else:
            prediction_method = (
                predict_proba
                if response_method == "predict_proba"
                else decision_function
            )
        if prediction_method is None:
            if response_method == "auto":
                raise ValueError(
                    "The estimator has no predict_proba and no "
                    "decision_function method."
                )
            elif response_method == "predict_proba":
                raise ValueError("The estimator has no predict_proba method.")
            else:
                raise ValueError("The estimator has no decision_function method.")

    X_eval = X.copy()
    try:
        # Note: predictions is of shape
        # (n_points,) for non-multioutput regressors
        # (n_points, n_tasks) for multioutput regressors
        # (n_points, 1) for the regressors in cross_decomposition (I think)
        # (n_points, 2) for binary classification
        # (n_points, n_classes) for multiclass classification
        pred = prediction_method(X_eval)
        predictions.append(pred)
    except NotFittedError as e:
        raise ValueError("'estimator' parameter must be a fitted estimator") from e

    n_samples = X.shape[0]

    # reshape to (n_targets, n_instances, n_points) where n_targets is:
    # - 1 for non-multioutput regression and binary classification (shape is
    #   already correct in those cases)
    # - n_tasks for multi-output regression
    # - n_classes for multiclass classification.
    predictions = np.array(predictions).T
    if is_regressor(est) and predictions.ndim == 2:
        # non-multioutput regression, shape is (n_instances, n_points,)
        predictions = predictions.reshape(n_samples, -1)
    elif is_classifier(est) and predictions.shape[0] == 2:
        # Binary classification, shape is (2, n_instances, n_points).
        # we output the effect of **positive** class
        predictions = predictions[1]
        predictions = predictions.reshape(n_samples, -1)

    return predictions
