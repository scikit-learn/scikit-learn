import numpy as np

from ..base import _is_pairwise
from ..metrics.pairwise import linear_kernel, pairwise_distances

from ._tags import _safe_tags


def _enforce_estimator_tags_y(estimator, y):
    """Modify `y` to be compatible with the available estimator tags.

    Parameters
    ----------
    estimator : object
        Estimator object to test.

    y : ndarray
        The data to be converted.

    Returns
    -------
    y : ndarray
        The converted data.
    """
    # Estimators with a `requires_positive_y` tag only accept strictly positive
    # data
    if _safe_tags(estimator, key="requires_positive_y"):
        # Create strictly positive y. The minimal increment above 0 is 1, as
        # y could be of integer dtype.
        y += 1 + abs(y.min())
    # Estimators with a `binary_only` tag only accept up to two unique y values
    if _safe_tags(estimator, key="binary_only") and y.size > 0:
        y = np.where(y == y.flat[0], y, y.flat[0] + 1)
    # Estimators in mono_output_task_error raise ValueError if y is of 1-D
    # Convert into a 2-D y for those estimators.
    if _safe_tags(estimator, key="multioutput_only"):
        return np.reshape(y, (-1, 1))
    return y


def _enforce_estimator_tags_x(estimator, X):
    """Modify `X` to be compatible with the available estimator tags.

    Parameters
    ----------
    estimator : object
        Estimator object to test.

    X : ndarray
        The data to be converted.

    Returns
    -------
    X : ndarray
        The converted data.
    """
    # Pairwise estimators only accept
    # X of shape (`n_samples`, `n_samples`)
    if _is_pairwise(estimator):
        X = X.dot(X.T)
    # Estimators with `1darray` in `X_types` tag only accept
    # X of shape (`n_samples`,)
    if "1darray" in _safe_tags(estimator, key="X_types"):
        X = X[:, 0]
    # Estimators with a `requires_positive_X` tag only accept
    # strictly positive data
    if _safe_tags(estimator, key="requires_positive_X"):
        X -= X.min()
    return X


def _is_pairwise_metric(estimator):
    """Returns True if estimator accepts pairwise metric.

    Parameters
    ----------
    estimator : object
        Estimator object to test.

    Returns
    -------
    out : bool
        True if _pairwise is set to True and False otherwise.
    """
    metric = getattr(estimator, "metric", None)

    return bool(metric == "precomputed")


def _pairwise_estimator_convert_X(X, estimator, kernel=linear_kernel):
    """Convert `X` so to be used by a pairwise estimator.

    Parameters
    ----------
    X : ndarray of shape (n_samples, n_features)
        The data to be converted.

    estimator : object
        An estimator to apply on `X`.

    kernel : callable, default=linear_kernel
        If `estimator` requires a kernel, this parameter will transform `X`
        into a kernel matrix.

    Returns
    -------
    X_new : ndarray of shape (n_samples, n_features) or (n_samples, n_samples)
        The converted `X`.
    """
    if _is_pairwise_metric(estimator):
        return pairwise_distances(X, metric="euclidean")
    if _is_pairwise(estimator):
        return kernel(X, X)
    return X
