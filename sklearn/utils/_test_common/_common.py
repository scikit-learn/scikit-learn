"""Utilities used in tests."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np

from sklearn.metrics import pairwise_distances
from sklearn.metrics.pairwise import linear_kernel
from sklearn.utils import get_tags


def _enforce_estimator_tags_X(estimator, X, kernel=linear_kernel):
    # Estimators with `1darray` in `X_types` tag only accept
    # X of shape (`n_samples`,)
    if get_tags(estimator).input_tags.one_d_array:
        X = X[:, 0]
    # Estimators with a `requires_positive_X` tag only accept
    # strictly positive data
    if get_tags(estimator).input_tags.positive_only:
        X = X - X.min()
    if get_tags(estimator).input_tags.categorical:
        dtype = np.float64 if get_tags(estimator).input_tags.allow_nan else np.int32
        X = np.round((X - X.min())).astype(dtype)

    if estimator.__class__.__name__ == "SkewedChi2Sampler":
        # SkewedChi2Sampler requires X > -skewdness in transform
        X = X - X.min()

    # Pairwise estimators only accept
    # X of shape (`n_samples`, `n_samples`)
    if _is_pairwise_metric(estimator):
        X = pairwise_distances(X, metric="euclidean")
    elif get_tags(estimator).input_tags.pairwise:
        X = kernel(X, X)
    return X


def _enforce_estimator_tags_y(estimator, y):
    # Estimators with a `requires_positive_y` tag only accept strictly positive
    # data
    tags = get_tags(estimator)
    if tags.target_tags.positive_only:
        # Create strictly positive y. The minimal increment above 0 is 1, as
        # y could be of integer dtype.
        y += 1 + abs(y.min())
    if (
        tags.classifier_tags is not None
        and not tags.classifier_tags.multi_class
        and y.size > 0
    ):
        y = np.where(y == y.flat[0], y, y.flat[0] + 1)
    # Estimators in mono_output_task_error raise ValueError if y is of 1-D
    # Convert into a 2-D y for those estimators.
    if tags.target_tags.multi_output and not tags.target_tags.single_output:
        return np.reshape(y, (-1, 1))
    return y


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


def _is_public_parameter(attr):
    return not (attr.startswith("_") or attr.endswith("_"))
