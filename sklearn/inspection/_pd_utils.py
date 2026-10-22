# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np


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


def _get_is_categorical(
    categorical_features, features_indices, feature_names, n_features
):
    """Tell whether each feature in `features_indices` is categorical.

    Parameters
    ----------
    categorical_features : array-like of shape (n_features,) or shape \
            (n_categorical_features,), dtype={bool, int, str} or None
        Boolean mask, integer indices or names of the categorical features.
        `None` means that no feature is categorical.

    features_indices : array-like of int
        Indices of the features of interest.

    feature_names : list of str
        All feature names from which to search the indices.

    n_features : int
        Number of features in `X`.

    Returns
    -------
    is_categorical : list of bool
        Whether each feature in `features_indices` is categorical.
    """
    if categorical_features is None:
        return [False] * len(features_indices)

    categorical_features = np.asarray(categorical_features)
    if categorical_features.size == 0:
        raise ValueError(
            "Passing an empty list (`[]`) to `categorical_features` is not "
            "supported. Use `None` instead to indicate that there are no "
            "categorical features."
        )
    if categorical_features.dtype.kind == "b":
        # categorical features provided as a list of boolean
        if categorical_features.size != n_features:
            raise ValueError(
                "When `categorical_features` is a boolean array-like, "
                "the array should be of shape (n_features,). Got "
                f"{categorical_features.size} elements while `X` contains "
                f"{n_features} features."
            )
        return [categorical_features[idx] for idx in features_indices]
    if categorical_features.dtype.kind in ("i", "O", "U"):
        # categorical features provided as a list of indices or feature names
        categorical_features_idx = [
            _get_feature_index(cat, feature_names=feature_names)
            for cat in categorical_features
        ]
        return [idx in categorical_features_idx for idx in features_indices]
    raise ValueError(
        "Expected `categorical_features` to be an array-like of boolean,"
        f" integer, or string. Got {categorical_features.dtype} instead."
    )
