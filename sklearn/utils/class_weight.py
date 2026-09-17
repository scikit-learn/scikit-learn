"""Utilities for handling weights based on class labels."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
from scipy import sparse

from sklearn.utils._array_api import (
    _bincount,
    _is_numpy_namespace,
    get_namespace_and_device,
    move_to,
    size,
    xpx,
)
from sklearn.utils._param_validation import StrOptions, validate_params
from sklearn.utils.validation import _check_sample_weight


@validate_params(
    {
        "class_weight": [dict, StrOptions({"balanced"}), None],
        "classes": [np.ndarray],
        "y": ["array-like"],
        "sample_weight": ["array-like", None],
    },
    prefer_skip_nested_validation=True,
)
def compute_class_weight(class_weight, *, classes, y, sample_weight=None):
    """Estimate class weights for unbalanced datasets.

    Parameters
    ----------
    class_weight : dict, "balanced" or None
        If "balanced", class weights will be given by
        `n_samples / (n_classes * np.bincount(y))` or their weighted equivalent if
        `sample_weight` is provided.
        If a dictionary is given, keys are classes and values are corresponding class
        weights.
        If `None` is given, the class weights will be uniform.

    classes : ndarray
        Array of the classes occurring in the data, as given by
        `np.unique(y_org)` with `y_org` the original class labels.

    y : array-like of shape (n_samples,)
        Array of original class labels per sample.

    sample_weight : array-like of shape (n_samples,), default=None
        Array of weights that are assigned to individual samples. Only used when
        `class_weight='balanced'`.

    Returns
    -------
    class_weight_vect : ndarray of shape (n_classes,)
        Array with `class_weight_vect[i]` the weight for i-th class.

    References
    ----------
    The "balanced" heuristic is inspired by
    Logistic Regression in Rare Events Data, King, Zen, 2001.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.utils.class_weight import compute_class_weight
    >>> y = [1, 1, 1, 1, 0, 0]
    >>> compute_class_weight(class_weight="balanced", classes=np.unique(y), y=y)
    array([1.5 , 0.75])
    """
    # Import error caused by circular imports.
    from sklearn.preprocessing import LabelEncoder

    xp, _, device = get_namespace_and_device(y, classes)
    unique_y = xp.unique_values(y)
    if set(move_to(unique_y, xp=np, device="cpu")) - set(
        move_to(classes, xp=np, device="cpu")
    ):
        raise ValueError("classes should include all valid labels that are in y")

    if class_weight is None or len(class_weight) == 0:
        # uniform class weights
        weight = xp.ones(size(classes), device=device)
    elif class_weight == "balanced":
        # Find the weight of each class as present in y.
        le = LabelEncoder()
        y_ind = le.fit_transform(y)
        if not all(xpx.isin(classes, xp.astype(le.classes_, classes.dtype), xp=xp)):
            raise ValueError("classes should have valid labels that are in y")

        if _is_numpy_namespace(xp) and sample_weight is not None:
            sample_weight = move_to(sample_weight, xp=np, device="cpu")

        sample_weight = _check_sample_weight(sample_weight, y)
        weighted_class_counts = _bincount(y_ind, weights=sample_weight, xp=xp)
        recip_freq = xp.sum(weighted_class_counts) / (
            size(le.classes_) * weighted_class_counts
        )
        weight = recip_freq[le.transform(classes)]
    else:
        # user-defined dictionary
        weight = xp.ones(size(classes), device=device)
        unweighted_classes = []
        for i, c in enumerate(classes):
            # Use the class label directly for dict lookup.
            # Previously, int(c) was attempted which coerced string labels
            # like "1" to integers, breaking lookup when class_weight has
            # string keys.
            if c in class_weight:
                weight[i] = class_weight[c]
            else:
                unweighted_classes.append(c)

        n_weighted_classes = size(classes) - len(unweighted_classes)
        if unweighted_classes and n_weighted_classes != len(class_weight):
            unweighted_classes_user_friendly_str = np.array(unweighted_classes).tolist()
            raise ValueError(
                f"The classes, {unweighted_classes_user_friendly_str}, are not in"
                " class_weight"
            )

    return weight


@validate_params(
    {
        "class_weight": [
            dict,
            list,
            StrOptions({"balanced"}),
            None,
        ],
        "y": ["array-like", "sparse matrix"],
        "indices": ["array-like", None],
    },
    prefer_skip_nested_validation=True,
)
def compute_sample_weight(class_weight, y, *, indices=None):
    """Estimate sample weights by class for unbalanced datasets.

    Parameters
    ----------
    class_weight : dict, list of dicts, "balanced", or None
        Weights associated with classes in the form `{class_label: weight}`.
        If not given, all classes are supposed to have weight one. For
        multi-output problems, a list of dicts can be provided in the same
        order as the columns of y.

        Note that for multioutput (including multilabel) weights should be
        defined for each class of every column in its own dict. For example,
        for four-class multilabel classification weights should be
        `[{0: 1, 1: 1}, {0: 1, 1: 5}, {0: 1, 1: 1}, {0: 1, 1: 1}]` instead of
        `[{1:1}, {2:5}, {3:1}, {4:1}]`.

        The `"balanced"` mode uses the values of y to automatically adjust
        weights inversely proportional to class frequencies in the input data:
        `n_samples / (n_classes * np.bincount(y))`.

        For multi-output, the weights of each column of y will be multiplied.

    y : {array-like, sparse matrix} of shape (n_samples,) or (n_samples, n_outputs)
        Array of original class labels per sample.

    indices : array-like of shape (n_subsample,), default=None
        Array of indices to be used in a subsample. Can be of length less
        than n_samples in the case of a subsample. If None, the full
        array is used.

    Returns
    -------
    sample_weight_vect : ndarray of shape (n_subsample,)
        Array with `sample_weight_vect[i]` the weight for i-th sample.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.utils.class_weight import compute_sample_weight
    >>> y = [1, 1, 1, 1, 0, 0]
    >>> compute_sample_weight(class_weight="balanced", y=y)
    array([0.75, 0.75, 0.75, 0.75, 1.5 , 1.5 ])
    """
    # Ensure that y is a list of arrays for multi-output
    if sparse.issparse(y):
        if y.format != "csr":
            y = y.tocsr()
        y = [np.array(y[:, i].todense()).ravel() for i in range(y.shape[1])]
        y = np.vstack(y).T
    elif isinstance(y, list) and len(y) > 0 and not isinstance(y[0], (list, tuple, np.ndarray)):
        # y is a list of scalars, treat as single output
        y = np.asarray(y)
    elif not hasattr(y, "shape"):
        y = np.asarray(y)

    if y.ndim == 1:
        y = np.reshape(y, (-1, 1))

    n_outputs = y.shape[1]

    if indices is not None:
        y = y[indices]

    if isinstance(class_weight, list):
        if len(class_weight) != n_outputs:
            raise ValueError(
                "For multi-output, number of class weights should match number of"
                f" outputs. Got {len(class_weight)} class weights and {n_outputs}"
                " outputs."
            )
    else:
        class_weight = [class_weight] * n_outputs

    weights = []
    for k in range(n_outputs):
        weight_k = compute_class_weight(
            class_weight[k], classes=np.unique(y[:, k]), y=y[:, k]
        )
        weight_k = weight_k[np.searchsorted(np.unique(y[:, k]), y[:, k])]
        weights.append(weight_k)

    return np.prod(weights, axis=0)
