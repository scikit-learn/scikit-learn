"""Input validation and categorical encoding shared by tree-based estimators.

Used by decision trees, forests and gradient boosting.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
from scipy.sparse import issparse

from sklearn.preprocessing import OrdinalEncoder
from sklearn.utils import _safe_indexing
from sklearn.utils._missing import is_scalar_nan
from sklearn.utils.validation import (
    _check_categorical_features,
    _check_n_features,
    _num_samples,
    check_array,
    validate_data,
)


def _encode_categorical_features(estimator, X, *, reset, dtype):
    """Ordinal-encode the categorical features of X, keeping the column order.

    At fit time (``reset=True``), also detect the categorical features and fit
    the encoder, stored in ``estimator._categorical_encoder`` (None if there is
    no categorical feature). X is returned unchanged if there is no categorical
    feature. Otherwise, a Fortran-ordered array of dtype ``dtype`` is returned.
    """
    if reset:
        estimator._categorical_encoder = None
        # RandomTreesEmbedding and gradient boosting have no
        # `categorical_features` parameter.
        if hasattr(estimator, "categorical_features"):
            estimator.is_categorical_ = _check_categorical_features(
                X, estimator.categorical_features
            )
        if getattr(estimator, "is_categorical_", None) is not None:
            estimator._categorical_encoder = OrdinalEncoder(
                dtype=dtype,
                categories="auto",
                handle_unknown="use_encoded_value",
                unknown_value=np.nan,
                encoded_missing_value=np.nan,
            )

    encoder = estimator._categorical_encoder
    if encoder is None:
        return X

    if issparse(X):
        raise NotImplementedError(
            "Categorical features not supported with sparse inputs"
        )

    if not (hasattr(X, "__array__") or hasattr(X, "__dataframe__")):
        # Lists of lists do not support column indexing.
        X = check_array(X, dtype=object, ensure_all_finite=False)
    if not reset:
        _check_n_features(estimator, X, reset=False)

    is_categorical = estimator.is_categorical_
    X_out = np.empty((_num_samples(X), is_categorical.shape[0]), dtype, order="F")

    X_cat = _safe_indexing(X, is_categorical, axis=1)
    X_out[:, is_categorical] = (
        encoder.fit_transform(X_cat) if reset else encoder.transform(X_cat)
    )
    if not is_categorical.all():
        X_out[:, ~is_categorical] = check_array(
            _safe_indexing(X, ~is_categorical, axis=1),
            input_name="X",
            estimator=estimator,
            dtype=dtype,
            ensure_all_finite=False,
        )
    return X_out


def _validate_X(
    estimator,
    X,
    y="no_validation",
    *,
    reset,
    dtype,
    accept_sparse=False,
    ensure_all_finite=True,
    check_input=True,
):
    """Validate X and ordinal-encode its categorical features.

    Used at fit time (``reset=True``) and predict time (``reset=False``). With
    ``check_input=False``, X is assumed to be already validated: only categorical
    features are encoded and the number of features is checked.

    ``y`` is only passed to :func:`validate_data` to raise an informative error
    if it is None at fit time. It is neither validated nor returned.
    """
    if check_input:
        # Check feature names on the original input before categorical
        # encoding converts it to a NumPy array and drops dataframe metadata.
        # The number of features is checked after check_array
        # (ensure_2d=False here), so that 1D inputs get check_array's
        # informative error.
        validate_data(
            estimator, X, y, reset=reset, skip_check_array=True, ensure_2d=False
        )

    X = _encode_categorical_features(estimator, X, reset=reset, dtype=dtype)

    if check_input:
        X = check_array(
            X,
            input_name="X",
            estimator=estimator,
            dtype=dtype,
            accept_sparse=accept_sparse,
            ensure_all_finite=ensure_all_finite,
        )
    _check_n_features(estimator, X, reset=reset)
    return X


def _get_n_categories(estimator):
    """Number of categories of each feature, -1 for numerical features.

    Missing values are not counted as a category.
    """
    n_categories = np.full(estimator.n_features_in_, -1, dtype=np.intp)
    encoder = estimator._categorical_encoder
    if encoder is not None:
        for idx, categories in zip(
            np.flatnonzero(estimator.is_categorical_), encoder.categories_
        ):
            # OrdinalEncoder places np.nan last if missing values reach fit.
            has_nan = len(categories) > 0 and is_scalar_nan(categories[-1])
            n_categories[idx] = len(categories) - has_nan
    return n_categories
