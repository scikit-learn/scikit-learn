"""Input validation and categorical encoding shared by tree-based estimators.

Used by decision trees, forests, gradient boosting and histogram-based gradient
boosting.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import narwhals.stable.v2 as nw
import numpy as np
from scipy.sparse import issparse

from sklearn.preprocessing import OrdinalEncoder
from sklearn.utils import _safe_indexing
from sklearn.utils._missing import is_scalar_nan
from sklearn.utils.validation import (
    _check_n_features,
    _num_samples,
    check_array,
    validate_data,
)


def _check_categorical_features(X, categorical_features):
    """Check and validate categorical features in X

    Parameters
    ----------
    X : {array-like, pandas DataFrame} of shape (n_samples, n_features)
        Input data.

    categorical_features : array-like of {bool, int, str} of shape (n_features) \
            or shape (n_categorical_features,), default='from_dtype'
        Indicates the categorical features in `X`.

        - None : no feature will be considered categorical.
        - boolean array-like : boolean mask indicating categorical features.
        - integer array-like : integer indices indicating categorical
          features.
        - str array-like: names of categorical features (assuming the training
          data has feature names).
        - `"from_dtype"`: dataframe columns with dtype "Categorical" and "Enum" are
          considered to be categorical features. The input must be a dataframe
          supported by narwhals, meaning :func:`narwhals.from_native` must work.
          This is the case, for instance, for pandas and polars DataFrames.

    Return
    ------
    is_categorical : ndarray of shape (n_features,) or None, dtype=bool
        Indicates whether a feature is categorical. If no feature is
        categorical, this is None.
    """
    if nw.dependencies.is_into_dataframe(X):
        X = nw.from_native(X)
        dtypes = X.schema.dtypes()
        X_is_dataframe = True
        categorical_columns_mask = np.asarray(
            [d in (nw.Categorical, nw.Enum) for d in dtypes]
        )
    else:
        X_is_dataframe = False
        categorical_columns_mask = None

    categorical_by_dtype = (
        isinstance(categorical_features, str) and categorical_features == "from_dtype"
    )
    no_categorical_dtype = categorical_features is None or (
        categorical_by_dtype and not X_is_dataframe
    )

    if no_categorical_dtype:
        return None

    if categorical_by_dtype and X_is_dataframe:
        categorical_features = categorical_columns_mask
    else:
        categorical_features = np.asarray(categorical_features)

    if categorical_features.size == 0:
        return None

    if categorical_features.dtype.kind not in ("i", "b", "U", "O"):
        raise ValueError(
            "categorical_features must be an array-like of bool, int or "
            f"str, got: {categorical_features.dtype.name}."
        )

    if categorical_features.dtype.kind == "O":
        types = set(type(f) for f in categorical_features)
        if types != {str}:
            raise ValueError(
                "categorical_features must be an array-like of bool, int or "
                f"str, got: {', '.join(sorted(t.__name__ for t in types))}."
            )

    n_features = X.shape[1]
    # At this point `validate_data` was not called yet because we use the original
    # dtypes to discover the categorical features. Thus `feature_names_in_`
    # is not defined yet.
    feature_names_in_ = getattr(X, "columns", None)

    if categorical_features.dtype.kind in ("U", "O"):
        # check for feature names
        if feature_names_in_ is None:
            raise ValueError(
                "categorical_features should be passed as an array of "
                "integers or as a boolean mask when the model is fitted "
                "on data without feature names."
            )
        is_categorical = np.zeros(n_features, dtype=bool)
        feature_names = list(feature_names_in_)
        for feature_name in categorical_features:
            try:
                is_categorical[feature_names.index(feature_name)] = True
            except ValueError as e:
                raise ValueError(
                    f"categorical_features has an item value '{feature_name}' "
                    "which is not a valid feature name of the training "
                    f"data. Observed feature names: {feature_names}"
                ) from e
    elif categorical_features.dtype.kind == "i":
        # check for categorical features as indices
        if (
            np.max(categorical_features) >= n_features
            or np.min(categorical_features) < 0
        ):
            raise ValueError(
                "categorical_features set as integer "
                "indices must be in [0, n_features - 1]"
            )
        is_categorical = np.zeros(n_features, dtype=bool)
        is_categorical[categorical_features] = True
    else:
        if categorical_features.shape[0] != n_features:
            raise ValueError(
                "categorical_features set as a boolean mask "
                "must have shape (n_features,), got: "
                f"{categorical_features.shape}"
            )
        is_categorical = categorical_features

    if not np.any(is_categorical):
        return None
    return is_categorical


def _encode_categorical_features(estimator, X, *, reset, dtype):
    """Ordinal-encode the categorical features of X, keeping the column order.

    At fit time (``reset=True``), also detect the categorical features and fit
    the encoder, stored in ``estimator._categorical_encoder`` (None if there is
    no categorical feature). X is returned unchanged if there is no categorical
    feature. Otherwise, a Fortran-ordered array of dtype ``dtype`` is returned.
    """
    if reset:
        estimator._categorical_encoder = None
        # RandomTreesEmbedding has no `categorical_features` parameter.
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
