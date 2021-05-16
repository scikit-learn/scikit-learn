from functools import wraps
from functools import partial
from inspect import signature

import scipy.sparse as sp_sparse


def _get_feature_names(X):
    """Get feature names of a dataframe."""
    if hasattr(X, "iloc"):  # duck-type dataframe
        return getattr(X, "columns", None)


def _get_index(X):
    """Get the index of a dataframe."""
    if hasattr(X, "iloc"):  # duck-type dataframe
        return getattr(X, "index", None)


def _check_get_feature_names_out(get_feature_names_out, estimator,
                                 n_feature_out):
    """Check and convert get_feature_names_out into a callable."""
    if callable(get_feature_names_out):
        get_feature_names_out_callable = get_feature_names_out(estimator)
    elif get_feature_names_out == 'one_to_one':
        def get_feature_names_out_callable(names):
            return names
    else:
        # get_feature_names_out == 'class_name'
        class_name = estimator.__class__.__name__.lower()

        def get_feature_names_out_callable():
            return [f"{class_name}{i}" for i in range(n_feature_out)]

    # feature names in can have zero or one argument. For one argument
    # it would be the input feature names
    parameters = signature(get_feature_names_out_callable).parameters

    if parameters:
        feature_names_in = getattr(estimator, "feature_names_in_", None)
        get_feature_names_out_callable = partial(
            get_feature_names_out_callable, feature_names_in)
    return get_feature_names_out_callable


def _make_array_out(X_out, index, get_feature_names_out, *,
                    array_out="default"):
    """Construct array container based on global configuration.

    Parameters
    ----------
    X_out: {ndarray, sparse matrix} of shape (n_samples, n_features_out)
        Output data to be wrapped.

    index: array-like of shape (n_samples,)
        Index of output data.

    get_features_names_out: callable
        Returns the feature names out. If the callable returns None, then
        the feature names will be ["X0", "X1", ...].

    array_out : {"default", "pandas"}, default="default"
        Specify the output array type. If "pandas", a pandas DataFrame is
        returned. If "default", an array-like without feature names is
        returned.

    Return
    ------
    array_out: {ndarray, sparse matrix, dataframe} of shape \
                (n_samples, n_features_out)
        Wrapped array with feature names.
    """
    if array_out not in {'default', 'pandas'}:
        raise ValueError("array_out must be 'default' or 'pandas'")

    if array_out == "default":
        return X_out

    feature_names_out = get_feature_names_out()
    if feature_names_out is None:
        feature_names_out = [f'X{i}' for i in range(X_out.shape[1])]

    # array_out == "pandas"
    import pandas as pd
    if sp_sparse.issparse(X_out):
        make_dataframe = pd.DataFrame.sparse.from_spmatrix
    else:
        make_dataframe = pd.DataFrame

    return make_dataframe(X_out, columns=feature_names_out, index=index)


def _array_out_wrap(get_feature_names_out):
    """Wrap around transform method to create array_out.

    Parameters
    ----------
    get_feature_names_out : callable or {"one_to_one", "class_name"}
        Called to get the feature names out. If `one_to_one`, then the
        feature_names_in will be used as feature name out. If `class_name`,
        then the class name will be used as prefixes for the feature names
        out.
    """
    def _wrapper_transform(transform):

        @wraps(transform)
        def inner_transform(*args, **kwargs):
            array_out = kwargs.get("array_out", "default")
            X_out = transform(*args, **kwargs)

            if array_out == "default":
                return X_out

            estimator, X_orig = args[0], args[1]
            index = _get_index(X_orig)
            get_features = _check_get_feature_names_out(
                get_feature_names_out, estimator, X_out.shape[1])

            return _make_array_out(X_out, index, get_features,
                                   array_out=array_out)
        return inner_transform
    return _wrapper_transform
