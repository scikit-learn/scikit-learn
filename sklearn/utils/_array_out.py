import numpy as np
import scipy.sparse as sp_sparse
from .._config import get_config


def _get_feature_names(X):
    """Get feature names of a dataframe or a dataarray."""
    if hasattr(X, "columns"):
        # pandas
        return np.array(X.columns, dtype=object)
    elif hasattr(X, "dims") and isinstance(X.dims, tuple) and len(X.dims) == 2:
        # xarray DataArray
        return np.array(X.coords[X.dims[1]], dtype=object)


def _make_array_out(X_out, X_orig, get_feature_names_out):
    """Construct array container based on global configuration.

    Parameters
    ----------
    X_out: {ndarray, sparse matrix} of shape (n_samples, n_features_out)
        Output data to be wrapped.

    X_orig: array-like of shape (n_samples, n_features)
        Original input data. For panda's DataFrames, this is used to get
        the index. For xarray's DataArrays, this is used to get the name
        of the dims and the coordinates for the first dims.

    get_features_names_out: callable
        Returns the feature names out. If the callable returns None, then
        the feature names will be ["X0", "X1", ...].

    Return
    ------
    array_out: {ndarray, sparse matrix, dataframe, dataarray} of shape \
                (n_samples, n_features_out)
        Wrapped array with feature names.
    """
    array_out = get_config()['array_out']
    if array_out not in ('default', 'pandas', 'xarray'):
        raise ValueError("array_out must be 'default', 'pandas' or 'xarray'")
    if array_out == 'default':
        return X_out

    feature_names_out = get_feature_names_out()
    if feature_names_out is None:
        feature_names_out = [f'X{i}' for i in range(X_out.shape[1])]

    if array_out == 'pandas':
        import pandas as pd
        if sp_sparse.issparse(X_out):
            make_dataframe = pd.DataFrame.sparse.from_spmatrix
        else:
            make_dataframe = pd.DataFrame

        return make_dataframe(X_out, columns=feature_names_out,
                              index=getattr(X_orig, "index", None))
    else:
        # xarray
        import xarray as xr
        dims = getattr(X_orig, "dims", ("index", "columns"))

        coords = {dims[1]: feature_names_out}
        if hasattr(X_orig, "coords") and dims[0] in X_orig.coords:
            coords[dims[0]] = X_orig.coords[dims[0]]

        if sp_sparse.issparse(X_out):
            # pydata/sparse
            import sparse as pydata_sparse
            X_out = pydata_sparse.COO.from_scipy_sparse(X_out)

        return xr.DataArray(X_out, dims=dims, coords=coords)
