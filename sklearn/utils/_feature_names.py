import numpy as np


def _get_feature_names(X):
    """Get feature names from X.

    Parameters
    ----------
    X : {dataframe, dataarray} of shape (n_samples, n_features)
        Array container to extract feature names.

        - pandas DataFrame : The columns will be considered to be feature
          names.
        - xarray DataArray : The coords of the second dimension will be
          considered to be feature names.
        - All other array containers will return `None`.

    Returns
    -------
    names: ndarray of shape (n_features,) or None
        Column names of `X`. Unrecognized array containers will return `None`.

    Raises
    ------
    ValueError
        If column names consist of a non-string data type.
    """
    if hasattr(X, "columns"):
        # pandas
        return np.array(X.columns, dtype=object)
    elif hasattr(X, "dims") and isinstance(X.dims, tuple) and len(X.dims) == 2:
        # xarray DataArray
        return np.array(X.coords[X.dims[1]], dtype=object)
