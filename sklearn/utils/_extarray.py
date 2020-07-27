import numpy as np


def _get_feature_names(X):
    """Get feature names from X.

    Supports:
       - pandas DataFrame
       - xarray DataArray
       - Return None for unrecognized array containers
    """
    if hasattr(X, "columns"):  # pandas
        return np.array(X.columns, dtype=object)
    elif hasattr(X, "dims") and isinstance(X.dims, tuple) and len(X.dims) == 2:
        # xarray DataArray
        return np.array(X.coords[X.dims[1]], dtype=object)
