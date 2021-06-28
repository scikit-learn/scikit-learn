import numpy as np


def _get_feature_names(X):
    """Get feature names from X.

    Parameters
    ----------
    X : {dataframe} of shape (n_samples, n_features)
        Array container to extract feature names.

        - pandas DataFrame : The columns will be considered to be feature
          names.
        - All other array containers will return `None`.

    Returns
    -------
    names: ndarray of shape (n_features,) or None
        Column names of `X`. Unrecognized array containers will return `None`.
    """
    if hasattr(X, "columns"):
        # pandas
        return np.array(X.columns, dtype=object)
