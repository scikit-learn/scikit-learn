"""Module for working with feature names."""
import numpy as np

import warnings


def _get_feature_names(X):
    """Get feature names from X.

    Support for other array containers should place its implementation here.

    Parameters
    ----------
    X : {ndarray, dataframe} of shape (n_samples, n_features)
        Array container to extract feature names.

        - pandas dataframe : The columns will be considered to be feature
          names. If the dataframe contains non-string feature names, `None` is
          returned.
        - All other array containers will return `None`.

    Returns
    -------
    names: ndarray or None
        Feature names of `X`. Unrecognized array containers will return `None`.
    """
    if hasattr(X, "columns"):
        feature_names = np.asarray(X.columns)
        # Only strings are supported
        if not all(isinstance(item, str) for item in feature_names):
            warnings.warn(
                "Feature name support requires all feature names to be strings. "
                "Passing non-str feature names will raise an error in 1.2",
                FutureWarning,
            )

            return
        return feature_names
