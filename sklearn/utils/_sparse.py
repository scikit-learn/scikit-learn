"""Control sparse interface based on config"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import scipy as sp

from sklearn._config import get_config


def _align_api_if_sparse(X):
    """
    Convert to sparse interface as set in config.

    Input can be dense or sparse.
    If sparse, convert to sparse_interface indicated by get_config.
    Otherwise, return X unchanged.
    """
    if not sp.sparse.issparse(X):
        return X

    config_sparse_interface = get_config()["sparse_interface"]

    # there are only two sparse interfaces: sparray and spmatrix
    if config_sparse_interface == "sparray":
        if sp.sparse.isspmatrix(X):
            # Fundamental code to switch to sparray in any format
            return getattr(sp.sparse, X.format + "_array")(X)
        return X
    else:  # config is spmatrix
        if sp.sparse.isspmatrix(X):
            return X
        # Fundamental code to switch to spmatrix in any format
        return getattr(sp.sparse, X.format + "_matrix")(X)
