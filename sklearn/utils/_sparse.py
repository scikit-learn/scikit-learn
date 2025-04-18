"""Control sparse interface based on config"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import scipy as sp

from .._config import get_config


def _as_sparse(X_sparse):
    """
    Convert to sparse interface as set in config. Input must be sparse.
    If you know the input is sparse, use e.g. `return _as_sparse(X)`.
    Otherwise you should check if sparse before calling.

        if sp.sparse.issparse(X):
            X = _as_sparse(X)
    """
    if not sp.sparse.issparse(X_sparse):
        raise TypeError("Input should be a SciPy sparse container")

    return _convert_sparse_to_config_chosen_interface(X_sparse)


def _select_interface_if_sparse(X):
    """
    Convert to sparse interface as set in config. Input X can be dense or sparse.
    If sparse, convert to sparse_interface indicated by get_config.
    Otherwise, return X unchanged.

        X = _select_interface_if_sparse(X)
    """
    if not sp.sparse.issparse(X):
        return X

    return _convert_sparse_to_config_chosen_interface(X)


def _convert_sparse_to_config_chosen_interface(X_sparse):
    # assume there are only two sparse interfaces: sparray and spmatrix
    X_is_sparray = not sp.sparse.isspmatrix(X_sparse)
    config_sparse_interface = get_config()["sparse_interface"]

    if config_sparse_interface == "sparray":
        if X_is_sparray:
            return X_sparse
        return _convert_from_spmatrix_to_sparray(X_sparse)
    else:  # global is spmatrix
        if not X_is_sparray:
            return X_sparse
        return _convert_from_sparray_to_spmatrix(X_sparse)


def _convert_from_spmatrix_to_sparray(X_sparse):
    """Fundamental code to switch to sparray in any format"""
    return getattr(sp.sparse, X_sparse.format + "_array")(X_sparse)


def _convert_from_sparray_to_spmatrix(X_sparse):
    """Fundamental code to switch to spmatrix in any format"""
    return getattr(sp.sparse, X_sparse.format + "_matrix")(X_sparse)
