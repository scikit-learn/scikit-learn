"""Control sparse interface based on config"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import scipy as sp

from .._config import get_config
from .fixes import parse_version, sp_base_version


def _align_api_if_sparse(X):
    """
    Convert to sparse interface as set in config. Input can be dense or sparse.
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


########### fixes for transitioning function names

# TODO: Replace when Scipy 1.12 is the minimum supported version
if sp_base_version >= parse_version("1.12.0"):
    _sparse_eye = sp.sparse.eye_array
    _sparse_block = sp.sparse.block_array
    _sparse_diags = sp.sparse.diags_array
    _sparse_random = sp.sparse.random_array
else:

    def _sparse_eye(m, n=None, *, k=0, dtype=float, format=None):
        return sp.sparse.eye(m, n, k=k, dtype=dtype, format=format)

    def _sparse_block(blocks, *, format=None, dtype=None):
        return sp.sparse.bmat(blocks, format=format, dtype=dtype)

    def _sparse_diags(diagonals, /, *, offsets=0, shape=None, format=None, dtype=None):
        return sp.sparse.diags(
            diagonals, offsets=offsets, shape=shape, format=format, dtype=dtype
        )

    def _sparse_random(
        shape,
        *,
        density=0.01,
        format="coo",
        dtype=None,
        random_state=None,
        rng=None,
        data_sampler=None,
    ):
        return sp.sparse.random(
            *shape,
            density=density,
            format=format,
            dtype=dtype,
            random_state=random_state,
            rng=rng,
            data_rvs=data_sampler,
        )
