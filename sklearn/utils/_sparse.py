"""Control sparse interface based on config"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import scipy as sp

from sklearn._config import get_config
from sklearn.utils.fixes import SCIPY_VERSION_BELOW_1_12


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
if not SCIPY_VERSION_BELOW_1_12:
    _sparse_eye = sp.sparse.eye_array
    _sparse_diags = sp.sparse.diags_array

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
        X = sp.sparse.random_array(
            shape,
            density=density,
            format=format,
            dtype=dtype,
            random_state=rng or random_state,
            data_sampler=data_sampler,
        )
        _ensure_sparse_index_int32(X)
        return X

else:

    def _sparse_eye(m, n=None, *, k=0, dtype=float, format=None):
        return sp.sparse.eye(m, n, k=k, dtype=dtype, format=format)

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
            random_state=rng or random_state,
            data_rvs=data_sampler,
        )


########### fixes for casting index arrays


def _ensure_sparse_index_int32(A):
    """Safely ensure that index arrays are int32."""
    if A.format in ("csc", "csr", "bsr"):
        A.indices, A.indptr = safely_cast_index_arrays(A)
    elif A.format == "coo":
        if hasattr(A, "coords"):
            A.coords = safely_cast_index_arrays(A)
        elif hasattr(A, "indices"):
            A.indices = safely_cast_index_arrays(A)
        else:
            A.row, A.col = safely_cast_index_arrays(A)
    elif A.format == "dia":
        A.offsets = safely_cast_index_arrays(A)


# TODO: remove when SciPy 1.15 is minimal supported version
#       (based on scipy.sparse._sputils.py function with same name)
def safely_cast_index_arrays(A, idx_dtype=np.int32, msg=""):
    """Safely cast sparse array indices to `idx_dtype`.

    Check the shape of `A` to determine if it is safe to cast its index
    arrays to dtype `idx_dtype`. If any dimension in shape is larger than
    fits in the dtype, casting is unsafe so raise ``ValueError``.
    If safe, cast the index arrays to `idx_dtype` and return the result
    without changing the input `A`. The caller can assign results to `A`
    attributes if desired or use the recast index arrays directly.

    Unless downcasting is needed, the original index arrays are returned.
    You can test e.g. ``A.indptr is new_indptr`` to see if downcasting occurred.

    See SciPy: scipy.sparse._sputils.py for more info on safely_cast_index_arrays()
    """
    max_value = np.iinfo(idx_dtype).max

    if A.format in ("csc", "csr"):
        if A.indptr[-1] > max_value:
            raise ValueError(f"indptr values too large for {msg}")
        # check shape vs dtype
        if max(*A.shape) > max_value:
            if (A.indices > max_value).any():
                raise ValueError(f"indices values too large for {msg}")

        indices = A.indices.astype(idx_dtype, copy=False)
        indptr = A.indptr.astype(idx_dtype, copy=False)
        return indices, indptr

    elif A.format == "coo":
        coords = getattr(A, "coords", None)
        if coords is None:
            coords = getattr(A, "indices", None)
            if coords is None:
                coords = (A.row, A.col)
        if max(*A.shape) > max_value:
            if any((co > max_value).any() for co in coords):
                raise ValueError(f"coords values too large for {msg}")
        return tuple(co.astype(idx_dtype, copy=False) for co in coords)

    elif A.format == "dia":
        if max(*A.shape) > max_value:
            if (A.offsets > max_value).any():
                raise ValueError(f"offsets values too large for {msg}")
        offsets = A.offsets.astype(idx_dtype, copy=False)
        return offsets

    elif A.format == "bsr":
        R, C = A.blocksize
        if A.indptr[-1] * R > max_value:
            raise ValueError("indptr values too large for {msg}")
        if max(*A.shape) > max_value:
            if (A.indices * C > max_value).any():
                raise ValueError(f"indices values too large for {msg}")
        indices = A.indices.astype(idx_dtype, copy=False)
        indptr = A.indptr.astype(idx_dtype, copy=False)
        return indices, indptr
    # DOK and LIL formats are not associated with index arrays.
