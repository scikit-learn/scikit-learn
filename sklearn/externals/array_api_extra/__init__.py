"""Extra array functions built on top of the array API standard."""

from . import testing
from ._agnostic._elementwise import angle, apply_where
from ._agnostic._inspection import default_dtype
from ._at import at
from ._creation import create_diagonal, one_hot
from ._elementwise import deg2rad, isclose, nan_to_num, rad2deg, sinc
from ._indexing import diag_indices, tril_indices, triu_indices, unravel_index
from ._lazy import lazy_apply
from ._linalg import kron
from ._manipulation import atleast_nd, broadcast_shapes, expand_dims, pad
from ._searching import searchsorted
from ._set import isin, nunique, setdiff1d, union1d
from ._sorting import argpartition, partition
from ._statistical import cov, nanmax, nanmean, nanmin, nansum

__version__ = "0.11.2"

__all__ = [
    "__version__",
    "angle",
    "apply_where",
    "argpartition",
    "at",
    "atleast_nd",
    "broadcast_shapes",
    "cov",
    "create_diagonal",
    "default_dtype",
    "deg2rad",
    "diag_indices",
    "expand_dims",
    "isclose",
    "isin",
    "kron",
    "lazy_apply",
    "nan_to_num",
    "nanmax",
    "nanmean",
    "nanmin",
    "nansum",
    "nunique",
    "one_hot",
    "pad",
    "partition",
    "rad2deg",
    "searchsorted",
    "setdiff1d",
    "sinc",
    "testing",
    "tril_indices",
    "triu_indices",
    "union1d",
    "unravel_index",
]
