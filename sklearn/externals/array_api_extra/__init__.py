"""Extra array functions built on top of the array API standard."""

from . import testing
from ._delegation import (
    argpartition,
    atleast_nd,
    broadcast_shapes,
    cov,
    create_diagonal,
    diag_indices,
    expand_dims,
    isclose,
    isin,
    kron,
    nan_to_num,
    nanmax,
    nanmin,
    one_hot,
    pad,
    partition,
    searchsorted,
    setdiff1d,
    sinc,
    tril_indices,
    triu_indices,
    union1d,
    unravel_index,
)
from ._lib._at import at
from ._lib._funcs import (
    angle,
    apply_where,
    default_dtype,
    nunique,
)
from ._lib._lazy import lazy_apply

__version__ = "0.11.1.dev0"

# pylint: disable=duplicate-code
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
    "diag_indices",
    "expand_dims",
    "isclose",
    "isin",
    "kron",
    "lazy_apply",
    "nan_to_num",
    "nanmax",
    "nanmin",
    "nunique",
    "one_hot",
    "pad",
    "partition",
    "searchsorted",
    "setdiff1d",
    "sinc",
    "testing",
    "tril_indices",
    "triu_indices",
    "union1d",
    "unravel_index",
]
