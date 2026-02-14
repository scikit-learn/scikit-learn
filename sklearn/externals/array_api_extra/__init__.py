"""Extra array functions built on top of the array API standard."""

from ._delegation import isclose, nan_to_num, one_hot, pad
from ._lib._at import at
from ._lib._funcs import (
    apply_where,
    atleast_nd,
    broadcast_shapes,
    cov,
    create_diagonal,
    default_dtype,
    expand_dims,
    kron,
    nunique,
    setdiff1d,
    sinc,
)
from ._lib._lazy import lazy_apply

__version__ = "0.8.2"

# pylint: disable=duplicate-code
__all__ = [
    "__version__",
    "apply_where",
    "at",
    "atleast_nd",
    "broadcast_shapes",
    "cov",
    "create_diagonal",
    "default_dtype",
    "expand_dims",
    "isclose",
    "kron",
    "lazy_apply",
    "nan_to_num",
    "nunique",
    "one_hot",
    "pad",
    "setdiff1d",
    "sinc",
]
