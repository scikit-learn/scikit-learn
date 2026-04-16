"""Extra array functions built on top of the array API standard."""

from ._delegation import (
    argpartition,
    atleast_nd,
    cov,
    create_diagonal,
    expand_dims,
    isclose,
    isin,
    nan_to_num,
    one_hot,
    pad,
    partition,
    setdiff1d,
    sinc,
)
from ._lib._at import at
from ._lib._funcs import (
    apply_where,
    broadcast_shapes,
    default_dtype,
    kron,
    nunique,
)
from ._lib._lazy import lazy_apply

__version__ = "0.9.1"

# pylint: disable=duplicate-code
__all__ = [
    "__version__",
    "apply_where",
    "argpartition",
    "at",
    "atleast_nd",
    "broadcast_shapes",
    "cov",
    "create_diagonal",
    "default_dtype",
    "expand_dims",
    "isclose",
    "isin",
    "kron",
    "lazy_apply",
    "nan_to_num",
    "nunique",
    "one_hot",
    "pad",
    "partition",
    "setdiff1d",
    "sinc",
]
