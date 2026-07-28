"""Extra array functions built on top of the array API standard."""

from . import testing
from ._delegation import (
    argpartition,
    atleast_nd,
    broadcast_shapes,
    cov,
    create_diagonal,
    expand_dims,
    isclose,
    isin,
    kron,
    nan_to_num,
    one_hot,
    pad,
    partition,
    searchsorted,
    setdiff1d,
    sinc,
    union1d,
)
from ._lib._at import at
from ._lib._funcs import (
    angle,
    apply_where,
    default_dtype,
    nunique,
)
from ._lib._lazy import lazy_apply

__version__ = "0.10.3"

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
    "searchsorted",
    "setdiff1d",
    "sinc",
    "testing",
    "union1d",
]
