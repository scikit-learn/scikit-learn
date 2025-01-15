"""Extra array functions built on top of the array API standard."""

from ._funcs import (
    at,
    atleast_nd,
    cov,
    create_diagonal,
    expand_dims,
    kron,
    pad,
    setdiff1d,
    sinc,
)

__version__ = "0.5.0"

# pylint: disable=duplicate-code
__all__ = [
    "__version__",
    "at",
    "atleast_nd",
    "cov",
    "create_diagonal",
    "expand_dims",
    "kron",
    "pad",
    "setdiff1d",
    "sinc",
]
