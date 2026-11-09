"""Array-agnostic implementations for linear algebra functions."""

import typing

from .._lib import _helpers
from .._lib._typing import Array, ArrayNamespace
from . import _manipulation

__all__ = ["kron"]


def kron(
    a: Array,
    b: Array,
    /,
    *,
    xp: ArrayNamespace,
) -> Array:  # numpydoc ignore=PR01,RT01
    """See docstring in `array_api_extra._linalg`."""

    singletons = (1,) * (b.ndim - a.ndim)
    a = typing.cast(Array, xp.broadcast_to(a, singletons + a.shape))

    nd_b, nd_a = b.ndim, a.ndim
    nd_max = max(nd_b, nd_a)
    if nd_a == 0 or nd_b == 0:
        return xp.multiply(a, b)

    a_shape = _helpers.eager_shape(a)
    b_shape = _helpers.eager_shape(b)

    # Equalise the shapes by prepending smaller one with 1s
    a_shape = (1,) * max(0, nd_b - nd_a) + a_shape
    b_shape = (1,) * max(0, nd_a - nd_b) + b_shape

    # Insert empty dimensions
    a_arr = _manipulation.expand_dims(a, axis=tuple(range(nd_b - nd_a)), xp=xp)
    b_arr = _manipulation.expand_dims(b, axis=tuple(range(nd_a - nd_b)), xp=xp)

    # Compute the product
    a_arr = _manipulation.expand_dims(a_arr, axis=tuple(range(1, nd_max * 2, 2)), xp=xp)
    b_arr = _manipulation.expand_dims(b_arr, axis=tuple(range(0, nd_max * 2, 2)), xp=xp)
    result = xp.multiply(a_arr, b_arr)

    # Reshape back and return
    res_shape = tuple(a_s * b_s for a_s, b_s in zip(a_shape, b_shape, strict=True))
    return xp.reshape(result, res_shape)
