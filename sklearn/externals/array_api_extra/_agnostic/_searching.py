"""Array-agnostic implementations for searching functions."""

import math
from typing import Literal

from .._lib import _compat
from .._lib._typing import Array, ArrayNamespace
from . import _inspection

__all__ = ["searchsorted"]


def searchsorted(
    x1: Array,
    x2: Array,
    /,
    *,
    side: Literal["left", "right"] = "left",
    xp: ArrayNamespace,
) -> Array:
    # numpydoc ignore=PR01,RT01
    """See docstring in `array_api_extra._searching`."""
    a = xp.full(x2.shape, 0, device=_compat.device(x1))

    if x1.shape[-1] == 0:
        return a

    n = xp.count_nonzero(~xp.isnan(x1), axis=-1, keepdims=True)
    b = xp.broadcast_to(n, x2.shape)

    compare = xp.less_equal if side == "left" else xp.less

    # while xp.any(b - a > 1):
    # refactored to for loop with ~log2(n) iterations for JAX JIT
    for _ in range(int(math.log2(x1.shape[-1])) + 1):  # type: ignore[arg-type]  # pyright: ignore[reportArgumentType]
        c = (a + b) // 2
        x0 = xp.take_along_axis(x1, c, axis=-1)
        j = compare(x2, x0)
        b = xp.where(j, c, b)
        a = xp.where(j, a, c)

    out = xp.where(compare(x2, xp.min(x1, axis=-1, keepdims=True)), 0, b)
    out = xp.where(xp.isnan(x2), x1.shape[-1], out) if side == "right" else out
    return xp.astype(out, _inspection.default_dtype(xp, kind="integral"), copy=False)
