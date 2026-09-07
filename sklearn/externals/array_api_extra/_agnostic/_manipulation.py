"""Array-agnostic implementations for manipulation functions."""

import math
import typing
from collections.abc import Sequence

from .._at import at
from .._lib import _compat, _helpers
from .._lib._typing import Array, ArrayNamespace

__all__ = ["atleast_nd", "broadcast_shapes", "expand_dims", "pad"]


def atleast_nd(x: Array, /, *, ndim: int, xp: ArrayNamespace) -> Array:
    # numpydoc ignore=PR01,RT01
    """See docstring in `array_api_extra._manipulation`."""

    if x.ndim < ndim:
        x = xp.expand_dims(x, axis=0)
        x = atleast_nd(x, ndim=ndim, xp=xp)
    return x


# `float` in signature to accept `math.nan` for Dask.
# `int`s are still accepted as `float` is a superclass of `int` in typing
def broadcast_shapes(  # numpydoc ignore=PR01,RT01
    *shapes: tuple[float | None, ...],
) -> tuple[int | None, ...]:
    """See docstring in `array_api_extra._manipulation`."""
    if not shapes:
        return ()  # Match NumPy output

    ndim = max(len(shape) for shape in shapes)
    out: list[int | None] = []
    for axis in range(-ndim, 0):
        sizes = {shape[axis] for shape in shapes if axis >= -len(shape)}
        # Dask uses NaN for unknown shape, which predates the Array API spec for None
        none_size = None in sizes or math.nan in sizes  # noqa: PLW0177
        sizes -= {1, None, math.nan}
        if len(sizes) > 1:
            msg = (
                "shape mismatch: objects cannot be broadcast to a single shape: "
                f"{shapes}."
            )
            raise ValueError(msg)
        out.append(None if none_size else typing.cast(int, sizes.pop()) if sizes else 1)

    return tuple(out)


def expand_dims(
    a: Array, /, *, axis: tuple[int, ...] = (0,), xp: ArrayNamespace
) -> Array:
    # numpydoc ignore=PR01,RT01
    """See docstring in `array_api_extra._manipulation`."""
    for i in sorted(axis):
        a = xp.expand_dims(a, axis=i)
    return a


def pad(
    x: Array,
    pad_width: int | tuple[int, int] | Sequence[tuple[int, int]],
    *,
    constant_values: complex = 0,
    xp: ArrayNamespace,
) -> Array:  # numpydoc ignore=PR01,RT01
    """See docstring in `array_api_extra._manipulation`."""
    pad_width_seq = _helpers.normalize_pad_width(pad_width, x.ndim)

    slices: list[slice] = []
    newshape: list[int] = []
    for ax, w_tpl in enumerate(pad_width_seq):
        if len(w_tpl) != 2:
            msg = f"expect a 2-tuple (before, after), got {w_tpl}."
            raise ValueError(msg)

        sh = _helpers.eager_shape(x)[ax]

        if w_tpl[0] == 0 and w_tpl[1] == 0:
            sl = slice(None, None, None)
        else:
            stop: int | None
            start, stop = w_tpl
            stop = None if stop == 0 else -stop

            sl = slice(start, stop, None)
            sh += w_tpl[0] + w_tpl[1]

        newshape.append(sh)
        slices.append(sl)

    padded = xp.full(
        tuple(newshape),
        fill_value=constant_values,
        dtype=x.dtype,
        device=_compat.device(x),
    )
    return at(padded, tuple(slices)).set(x)
