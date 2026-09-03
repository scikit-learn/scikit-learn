"""Array-agnostic implementations for sorting functions."""

from .._lib._typing import Array, ArrayNamespace

__all__ = ["argpartition", "partition"]


def partition(  # numpydoc ignore=PR01,RT01
    x: Array,
    kth: int,  # noqa: ARG001
    /,
    axis: int = -1,
    *,
    xp: ArrayNamespace,
) -> Array:
    """See docstring in `array_api_extra._sorting`."""
    return xp.sort(x, axis=axis, stable=False)


def argpartition(  # numpydoc ignore=PR01,RT01
    x: Array,
    kth: int,  # noqa: ARG001
    /,
    axis: int = -1,
    *,
    xp: ArrayNamespace,
) -> Array:
    """See docstring in `array_api_extra._sorting`."""
    return xp.argsort(x, axis=axis, stable=False)
