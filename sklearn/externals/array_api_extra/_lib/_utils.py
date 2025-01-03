"""Utility functions used by `array_api_extra/_funcs.py`."""

# https://github.com/scikit-learn/scikit-learn/pull/27910#issuecomment-2568023972
from __future__ import annotations

from . import _compat
from ._typing import Array, ModuleType

__all__ = ["in1d", "mean"]


def in1d(
    x1: Array,
    x2: Array,
    /,
    *,
    assume_unique: bool = False,
    invert: bool = False,
    xp: ModuleType | None = None,
) -> Array:  # numpydoc ignore=PR01,RT01
    """
    Check whether each element of an array is also present in a second array.

    Returns a boolean array the same length as `x1` that is True
    where an element of `x1` is in `x2` and False otherwise.

    This function has been adapted using the original implementation
    present in numpy:
    https://github.com/numpy/numpy/blob/v1.26.0/numpy/lib/arraysetops.py#L524-L758
    """
    if xp is None:
        xp = _compat.array_namespace(x1, x2)

    # This code is run to make the code significantly faster
    if x2.shape[0] < 10 * x1.shape[0] ** 0.145:
        if invert:
            mask = xp.ones(x1.shape[0], dtype=xp.bool, device=_compat.device(x1))
            for a in x2:
                mask &= x1 != a
        else:
            mask = xp.zeros(x1.shape[0], dtype=xp.bool, device=_compat.device(x1))
            for a in x2:
                mask |= x1 == a
        return mask

    rev_idx = xp.empty(0)  # placeholder
    if not assume_unique:
        x1, rev_idx = xp.unique_inverse(x1)
        x2 = xp.unique_values(x2)

    ar = xp.concat((x1, x2))
    device_ = _compat.device(ar)
    # We need this to be a stable sort.
    order = xp.argsort(ar, stable=True)
    reverse_order = xp.argsort(order, stable=True)
    sar = xp.take(ar, order, axis=0)
    if sar.size >= 1:
        bool_ar = sar[1:] != sar[:-1] if invert else sar[1:] == sar[:-1]
    else:
        bool_ar = xp.asarray([False]) if invert else xp.asarray([True])
    flag = xp.concat((bool_ar, xp.asarray([invert], device=device_)))
    ret = xp.take(flag, reverse_order, axis=0)

    if assume_unique:
        return ret[: x1.shape[0]]
    return xp.take(ret, rev_idx, axis=0)


def mean(
    x: Array,
    /,
    *,
    axis: int | tuple[int, ...] | None = None,
    keepdims: bool = False,
    xp: ModuleType | None = None,
) -> Array:  # numpydoc ignore=PR01,RT01
    """
    Complex mean, https://github.com/data-apis/array-api/issues/846.
    """
    if xp is None:
        xp = _compat.array_namespace(x)

    if xp.isdtype(x.dtype, "complex floating"):
        x_real = xp.real(x)
        x_imag = xp.imag(x)
        mean_real = xp.mean(x_real, axis=axis, keepdims=keepdims)
        mean_imag = xp.mean(x_imag, axis=axis, keepdims=keepdims)
        return mean_real + (mean_imag * xp.asarray(1j))
    return xp.mean(x, axis=axis, keepdims=keepdims)
