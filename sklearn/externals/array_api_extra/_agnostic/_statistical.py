"""Array-agnostic implementations for statistical functions."""

import math
import typing
import warnings

from .._lib import _compat, _helpers
from .._lib._typing import Array, ArrayNamespace
from . import _manipulation

__all__ = ["cov", "nanmax", "nanmean", "nanmin", "nansum"]


def cov(
    m: Array,
    /,
    *,
    correction: float = 1,
    fweights: Array | None = None,
    aweights: Array | None = None,
    xp: ArrayNamespace,
) -> Array:  # numpydoc ignore=PR01,RT01
    """See docstring in `array_api_extra._statistical`."""
    # NB: no `xp.asarray(m)` here. The delegation layer already guarantees `m`
    # is an array (it calls `array_namespace(m)` and reads `m.ndim`), and on
    # torch `xp.asarray` detaches gradients and mutates the caller's tensor.
    dtype = (
        xp.float64 if xp.isdtype(m.dtype, "integral") else xp.result_type(m, xp.float64)
    )

    m = _manipulation.atleast_nd(m, ndim=2, xp=xp)
    # Preserve the historical no-alias guarantee even when the dtype already matches.
    m = xp.astype(m, dtype, copy=True)

    # Validate weight shapes (eager metadata, lazy-safe).
    n_obs = m.shape[-1]
    for name, w_in in (("fweights", fweights), ("aweights", aweights)):
        if w_in is None:
            continue
        if w_in.ndim != 1:
            msg = f"`{name}` must be 1-D, got ndim={w_in.ndim}"
            raise ValueError(msg)
        weight_length = w_in.shape[0]
        # Unknown dims are `None` per the standard; Dask non-standardly
        # reports them as NaN, hence the `isnan` checks below.
        if (
            weight_length is not None
            and n_obs is not None
            and not math.isnan(weight_length)
            and not math.isnan(n_obs)
            and weight_length != n_obs
        ):
            msg = (
                f"`{name}` has length {weight_length} but `m` has {n_obs} observations"
            )
            raise ValueError(msg)

    fw = None
    if fweights is not None:
        fw = xp.astype(xp.asarray(fweights), dtype)
    aw = None
    if aweights is not None:
        aw = xp.astype(xp.asarray(aweights), dtype)
    if fw is None and aw is None:
        w = None
    elif fw is None:
        w = aw
    elif aw is None:
        w = fw
    else:
        w = fw * aw

    if w is None:
        avg = xp.mean(m, axis=-1, keepdims=True)
        fact = _helpers.eager_shape(m, axis=-1)[0] - correction
    else:
        v1 = xp.sum(w, axis=-1)
        avg = xp.sum(m * w, axis=-1, keepdims=True) / v1
        if aw is None:
            fact = v1 - correction
        else:
            fact = v1 - correction * xp.sum(w * aw, axis=-1) / v1

    if not _compat.is_lazy_array(fact):
        # Weights are cast to `dtype`, so a complex input produces a complex
        # normalizer with a zero imaginary part. Complex ordering is undefined;
        # compare its real component instead.
        if w is not None:
            fact_array = typing.cast(Array, fact)
            fact_to_check = (
                xp.real(fact_array)
                if xp.isdtype(fact_array.dtype, "complex floating")
                else fact_array
            )
        else:
            fact_to_check = fact
        if fact_to_check <= 0:
            warnings.warn(
                "Degrees of freedom <= 0 for slice", RuntimeWarning, stacklevel=2
            )
            fact = 0

    m_c = m - avg
    m_w = m_c if w is None else m_c * w
    m_cT = xp.matrix_transpose(m_c)
    if xp.isdtype(m_cT.dtype, "complex floating"):
        m_cT = xp.conj(m_cT)
    c = m_w @ m_cT / fact
    axes = tuple(axis for axis, length in enumerate(c.shape) if length == 1)
    return xp.squeeze(c, axis=axes)


def nanmin(  # numpydoc ignore=PR01,RT01
    a: Array,
    /,
    *,
    axis: int | tuple[int, ...] | None,
    xp: ArrayNamespace,
) -> Array:
    """See docstring in `array_api_extra._statistical`."""
    mask = xp.isnan(a)
    x = xp.min(xp.where(mask, xp.full_like(a, +xp.inf), a), axis=axis)
    # Replace Infs from all NaN slices with NaN again
    mask = xp.all(mask, axis=axis)
    if xp.any(mask):
        x = xp.where(mask, xp.full_like(x, xp.nan), x)
    return x


def nanmax(  # numpydoc ignore=PR01,RT01
    a: Array,
    /,
    *,
    axis: int | tuple[int, ...] | None,
    xp: ArrayNamespace,
) -> Array:
    """See docstring in `array_api_extra._statistical`."""
    mask = xp.isnan(a)
    x = xp.max(xp.where(mask, xp.full_like(a, -xp.inf), a), axis=axis)
    # Replace Infs from all NaN slices with NaN again
    mask = xp.all(mask, axis=axis)
    if xp.any(mask):
        x = xp.where(mask, xp.full_like(x, xp.nan), x)
    return x


def nansum(  # numpydoc ignore=PR01,RT01
    a: Array,
    /,
    *,
    axis: int | tuple[int, ...] | None,
    xp: ArrayNamespace,
) -> Array:
    """See docstring in `array_api_extra._statistical`."""
    mask = xp.isnan(a)
    return xp.sum(xp.where(mask, xp.zeros_like(a), a), axis=axis)


def nanmean(  # numpydoc ignore=PR01,RT01
    a: Array,
    /,
    *,
    axis: int | tuple[int, ...] | None,
    xp: ArrayNamespace,
) -> Array:
    """See docstring in `array_api_extra._statistical`."""
    mask = xp.isnan(a)
    sum_ = nansum(a, axis=axis, xp=xp)
    count = xp.sum(xp.where(mask, xp.zeros_like(a), xp.ones_like(a)), axis=axis)
    safe_count = xp.astype(
        xp.where(count == 0, xp.ones_like(count), count),
        sum_.dtype,
        copy=False,
    )
    result = sum_ / safe_count
    return xp.where(
        count == 0,
        xp.full_like(result, xp.nan),
        result,
    )
