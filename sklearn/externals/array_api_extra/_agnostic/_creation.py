"""Array-agnostic implementations for creation functions."""

from .._at import at
from .._lib import _compat, _helpers
from .._lib._typing import Array, ArrayNamespace

__all__ = ["create_diagonal", "one_hot"]


def create_diagonal(
    x: Array, /, *, offset: int = 0, xp: ArrayNamespace
) -> Array:  # numpydoc ignore=PR01,RT01
    """See docstring in `array_api_extra._creation`."""
    x_shape = _helpers.eager_shape(x)
    batch_dims = x_shape[:-1]
    n = x_shape[-1] + abs(offset)
    diag = xp.zeros((*batch_dims, n**2), dtype=x.dtype, device=_compat.device(x))

    target_slice = slice(
        offset if offset >= 0 else abs(offset) * n,
        min(n * (n - offset), diag.shape[-1]),
        n + 1,
    )
    for index in _helpers.ndindex(*batch_dims):
        diag = at(diag)[(*index, target_slice)].set(x[(*index, slice(None))])
    return xp.reshape(diag, (*batch_dims, n, n))


def one_hot(
    x: Array,
    /,
    num_classes: int,
    *,
    xp: ArrayNamespace,
) -> Array:  # numpydoc ignore=PR01,RT01
    """See docstring in `array_api_extra._creation`."""
    # TODO: Benchmark whether this is faster on the NumPy backend:
    # if is_numpy_array(x):
    #     out = xp.zeros((x.size, num_classes), dtype=dtype)
    #     out[xp.arange(x.size), xp.reshape(x, (-1,))] = 1
    #     return xp.reshape(out, (*x.shape, num_classes))
    range_num_classes = xp.arange(num_classes, dtype=x.dtype, device=_compat.device(x))
    return x[..., xp.newaxis] == range_num_classes
