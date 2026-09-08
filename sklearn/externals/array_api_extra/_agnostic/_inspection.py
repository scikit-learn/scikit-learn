"""Array-agnostic implementations for inspection functions."""

from typing import Literal

from .._lib._typing import ArrayNamespace, Device, DType

__all__ = ["default_dtype"]


def default_dtype(
    xp: ArrayNamespace,
    kind: Literal[
        "real floating", "complex floating", "integral", "indexing"
    ] = "real floating",
    *,
    device: Device | None = None,
) -> DType:
    """
    Return the default dtype for the given namespace and device.

    This is a convenience shorthand for
    ``xp.__array_namespace_info__().default_dtypes(device=device)[kind]``.

    Parameters
    ----------
    xp : array_namespace
        The standard-compatible namespace for which to get the default dtype.
    kind : {'real floating', 'complex floating', 'integral', 'indexing'}, optional
        The kind of dtype to return. Default is 'real floating'.
    device : Device, optional
        The device for which to get the default dtype. Default: current device.

    Returns
    -------
    dtype
        The default dtype for the given namespace, kind, and device.
    """
    dtypes = xp.__array_namespace_info__().default_dtypes(device=device)
    try:
        return dtypes[kind]
    except KeyError as e:
        domain = ("real floating", "complex floating", "integral", "indexing")
        assert set(dtypes) == set(domain), f"Non-compliant namespace: {dtypes}"
        msg = f"Unknown kind '{kind}'. Expected one of {domain}."
        raise ValueError(msg) from e
