from __future__ import annotations

from collections.abc import Sequence
from typing import Literal

import torch  # noqa: F401
import torch.fft

from ._typing import Array
from .._internal import clone_module

__all__ = clone_module("torch.fft", globals())

# Several torch fft functions do not map axes to dim

def fftn(
    x: Array,
    /,
    *,
    s: Sequence[int] = None,
    axes: Sequence[int] = None,
    norm: Literal["backward", "ortho", "forward"] = "backward",
    **kwargs: object,
) -> Array:
    return torch.fft.fftn(x, s=s, dim=axes, norm=norm, **kwargs)

def ifftn(
    x: Array,
    /,
    *,
    s: Sequence[int] = None,
    axes: Sequence[int] = None,
    norm: Literal["backward", "ortho", "forward"] = "backward",
    **kwargs: object,
) -> Array:
    return torch.fft.ifftn(x, s=s, dim=axes, norm=norm, **kwargs)

def rfftn(
    x: Array,
    /,
    *,
    s: Sequence[int] = None,
    axes: Sequence[int] = None,
    norm: Literal["backward", "ortho", "forward"] = "backward",
    **kwargs: object,
) -> Array:
    return torch.fft.rfftn(x, s=s, dim=axes, norm=norm, **kwargs)

def irfftn(
    x: Array,
    /,
    *,
    s: Sequence[int] = None,
    axes: Sequence[int] = None,
    norm: Literal["backward", "ortho", "forward"] = "backward",
    **kwargs: object,
) -> Array:
    return torch.fft.irfftn(x, s=s, dim=axes, norm=norm, **kwargs)

def fftshift(
    x: Array,
    /,
    *,
    axes: int | Sequence[int] = None,
    **kwargs: object,
) -> Array:
    return torch.fft.fftshift(x, dim=axes, **kwargs)

def ifftshift(
    x: Array,
    /,
    *,
    axes: int | Sequence[int] = None,
    **kwargs: object,
) -> Array:
    return torch.fft.ifftshift(x, dim=axes, **kwargs)


__all__ += ["fftn", "ifftn", "rfftn", "irfftn", "fftshift", "ifftshift"]

def __dir__() -> list[str]:
    return __all__
