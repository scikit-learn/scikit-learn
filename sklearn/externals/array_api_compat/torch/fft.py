from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    import torch
    array = torch.Tensor
    from typing import Union, Sequence, Literal

from torch.fft import * # noqa: F403
import torch.fft

# Several torch fft functions do not map axes to dim

def fftn(
    x: array,
    /,
    *,
    s: Sequence[int] = None,
    axes: Sequence[int] = None,
    norm: Literal["backward", "ortho", "forward"] = "backward",
    **kwargs,
) -> array:
    return torch.fft.fftn(x, s=s, dim=axes, norm=norm, **kwargs)

def ifftn(
    x: array,
    /,
    *,
    s: Sequence[int] = None,
    axes: Sequence[int] = None,
    norm: Literal["backward", "ortho", "forward"] = "backward",
    **kwargs,
) -> array:
    return torch.fft.ifftn(x, s=s, dim=axes, norm=norm, **kwargs)

def rfftn(
    x: array,
    /,
    *,
    s: Sequence[int] = None,
    axes: Sequence[int] = None,
    norm: Literal["backward", "ortho", "forward"] = "backward",
    **kwargs,
) -> array:
    return torch.fft.rfftn(x, s=s, dim=axes, norm=norm, **kwargs)

def irfftn(
    x: array,
    /,
    *,
    s: Sequence[int] = None,
    axes: Sequence[int] = None,
    norm: Literal["backward", "ortho", "forward"] = "backward",
    **kwargs,
) -> array:
    return torch.fft.irfftn(x, s=s, dim=axes, norm=norm, **kwargs)

def fftshift(
    x: array,
    /,
    *,
    axes: Union[int, Sequence[int]] = None,
    **kwargs,
) -> array:
    return torch.fft.fftshift(x, dim=axes, **kwargs)

def ifftshift(
    x: array,
    /,
    *,
    axes: Union[int, Sequence[int]] = None,
    **kwargs,
) -> array:
    return torch.fft.ifftshift(x, dim=axes, **kwargs)


__all__ = torch.fft.__all__ + [
    "fftn",
    "ifftn",
    "rfftn",
    "irfftn",
    "fftshift",
    "ifftshift",
]

_all_ignore = ['torch']
