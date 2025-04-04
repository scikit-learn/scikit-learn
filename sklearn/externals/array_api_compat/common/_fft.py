from __future__ import annotations

from typing import TYPE_CHECKING, Union, Optional, Literal

if TYPE_CHECKING:
    from ._typing import Device, ndarray, DType
    from collections.abc import Sequence

# Note: NumPy fft functions improperly upcast float32 and complex64 to
# complex128, which is why we require wrapping them all here.

def fft(
    x: ndarray,
    /,
    xp,
    *,
    n: Optional[int] = None,
    axis: int = -1,
    norm: Literal["backward", "ortho", "forward"] = "backward",
) -> ndarray:
    res = xp.fft.fft(x, n=n, axis=axis, norm=norm)
    if x.dtype in [xp.float32, xp.complex64]:
        return res.astype(xp.complex64)
    return res

def ifft(
    x: ndarray,
    /,
    xp,
    *,
    n: Optional[int] = None,
    axis: int = -1,
    norm: Literal["backward", "ortho", "forward"] = "backward",
) -> ndarray:
    res = xp.fft.ifft(x, n=n, axis=axis, norm=norm)
    if x.dtype in [xp.float32, xp.complex64]:
        return res.astype(xp.complex64)
    return res

def fftn(
    x: ndarray,
    /,
    xp,
    *,
    s: Sequence[int] = None,
    axes: Sequence[int] = None,
    norm: Literal["backward", "ortho", "forward"] = "backward",
) -> ndarray:
    res = xp.fft.fftn(x, s=s, axes=axes, norm=norm)
    if x.dtype in [xp.float32, xp.complex64]:
        return res.astype(xp.complex64)
    return res

def ifftn(
    x: ndarray,
    /,
    xp,
    *,
    s: Sequence[int] = None,
    axes: Sequence[int] = None,
    norm: Literal["backward", "ortho", "forward"] = "backward",
) -> ndarray:
    res = xp.fft.ifftn(x, s=s, axes=axes, norm=norm)
    if x.dtype in [xp.float32, xp.complex64]:
        return res.astype(xp.complex64)
    return res

def rfft(
    x: ndarray,
    /,
    xp,
    *,
    n: Optional[int] = None,
    axis: int = -1,
    norm: Literal["backward", "ortho", "forward"] = "backward",
) -> ndarray:
    res = xp.fft.rfft(x, n=n, axis=axis, norm=norm)
    if x.dtype == xp.float32:
        return res.astype(xp.complex64)
    return res

def irfft(
    x: ndarray,
    /,
    xp,
    *,
    n: Optional[int] = None,
    axis: int = -1,
    norm: Literal["backward", "ortho", "forward"] = "backward",
) -> ndarray:
    res = xp.fft.irfft(x, n=n, axis=axis, norm=norm)
    if x.dtype == xp.complex64:
        return res.astype(xp.float32)
    return res

def rfftn(
    x: ndarray,
    /,
    xp,
    *,
    s: Sequence[int] = None,
    axes: Sequence[int] = None,
    norm: Literal["backward", "ortho", "forward"] = "backward",
) -> ndarray:
    res = xp.fft.rfftn(x, s=s, axes=axes, norm=norm)
    if x.dtype == xp.float32:
        return res.astype(xp.complex64)
    return res

def irfftn(
    x: ndarray,
    /,
    xp,
    *,
    s: Sequence[int] = None,
    axes: Sequence[int] = None,
    norm: Literal["backward", "ortho", "forward"] = "backward",
) -> ndarray:
    res = xp.fft.irfftn(x, s=s, axes=axes, norm=norm)
    if x.dtype == xp.complex64:
        return res.astype(xp.float32)
    return res

def hfft(
    x: ndarray,
    /,
    xp,
    *,
    n: Optional[int] = None,
    axis: int = -1,
    norm: Literal["backward", "ortho", "forward"] = "backward",
) -> ndarray:
    res = xp.fft.hfft(x, n=n, axis=axis, norm=norm)
    if x.dtype in [xp.float32, xp.complex64]:
        return res.astype(xp.float32)
    return res

def ihfft(
    x: ndarray,
    /,
    xp,
    *,
    n: Optional[int] = None,
    axis: int = -1,
    norm: Literal["backward", "ortho", "forward"] = "backward",
) -> ndarray:
    res = xp.fft.ihfft(x, n=n, axis=axis, norm=norm)
    if x.dtype in [xp.float32, xp.complex64]:
        return res.astype(xp.complex64)
    return res

def fftfreq(
    n: int,
    /,
    xp,
    *,
    d: float = 1.0,
    dtype: Optional[DType] = None,
    device: Optional[Device] = None
) -> ndarray:
    if device not in ["cpu", None]:
        raise ValueError(f"Unsupported device {device!r}")
    res = xp.fft.fftfreq(n, d=d)
    if dtype is not None:
        return res.astype(dtype)
    return res

def rfftfreq(
    n: int,
    /,
    xp,
    *,
    d: float = 1.0,
    dtype: Optional[DType] = None,
    device: Optional[Device] = None
) -> ndarray:
    if device not in ["cpu", None]:
        raise ValueError(f"Unsupported device {device!r}")
    res = xp.fft.rfftfreq(n, d=d)
    if dtype is not None:
        return res.astype(dtype)
    return res

def fftshift(x: ndarray, /, xp, *, axes: Union[int, Sequence[int]] = None) -> ndarray:
    return xp.fft.fftshift(x, axes=axes)

def ifftshift(x: ndarray, /, xp, *, axes: Union[int, Sequence[int]] = None) -> ndarray:
    return xp.fft.ifftshift(x, axes=axes)

__all__ = [
    "fft",
    "ifft",
    "fftn",
    "ifftn",
    "rfft",
    "irfft",
    "rfftn",
    "irfftn",
    "hfft",
    "ihfft",
    "fftfreq",
    "rfftfreq",
    "fftshift",
    "ifftshift",
]
