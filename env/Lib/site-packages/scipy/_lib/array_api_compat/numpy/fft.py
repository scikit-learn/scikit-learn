from numpy.fft import * # noqa: F403
from numpy.fft import __all__ as fft_all

from ..common import _fft
from .._internal import get_xp

import numpy as np

fft = get_xp(np)(_fft.fft)
ifft = get_xp(np)(_fft.ifft)
fftn = get_xp(np)(_fft.fftn)
ifftn = get_xp(np)(_fft.ifftn)
rfft = get_xp(np)(_fft.rfft)
irfft = get_xp(np)(_fft.irfft)
rfftn = get_xp(np)(_fft.rfftn)
irfftn = get_xp(np)(_fft.irfftn)
hfft = get_xp(np)(_fft.hfft)
ihfft = get_xp(np)(_fft.ihfft)
fftfreq = get_xp(np)(_fft.fftfreq)
rfftfreq = get_xp(np)(_fft.rfftfreq)
fftshift = get_xp(np)(_fft.fftshift)
ifftshift = get_xp(np)(_fft.ifftshift)

__all__ = fft_all + _fft.__all__

del get_xp
del np
del fft_all
del _fft
