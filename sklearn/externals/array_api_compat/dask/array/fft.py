from ..._internal import clone_module

__all__ = clone_module("dask.array.fft", globals())

from ...common import _fft
from ..._internal import get_xp

import dask.array as da

fftfreq = get_xp(da)(_fft.fftfreq)
rfftfreq = get_xp(da)(_fft.rfftfreq)

__all__ += ["fftfreq", "rfftfreq"]

def __dir__() -> list[str]:
    return __all__
