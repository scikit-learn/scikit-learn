from dask.array.fft import * # noqa: F403
# dask.array.fft doesn't have __all__. If it is added, replace this with
#
# from dask.array.fft import __all__ as linalg_all
_n = {}
exec('from dask.array.fft import *', _n)
del _n['__builtins__']
fft_all = list(_n)
del _n

from ...common import _fft
from ..._internal import get_xp

import dask.array as da

fftfreq = get_xp(da)(_fft.fftfreq)
rfftfreq = get_xp(da)(_fft.rfftfreq)

__all__ = [elem for elem in fft_all if elem != "annotations"] + ["fftfreq", "rfftfreq"]

del get_xp
del da
del fft_all
del _fft
