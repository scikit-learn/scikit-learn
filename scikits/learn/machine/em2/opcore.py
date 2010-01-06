import numpy as np

from scikits.learn.machine.em2._lk import logresp as _logresp

def logresp(x, w, mu, va, y=None):
    k = w.size
    n = x.shape[0]

    if not y:
        y = np.empty((n, k), x.dtype)

    _logresp(x, w, mu, va, y)
    return y
