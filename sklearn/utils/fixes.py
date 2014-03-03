"""Compatibility fixes for older version of python, numpy and scipy

If you add content to this file, please give the version of the package
at which the fixe is no longer needed.
"""
# Authors: Emmanuelle Gouillart <emmanuelle.gouillart@normalesup.org>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Fabian Pedregosa <fpedregosa@acm.org>
#          Lars Buitinck
#
# License: BSD 3 clause

import inspect

import numpy as np


np_version = []
for x in np.__version__.split('.'):
    try:
        np_version.append(int(x))
    except ValueError:
        # x may be of the form dev-1ea1592
        np_version.append(x)
np_version = tuple(np_version)


try:
    from scipy.special import expit     # SciPy >= 0.10
except ImportError:
    def expit(x, out=None):
        """Logistic sigmoid function, ``1 / (1 + exp(-x))``.

        See sklearn.utils.extmath.log_logistic for the log of this function.
        """
        if out is None:
            out = np.copy(x)

        # 1 / (1 + exp(-x)) = (1 + tanh(x / 2)) / 2
        # This way of computing the logistic is both fast and stable.
        out *= .5
        np.tanh(out, out)
        out += 1
        out *= .5

        return out


# little danse to see if np.copy has an 'order' keyword argument
if 'order' in inspect.getargspec(np.copy)[0]:
    def safe_copy(X):
        # Copy, but keep the order
        return np.copy(X, order='K')
else:
    # Before an 'order' argument was introduced, numpy wouldn't muck with
    # the ordering
    safe_copy = np.copy

try:
    if (not np.allclose(np.divide(.4, 1), np.divide(.4, 1, dtype=np.float))
            or not np.allclose(np.divide(.4, 1), .4)):
        raise TypeError('Divide not working with dtype: '
                        'https://github.com/numpy/numpy/issues/3484')
    divide = np.divide

except TypeError:
    # Compat for old versions of np.divide that do not provide support for
    # the dtype args
    def divide(x1, x2, out=None, dtype=None):
        out_orig = out
        if out is None:
            out = np.asarray(x1, dtype=dtype)
            if out is x1:
                out = x1.copy()
        else:
            if out is not x1:
                out[:] = x1
        if dtype is not None and out.dtype != dtype:
            out = out.astype(dtype)
        out /= x2
        if out_orig is None and np.isscalar(x1):
            out = np.asscalar(out)
        return out
