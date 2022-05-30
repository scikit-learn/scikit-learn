"""Compatibility fixes for older version of python, numpy and scipy

If you add content to this file, please give the version of the package
at which the fix is no longer needed.
"""
# Authors: Emmanuelle Gouillart <emmanuelle.gouillart@normalesup.org>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Fabian Pedregosa <fpedregosa@acm.org>
#          Lars Buitinck
#
# License: BSD 3 clause

from functools import update_wrapper
import functools

import sklearn
import numpy as np
import scipy
import scipy.stats
import threadpoolctl
from .._config import config_context, get_config
from ..externals._packaging.version import parse as parse_version


np_version = parse_version(np.__version__)
sp_version = parse_version(scipy.__version__)


if sp_version >= parse_version("1.4"):
    from scipy.sparse.linalg import lobpcg
else:
    # Backport of lobpcg functionality from scipy 1.4.0, can be removed
    # once support for sp_version < parse_version('1.4') is dropped
    # mypy error: Name 'lobpcg' already defined (possibly by an import)
    from ..externals._lobpcg import lobpcg  # type: ignore  # noqa

try:
    from scipy.optimize._linesearch import line_search_wolfe2, line_search_wolfe1
except ImportError:  # SciPy < 1.8
    from scipy.optimize.linesearch import line_search_wolfe2, line_search_wolfe1  # type: ignore  # noqa


def _object_dtype_isnan(X):
    return X != X


class loguniform(scipy.stats.reciprocal):
    """A class supporting log-uniform random variables.

    Parameters
    ----------
    low : float
        The minimum value
    high : float
        The maximum value

    Methods
    -------
    rvs(self, size=None, random_state=None)
        Generate log-uniform random variables

    The most useful method for Scikit-learn usage is highlighted here.
    For a full list, see
    `scipy.stats.reciprocal
    <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.reciprocal.html>`_.
    This list includes all functions of ``scipy.stats`` continuous
    distributions such as ``pdf``.

    Notes
    -----
    This class generates values between ``low`` and ``high`` or

        low <= loguniform(low, high).rvs() <= high

    The logarithmic probability density function (PDF) is uniform. When
    ``x`` is a uniformly distributed random variable between 0 and 1, ``10**x``
    are random variables that are equally likely to be returned.

    This class is an alias to ``scipy.stats.reciprocal``, which uses the
    reciprocal distribution:
    https://en.wikipedia.org/wiki/Reciprocal_distribution

    Examples
    --------

    >>> from sklearn.utils.fixes import loguniform
    >>> rv = loguniform(1e-3, 1e1)
    >>> rvs = rv.rvs(random_state=42, size=1000)
    >>> rvs.min()  # doctest: +SKIP
    0.0010435856341129003
    >>> rvs.max()  # doctest: +SKIP
    9.97403052786026
    """


# remove when https://github.com/joblib/joblib/issues/1071 is fixed
def delayed(function):
    """Decorator used to capture the arguments of a function."""

    @functools.wraps(function)
    def delayed_function(*args, **kwargs):
        return _FuncWrapper(function), args, kwargs

    return delayed_function


class _FuncWrapper:
    """ "Load the global configuration before calling the function."""

    def __init__(self, function):
        self.function = function
        self.config = get_config()
        update_wrapper(self, self.function)

    def __call__(self, *args, **kwargs):
        with config_context(**self.config):
            return self.function(*args, **kwargs)


# Rename the `method` kwarg to `interpolation` for NumPy < 1.22, because
# `interpolation` kwarg was deprecated in favor of `method` in NumPy >= 1.22.
def _percentile(a, q, *, method="linear", **kwargs):
    return np.percentile(a, q, interpolation=method, **kwargs)


if np_version < parse_version("1.22"):
    percentile = _percentile
else:  # >= 1.22
    from numpy import percentile  # type: ignore  # noqa


# compatibility fix for threadpoolctl >= 3.0.0
# since version 3 it's possible to setup a global threadpool controller to avoid
# looping through all loaded shared libraries each time.
# the global controller is created during the first call to threadpoolctl.
def _get_threadpool_controller():
    if not hasattr(threadpoolctl, "ThreadpoolController"):
        return None

    if not hasattr(sklearn, "_sklearn_threadpool_controller"):
        sklearn._sklearn_threadpool_controller = threadpoolctl.ThreadpoolController()

    return sklearn._sklearn_threadpool_controller


def threadpool_limits(limits=None, user_api=None):
    controller = _get_threadpool_controller()
    if controller is not None:
        return controller.limit(limits=limits, user_api=user_api)
    else:
        return threadpoolctl.threadpool_limits(limits=limits, user_api=user_api)


threadpool_limits.__doc__ = threadpoolctl.threadpool_limits.__doc__


def threadpool_info():
    controller = _get_threadpool_controller()
    if controller is not None:
        return controller.info()
    else:
        return threadpoolctl.threadpool_info()


threadpool_info.__doc__ = threadpoolctl.threadpool_info.__doc__
