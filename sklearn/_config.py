"""Global configuration state and functions for management
"""
import os
from contextlib import contextmanager as contextmanager

_ASSUME_FINITE = bool(os.environ.get('SKLEARN_ASSUME_FINITE', False))
_WORKING_MEMORY = int(os.environ.get('SKLEARN_WORKING_MEMORY', 64))


def get_config():
    """Retrieve current values for configuration set by :func:`set_config`

    Returns
    -------
    config : dict
        Keys are parameter names that can be passed to :func:`set_config`.
    """
    return {'assume_finite': _ASSUME_FINITE,
            'working_memory': _WORKING_MEMORY}


def set_config(assume_finite=None, working_memory=None):
    """Set global scikit-learn configuration

    Parameters
    ----------
    assume_finite : bool, optional
        If True, validation for finiteness will be skipped,
        saving time, but leading to potential crashes. If
        False, validation for finiteness will be performed,
        avoiding error.  Global default: False.

    working_memory : int, optional
        If set, scikit-learn will attempt to limit the size of temporary arrays
        to this number of MiB (per job when parallelised), often saving both
        computation time and memory on expensive operations that can be
        performed in chunks. Global default: 64.
    """
    global _ASSUME_FINITE, _WORKING_MEMORY
    if assume_finite is not None:
        _ASSUME_FINITE = assume_finite
    if working_memory is not None:
        _WORKING_MEMORY = working_memory


@contextmanager
def config_context(**new_config):
    """Context manager for global scikit-learn configuration

    Parameters
    ----------
    assume_finite : bool, optional
        If True, validation for finiteness will be skipped,
        saving time, but leading to potential crashes. If
        False, validation for finiteness will be performed,
        avoiding error.  Global default: False.

    working_memory : int, optional
        If set, scikit-learn will attempt to limit the size of temporary arrays
        to this number of MiB (per job when parallelised), often saving both
        computation time and memory on expensive operations that can be
        performed in chunks. Global default: 64.

    Notes
    -----
    All settings, not just those presently modified, will be returned to
    their previous values when the context manager is exited. This is not
    thread-safe.

    Examples
    --------
    >>> import sklearn
    >>> from sklearn.utils.validation import assert_all_finite
    >>> with sklearn.config_context(assume_finite=True):
    ...     assert_all_finite([float('nan')])
    >>> with sklearn.config_context(assume_finite=True):
    ...     with sklearn.config_context(assume_finite=False):
    ...         assert_all_finite([float('nan')])
    ... # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    ValueError: Input contains NaN, ...
    """
    old_config = get_config().copy()
    set_config(**new_config)

    try:
        yield
    finally:
        set_config(**old_config)
