# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import math
import numbers
from contextlib import suppress

import numpy as np


def is_scalar_nan(x):
    """Test if x is NaN.

    This function is meant to overcome the issue that np.isnan does not allow
    non-numerical types as input, and that np.nan is not float('nan').

    Parameters
    ----------
    x : any type
        Any scalar value.

    Returns
    -------
    bool
        Returns true if x is NaN, and false otherwise.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.utils._missing import is_scalar_nan
    >>> is_scalar_nan(np.nan)
    True
    >>> is_scalar_nan(float("nan"))
    True
    >>> is_scalar_nan(None)
    False
    >>> is_scalar_nan("")
    False
    >>> is_scalar_nan([np.nan])
    False
    """
    return (
        not isinstance(x, numbers.Integral)
        and isinstance(x, numbers.Real)
        and math.isnan(x)
    )


def is_pandas_na(x):
    """Test if x is pandas.NA.

    We intentionally do not use this function to return `True` for `pd.NA` in
    `is_scalar_nan`, because estimators that support `pd.NA` are the exception
    rather than the rule at the moment. When `pd.NA` is more universally
    supported, we may reconsider this decision.

    Parameters
    ----------
    x : any type

    Returns
    -------
    boolean
    """
    with suppress(ImportError):
        from pandas import NA

        return x is NA

    return False


def array_like_has_any_na(x):
    """Test if x is list, tuple or set and have pd.NA o NaN

    Parameters
    ----------
    x : any type

    Returns
    -------
    boolean
    """
    if isinstance(x, (list, set, tuple)):
        return any(is_scalar_nan(v) or is_pandas_na(v) for v in x)
    return False


def has_only_number_or_none(x):
    """Test if x is number or none, or if x is list, tuple or set
    that has only numbers or None

    Parameters
    ----------
    x : any type

    Returns
    -------
    boolean
    """
    if isinstance(x, (list, set, tuple)):
        return any(isinstance(v, numbers.Real) or v is None for v in x)
    return isinstance(x, numbers.Real) or x is None


try:
    import pandas

    is_pandas_na_only = np.frompyfunc(lambda a: a is pandas.NA or a is pandas.NaT, 1, 1)
except ImportError:
    is_pandas_na_only = lambda a: False
