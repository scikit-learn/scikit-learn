# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import math
import numbers
from contextlib import suppress
from enum import Enum, unique

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
        The input value to test.

    Returns
    -------
    boolean
        True if `x` is `pandas.NA`, False otherwise.
    """
    with suppress(ImportError):
        from pandas import NA

        return x is NA

    return False


def is_pandas_nat(x):
    """Test if x is pandas.NaT.

    We intentionally do not make this function return `True` for NumPy NaT
    values (e.g., `np.datetime64("NaT")` or `np.timedelta64("NaT")`). Estimators
    that handle NumPy NaT should use `numpy.isnat` instead.

    Parameters
    ----------
    x : any type
        The input value to test.

    Returns
    -------
    boolean
        True if `x` is `pandas.NaT`, False otherwise.
    """
    with suppress(ImportError):
        from pandas import NaT

        return x is NaT
    return False


@unique
class _NAKey(Enum):
    NONE = 0
    NAN = 1
    PD_NA = 2
    PD_NAT = 3
    NAT = 4


def _normalize_na_key(x):
    if x is None:
        return _NAKey.NONE
    if is_scalar_nan(x):
        return _NAKey.NAN
    if is_pandas_na(x):
        return _NAKey.PD_NA
    if is_pandas_nat(x):
        return _NAKey.PD_NAT
    if isinstance(x, (np.datetime64, np.timedelta64)) and np.isnat(x):
        return _NAKey.NAT
    return x
