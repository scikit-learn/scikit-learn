"""Test how the ufuncs in special handle nan inputs.

"""
import warnings

from collections.abc import Callable

import numpy as np
from numpy.testing import assert_array_equal, assert_
import pytest
import scipy.special as sc


KNOWNFAILURES: dict[str, Callable] = {}

POSTPROCESSING: dict[str, Callable] = {}


def _get_ufuncs():
    ufuncs = []
    ufunc_names = []
    for name in sorted(sc.__dict__):
        obj = sc.__dict__[name]
        if not isinstance(obj, np.ufunc):
            continue
        msg = KNOWNFAILURES.get(obj)
        if msg is None:
            ufuncs.append(obj)
            ufunc_names.append(name)
        else:
            fail = pytest.mark.xfail(run=False, reason=msg)
            ufuncs.append(pytest.param(obj, marks=fail))
            ufunc_names.append(name)
    return ufuncs, ufunc_names


UFUNCS, UFUNC_NAMES = _get_ufuncs()


@pytest.mark.parametrize("func", UFUNCS, ids=UFUNC_NAMES)
def test_nan_inputs(func):
    args = (np.nan,)*func.nin
    with warnings.catch_warnings():
        # Ignore warnings about unsafe casts from legacy wrappers
        warnings.filterwarnings(
            "ignore",
            "floating point number truncated to an integer",
            RuntimeWarning
        )
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", DeprecationWarning)
                res = func(*args)
        except TypeError:
            # One of the arguments doesn't take real inputs
            return
    if func in POSTPROCESSING:
        res = POSTPROCESSING[func](*res)

    msg = f"got {res} instead of nan"
    assert_array_equal(np.isnan(res), True, err_msg=msg)


def test_legacy_cast():
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            "floating point number truncated to an integer",
            RuntimeWarning
        )
        res = sc.bdtrc(np.nan, 1, 0.5)
        assert_(np.isnan(res))
