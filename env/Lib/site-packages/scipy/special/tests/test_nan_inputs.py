"""Test how the ufuncs in special handle nan inputs.

"""
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_array_equal, assert_
import pytest

import scipy.special as sc
from scipy._lib._numpy_compat import suppress_warnings


KNOWNFAILURES = {}

POSTPROCESSING = {}


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
    with suppress_warnings() as sup:
        # Ignore warnings about unsafe casts from legacy wrappers
        sup.filter(RuntimeWarning,
                   "floating point number truncated to an integer")
        try:
            res = func(*args)
        except TypeError:
            # One of the arguments doesn't take real inputs
            return
    if func in POSTPROCESSING:
        res = POSTPROCESSING[func](*res)

    msg = "got {} instead of nan".format(res)
    assert_array_equal(np.isnan(res), True, err_msg=msg)


def test_legacy_cast():
    with suppress_warnings() as sup:
        sup.filter(RuntimeWarning,
                   "floating point number truncated to an integer")
        res = sc.bdtrc(np.nan, 1, 0.5)
        assert_(np.isnan(res))
