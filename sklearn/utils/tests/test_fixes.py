# Authors: Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Justin Vincent
#          Lars Buitinck
# License: BSD 3 clause

import math

import numpy as np
import pytest
import scipy.stats

from sklearn.utils._testing import assert_array_equal

from sklearn.utils.fixes import _object_dtype_isnan
from sklearn.utils.fixes import loguniform
from sklearn.utils.fixes import mode


@pytest.mark.parametrize("dtype, val", ([object, 1], [object, "a"], [float, 1]))
def test_object_dtype_isnan(dtype, val):
    X = np.array([[val, np.nan], [np.nan, val]], dtype=dtype)

    expected_mask = np.array([[False, True], [True, False]])

    mask = _object_dtype_isnan(X)

    assert_array_equal(mask, expected_mask)


@pytest.mark.parametrize("low,high,base", [(-1, 0, 10), (0, 2, np.exp(1)), (-1, 1, 2)])
def test_loguniform(low, high, base):
    rv = loguniform(base**low, base**high)
    assert isinstance(rv, scipy.stats._distn_infrastructure.rv_frozen)
    rvs = rv.rvs(size=2000, random_state=0)

    # Test the basics; right bounds, right size
    assert (base**low <= rvs).all() and (rvs <= base**high).all()
    assert len(rvs) == 2000

    # Test that it's actually (fairly) uniform
    log_rvs = np.array([math.log(x, base) for x in rvs])
    counts, _ = np.histogram(log_rvs)
    assert counts.mean() == 200
    assert np.abs(counts - counts.mean()).max() <= 40

    # Test that random_state works
    assert loguniform(base**low, base**high).rvs(random_state=0) == loguniform(
        base**low, base**high
    ).rvs(random_state=0)


def test_mode():
    """Check that `mode` returns the expected type of data and array dimensions."""
    a = np.array([[6, 8, 3, 0], [3, 2, 1, 7], [8, 1, 8, 4], [5, 3, 0, 5], [4, 7, 5, 9]])
    mode_result = mode(a)
    assert_array_equal(mode_result.mode, np.array([3, 1, 0, 0]))
    assert_array_equal(mode_result.count, np.array([1, 1, 1, 1]))

    mode_result = mode(a, axis=None)
    assert mode_result.mode == 3
    assert mode_result.count == 3
