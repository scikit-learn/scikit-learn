# Authors: Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Justin Vincent
#          Lars Buitinck
# License: BSD 3 clause

import math
import threading

import numpy as np
import pytest
import scipy.stats

from joblib import Parallel

import sklearn
from sklearn.utils._testing import assert_array_equal

from sklearn.utils.fixes import _object_dtype_isnan, delayed, loguniform


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


def test_delayed_fetching_right_config():
    """Check that `delayed` function fetches the right config associated to
    the main thread.

    Non-regression test for:
    https://github.com/scikit-learn/scikit-learn/issues/25239
    """

    def get_working_memory():
        return sklearn.get_config()["working_memory"]

    n_iter = 10

    # by default, we register the main thread and we should retrieve the
    # parameters defined within the context manager
    with sklearn.config_context(working_memory=123):
        results = Parallel(n_jobs=2, pre_dispatch=4)(
            delayed(get_working_memory)() for _ in range(n_iter)
        )

    assert results == [123] * n_iter

    # simulate that we refer to another thread
    local_thread = threading.Thread(target=sklearn.get_config)
    local_thread.start()
    local_thread.join()
    with sklearn.config_context(working_memory=123):
        results = Parallel(n_jobs=2, pre_dispatch=4)(
            delayed(get_working_memory, thread=local_thread)() for _ in range(n_iter)
        )

    assert results == [get_working_memory()] * n_iter
