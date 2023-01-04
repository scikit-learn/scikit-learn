# Authors: Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Justin Vincent
#          Lars Buitinck
# License: BSD 3 clause

import math
import time

import numpy as np
import pytest
import scipy.stats
from joblib import Parallel

from sklearn.exceptions import ConfigPropagationWarning
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


def _sleep(duration):
    time.sleep(duration)


def test_delayed_warning_config():
    """Check that we raise an informing warning when the user does not pass
    the configuration to be used within the threads executing the tasks.
    """

    n_tasks = 10
    warning_msg = (
        "scikit-learn 1.3 and above will require a config parameter to be passed"
    )
    with pytest.warns(ConfigPropagationWarning, match=warning_msg) as record:
        Parallel(n_jobs=2, pre_dispatch=1)(
            delayed(_sleep)(1e-5) for _ in range(n_tasks)
        )

    assert len(record) == n_tasks
