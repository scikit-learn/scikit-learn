import numpy as np
from numpy.testing import assert_array_equal

from sklearn.utils._sorting import _py_sort


def test_sort_log2_build():
    """Non-regression test for gh-30554.

    Using log2 and log in sort correctly sorts feature_values, but the tie breaking is
    different which can results in placing samples in a different order.
    """
    rng = np.random.default_rng(75)
    some = rng.normal(loc=0.0, scale=10.0, size=10).astype(np.float32)
    feature_values = np.concatenate([some] * 5)
    samples = np.arange(50, dtype=np.intp)

    _py_sort(feature_values, samples, 50)

    # fmt: off
    # no black reformatting for this specific array
    expected_samples = [
        0, 40, 30, 20, 10, 29, 39, 19, 49,  9, 45, 15, 35,  5, 25, 11, 31,
        41,  1, 21, 22, 12,  2, 42, 32, 23, 13, 43,  3, 33,  6, 36, 46, 16,
        26,  4, 14, 24, 34, 44, 27, 47,  7, 37, 17,  8, 38, 48, 28, 18
    ]
    # fmt: on
    assert_array_equal(samples, expected_samples)
