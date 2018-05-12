"""Test the california_housing loader.

Skipped if california_housing is not already downloaded to data_home.
"""

from sklearn.datasets import fetch_california_housing
from sklearn.utils.testing import SkipTest
from sklearn.datasets.tests.test_common import check_return_X_y
from functools import partial


def fetch(*args, **kwargs):
    return fetch_california_housing(*args, download_if_missing=False, **kwargs)


def test_fetch():
    try:
        data = fetch()
    except IOError:
        raise SkipTest("California housing dataset can not be loaded.")
    assert((20640, 8) == data.data.shape)
    assert((20640, ) == data.target.shape)

    # test return_X_y option
    fetch_func = partial(fetch)
    check_return_X_y(data, fetch_func)
