"""Test the covtype loader, if the data is available,
or if specifically requested via environment variable
(e.g. for travis cron job)."""

import os
from sklearn.datasets import fetch_covtype
from sklearn.utils._testing import SkipTest
from sklearn.datasets.tests.test_common import check_return_X_y
from functools import partial


def fetch(*args, **kwargs):
    # Do not download data, unless explicitly requested via environment var
    download_if_missing = False
    if int(os.environ.get('SKLEARN_SKIP_NETWORK_TESTS', 1)) == 0:
        download_if_missing = True
    return fetch_covtype(*args, download_if_missing=download_if_missing,
                         **kwargs)


def test_fetch():
    try:
        data1 = fetch(shuffle=True, random_state=42)
    except IOError:
        raise SkipTest("Covertype dataset can not be loaded.")

    data2 = fetch(shuffle=True, random_state=37)

    X1, X2 = data1['data'], data2['data']
    assert (581012, 54) == X1.shape
    assert X1.shape == X2.shape

    assert X1.sum() == X2.sum()

    y1, y2 = data1['target'], data2['target']
    assert (X1.shape[0],) == y1.shape
    assert (X1.shape[0],) == y2.shape

    # test return_X_y option
    fetch_func = partial(fetch)
    check_return_X_y(data1, fetch_func)
