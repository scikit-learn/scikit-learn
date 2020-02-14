"""Test the covtype loader, if the data is available,
or if specifically requested via environment variable
(e.g. for travis cron job)."""

from os import environ
import pytest
from sklearn.datasets import fetch_covtype
from sklearn.datasets.tests.test_common import check_return_X_y
from functools import partial


def fetch(*args, **kwargs):
    # Do not download data, unless explicitly requested via environment var
    download_if_missing = environ.get('SKLEARN_SKIP_NETWORK_TESTS', '1') == '0'
    try:
        return fetch_covtype(*args, download_if_missing=download_if_missing,
                             **kwargs)
    except IOError:
        return None


@pytest.mark.skipif(fetch() is None,
                    reason="Download covtype to run this test")
def test_fetch():
    data1 = fetch(shuffle=True, random_state=42)
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
