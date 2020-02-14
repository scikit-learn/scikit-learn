"""Test  kddcup99 loader, if the data is available,
or if specifically requested via environment variable
(e.g. for travis cron job).

Only 'percent10' mode is tested, as the full data
is too big to use in unit-testing.
"""

from os import environ
import pytest
from sklearn.datasets import fetch_kddcup99
from sklearn.datasets.tests.test_common import check_return_X_y
from functools import partial


def _fetch_dataset(*args, **kwargs):
    # Do not download data, unless explicitly requested via environment var
    download_if_missing = environ.get('SKLEARN_SKIP_NETWORK_TESTS', '1') == '0'
    try:
        return fetch_kddcup99(*args, download_if_missing=download_if_missing,
                              **kwargs)
    except IOError:
        return None


@pytest.mark.skipif(_fetch_dataset() is None,
                    reason="Download kddcup99 to run this test")
def test_percent10():
    data = _fetch_dataset()

    assert data.data.shape == (494021, 41)
    assert data.target.shape == (494021,)

    data_shuffled = _fetch_dataset(shuffle=True, random_state=0)
    assert data.data.shape == data_shuffled.data.shape
    assert data.target.shape == data_shuffled.target.shape

    data = _fetch_dataset('SA')
    assert data.data.shape == (100655, 41)
    assert data.target.shape == (100655,)

    data = _fetch_dataset('SF')
    assert data.data.shape == (73237, 4)
    assert data.target.shape == (73237,)

    data = _fetch_dataset('http')
    assert data.data.shape == (58725, 3)
    assert data.target.shape == (58725,)

    data = _fetch_dataset('smtp')
    assert data.data.shape == (9571, 3)
    assert data.target.shape == (9571,)

    fetch_func = partial(_fetch_dataset, 'smtp')
    check_return_X_y(data, fetch_func)


@pytest.mark.skipif(_fetch_dataset() is None,
                    reason="Download kddcup99 to run this test")
def test_shuffle():
    dataset = _fetch_dataset(random_state=0, subset='SA', shuffle=True,
                             percent10=True)
    assert(any(dataset.target[-100:] == b'normal.'))
