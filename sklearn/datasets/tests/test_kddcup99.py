"""Test  kddcup99 loader, if the data is available,
or if specifically requested via environment variable
(e.g. for travis cron job).

Only 'percent10' mode is tested, as the full data
is too big to use in unit-testing.
"""

import os
from sklearn.datasets import fetch_kddcup99
from sklearn.datasets.tests.test_common import check_return_X_y
from sklearn.utils._testing import SkipTest
from functools import partial


def test_percent10():
    # Do not download data, unless explicitly requested via environment var
    download_if_missing = False
    if int(os.environ.get('SKLEARN_SKIP_NETWORK_TESTS', 1)) == 0:
        download_if_missing = True
    try:
        data = fetch_kddcup99(download_if_missing=False)
    except IOError:
        raise SkipTest("kddcup99 dataset can not be loaded.")

    assert data.data.shape == (494021, 41)
    assert data.target.shape == (494021,)

    data_shuffled = fetch_kddcup99(shuffle=True, random_state=0)
    assert data.data.shape == data_shuffled.data.shape
    assert data.target.shape == data_shuffled.target.shape

    data = fetch_kddcup99('SA')
    assert data.data.shape == (100655, 41)
    assert data.target.shape == (100655,)

    data = fetch_kddcup99('SF')
    assert data.data.shape == (73237, 4)
    assert data.target.shape == (73237,)

    data = fetch_kddcup99('http')
    assert data.data.shape == (58725, 3)
    assert data.target.shape == (58725,)

    data = fetch_kddcup99('smtp')
    assert data.data.shape == (9571, 3)
    assert data.target.shape == (9571,)

    fetch_func = partial(fetch_kddcup99, 'smtp')
    check_return_X_y(data, fetch_func)


def test_shuffle():
    try:
        dataset = fetch_kddcup99(random_state=0, subset='SA', shuffle=True,
                                 percent10=True, download_if_missing=False)
    except IOError:
        raise SkipTest("kddcup99 dataset can not be loaded.")

    assert(any(dataset.target[-100:] == b'normal.'))
