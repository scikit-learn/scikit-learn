"""Test the landsat loader.

Skipped if landsat is not already downloaded to data_home.
"""

import errno
from sklearn.datasets import fetch_landsat
from sklearn.utils.testing import assert_equal, SkipTest


def fetch(*args, **kwargs):
    return fetch_landsat(*args, download_if_missing=False, **kwargs)


def test_fetch():
    try:
        data1 = fetch(shuffle=True, random_state=42)
    except IOError as e:
        if e.errno == errno.ENOENT:
            raise SkipTest("Covertype dataset can not be loaded.")

    data2 = fetch(shuffle=True, random_state=37)

    X1, X2 = data1['data'], data2['data']
    assert_equal((6435, 36), X1.shape)
    assert_equal(X1.shape, X2.shape)

    assert_equal(X1.sum(), X2.sum())

    y1, y2 = data1['target'], data2['target']
    assert_equal((X1.shape[0],), y1.shape)
    assert_equal((X1.shape[0],), y2.shape)


if __name__ == '__main__':
    test_fetch()