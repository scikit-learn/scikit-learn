"""Test Olivetti faces fetcher, if the data is available,
or if specifically requested via environment variable
(e.g. for travis cron job)."""

import os
import pytest
import numpy as np

from sklearn import datasets
from sklearn.utils import Bunch
from sklearn.datasets.tests.test_common import check_return_X_y

from sklearn.utils._testing import assert_array_equal


def _is_olivetti_faces_not_available():
    # Do not download data, unless explicitly requested via environment var
    download_if_missing = False
    if int(os.environ.get('SKLEARN_SKIP_NETWORK_TESTS', 1)) == 0:
        download_if_missing = True
    try:
        datasets.fetch_olivetti_faces(download_if_missing=download_if_missing)
        return False
    except IOError:
        return True


@pytest.mark.skipif(
    _is_olivetti_faces_not_available(),
    reason='Download Olivetti faces dataset to run this test'
)
def test_olivetti_faces():
    data = datasets.fetch_olivetti_faces(shuffle=True, random_state=0)

    assert isinstance(data, Bunch)
    for expected_keys in ('data', 'images', 'target', 'DESCR'):
        assert expected_keys in data.keys()

    assert data.data.shape == (400, 4096)
    assert data.images.shape == (400, 64, 64)
    assert data.target.shape == (400,)
    assert_array_equal(np.unique(np.sort(data.target)), np.arange(40))

    # test the return_X_y option
    check_return_X_y(data, datasets.fetch_olivetti_faces)
