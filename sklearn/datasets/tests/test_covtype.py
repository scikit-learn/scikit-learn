"""Test the covtype loader, if the data is available,
or if specifically requested via environment variable
(e.g. for travis cron job)."""

from sklearn.datasets.tests.test_common import check_return_X_y
from functools import partial


def test_fetch(fetch_covtype_fxt):
    data1 = fetch_covtype_fxt(shuffle=True, random_state=42)
    data2 = fetch_covtype_fxt(shuffle=True, random_state=37)

    X1, X2 = data1['data'], data2['data']
    assert (581012, 54) == X1.shape
    assert X1.shape == X2.shape

    assert X1.sum() == X2.sum()

    y1, y2 = data1['target'], data2['target']
    assert (X1.shape[0],) == y1.shape
    assert (X1.shape[0],) == y2.shape

    # test return_X_y option
    fetch_func = partial(fetch_covtype_fxt)
    check_return_X_y(data1, fetch_func)
