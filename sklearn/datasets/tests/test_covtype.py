"""Test the covtype loader, if the data is available,
or if specifically requested via environment variable
(e.g. for travis cron job)."""
import pytest
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

def test_fetch_asframe(fetch_covtype_fxt):
    pd = pytest.importorskip('pandas')
    bunch = fetch_covtype_fxt(as_frame=True)
    frame = bunch.frame
    assert hasattr(bunch, frame) is True
    assert frame.shape == (581012, 55)
    assert isinstance(bunch.data, pd.DataFrame)
    assert isinstance(bunch.target, pd.Series)

def test_pandas_dependency_message(fetch_covtype_fxt,
                                   hide_available_pandas):
    expected_msg = ('fetch_covtype_fxt with as_frame=True'
                    ' requires pandas')
    with pytest.raises(ImportError, match=expected_msg):
        fetch_covtype_fxt(as_frame=True)
