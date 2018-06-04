import pytest
import numpy as np
from pandas import SparseDataFrame, DataFrame
from pandas.util import testing as tm


pytestmark = pytest.mark.skip("Wrong SparseBlock initialization (GH 17386)")


@pytest.mark.parametrize('data', [
    [[1, 1], [2, 2], [3, 3], [4, 4], [0, 0]],
    [[1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [4.0, 4.0], [np.nan, np.nan]],
    [
        [1.0, 1.0 + 1.0j],
        [2.0 + 2.0j, 2.0],
        [3.0, 3.0 + 3.0j],
        [4.0 + 4.0j, 4.0],
        [np.nan, np.nan]
    ]
])
@pytest.mark.xfail(reason='Wrong SparseBlock initialization '
                          '(GH 17386)')
def test_where_with_numeric_data(data):
    # GH 17386
    lower_bound = 1.5

    sparse = SparseDataFrame(data)
    result = sparse.where(sparse > lower_bound)

    dense = DataFrame(data)
    dense_expected = dense.where(dense > lower_bound)
    sparse_expected = SparseDataFrame(dense_expected)

    tm.assert_frame_equal(result, dense_expected)
    tm.assert_sp_frame_equal(result, sparse_expected)


@pytest.mark.parametrize('data', [
    [[1, 1], [2, 2], [3, 3], [4, 4], [0, 0]],
    [[1.0, 1.0], [2.0, 2.0], [3.0, 3.0], [4.0, 4.0], [np.nan, np.nan]],
    [
        [1.0, 1.0 + 1.0j],
        [2.0 + 2.0j, 2.0],
        [3.0, 3.0 + 3.0j],
        [4.0 + 4.0j, 4.0],
        [np.nan, np.nan]
    ]
])
@pytest.mark.parametrize('other', [
    True,
    -100,
    0.1,
    100.0 + 100.0j
])
@pytest.mark.xfail(reason='Wrong SparseBlock initialization '
                          '(GH 17386)')
def test_where_with_numeric_data_and_other(data, other):
    # GH 17386
    lower_bound = 1.5

    sparse = SparseDataFrame(data)
    result = sparse.where(sparse > lower_bound, other)

    dense = DataFrame(data)
    dense_expected = dense.where(dense > lower_bound, other)
    sparse_expected = SparseDataFrame(dense_expected,
                                      default_fill_value=other)

    tm.assert_frame_equal(result, dense_expected)
    tm.assert_sp_frame_equal(result, sparse_expected)


@pytest.mark.xfail(reason='Wrong SparseBlock initialization '
                          '(GH 17386)')
def test_where_with_bool_data():
    # GH 17386
    data = [[False, False], [True, True], [False, False]]
    cond = True

    sparse = SparseDataFrame(data)
    result = sparse.where(sparse == cond)

    dense = DataFrame(data)
    dense_expected = dense.where(dense == cond)
    sparse_expected = SparseDataFrame(dense_expected)

    tm.assert_frame_equal(result, dense_expected)
    tm.assert_sp_frame_equal(result, sparse_expected)


@pytest.mark.parametrize('other', [
    True,
    0,
    0.1,
    100.0 + 100.0j
])
@pytest.mark.xfail(reason='Wrong SparseBlock initialization '
                          '(GH 17386)')
def test_where_with_bool_data_and_other(other):
    # GH 17386
    data = [[False, False], [True, True], [False, False]]
    cond = True

    sparse = SparseDataFrame(data)
    result = sparse.where(sparse == cond, other)

    dense = DataFrame(data)
    dense_expected = dense.where(dense == cond, other)
    sparse_expected = SparseDataFrame(dense_expected,
                                      default_fill_value=other)

    tm.assert_frame_equal(result, dense_expected)
    tm.assert_sp_frame_equal(result, sparse_expected)
