import pytest
import numpy as np
from pandas import SparseSeries, Series
from pandas.util import testing as tm


pytestmark = pytest.mark.skip("Wrong SparseBlock initialization (GH 17386)")


@pytest.mark.parametrize('data', [
    [1, 1, 2, 2, 3, 3, 4, 4, 0, 0],
    [1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, np.nan, np.nan],
    [
        1.0, 1.0 + 1.0j,
        2.0 + 2.0j, 2.0,
        3.0, 3.0 + 3.0j,
        4.0 + 4.0j, 4.0,
        np.nan, np.nan
    ]
])
@pytest.mark.xfail(reason='Wrong SparseBlock initialization '
                          '(GH 17386)')
def test_where_with_numeric_data(data):
    # GH 17386
    lower_bound = 1.5

    sparse = SparseSeries(data)
    result = sparse.where(sparse > lower_bound)

    dense = Series(data)
    dense_expected = dense.where(dense > lower_bound)
    sparse_expected = SparseSeries(dense_expected)

    tm.assert_series_equal(result, dense_expected)
    tm.assert_sp_series_equal(result, sparse_expected)


@pytest.mark.parametrize('data', [
    [1, 1, 2, 2, 3, 3, 4, 4, 0, 0],
    [1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, np.nan, np.nan],
    [
        1.0, 1.0 + 1.0j,
        2.0 + 2.0j, 2.0,
        3.0, 3.0 + 3.0j,
        4.0 + 4.0j, 4.0,
        np.nan, np.nan
    ]
])
@pytest.mark.parametrize('other', [
    True,
    -100,
    0.1,
    100.0 + 100.0j
])
@pytest.mark.skip(reason='Wrong SparseBlock initialization '
                         '(Segfault) '
                         '(GH 17386)')
def test_where_with_numeric_data_and_other(data, other):
    # GH 17386
    lower_bound = 1.5

    sparse = SparseSeries(data)
    result = sparse.where(sparse > lower_bound, other)

    dense = Series(data)
    dense_expected = dense.where(dense > lower_bound, other)
    sparse_expected = SparseSeries(dense_expected, fill_value=other)

    tm.assert_series_equal(result, dense_expected)
    tm.assert_sp_series_equal(result, sparse_expected)


@pytest.mark.xfail(reason='Wrong SparseBlock initialization '
                          '(GH 17386)')
def test_where_with_bool_data():
    # GH 17386
    data = [False, False, True, True, False, False]
    cond = True

    sparse = SparseSeries(data)
    result = sparse.where(sparse == cond)

    dense = Series(data)
    dense_expected = dense.where(dense == cond)
    sparse_expected = SparseSeries(dense_expected)

    tm.assert_series_equal(result, dense_expected)
    tm.assert_sp_series_equal(result, sparse_expected)


@pytest.mark.parametrize('other', [
    True,
    0,
    0.1,
    100.0 + 100.0j
])
@pytest.mark.skip(reason='Wrong SparseBlock initialization '
                         '(Segfault) '
                         '(GH 17386)')
def test_where_with_bool_data_and_other(other):
    # GH 17386
    data = [False, False, True, True, False, False]
    cond = True

    sparse = SparseSeries(data)
    result = sparse.where(sparse == cond, other)

    dense = Series(data)
    dense_expected = dense.where(dense == cond, other)
    sparse_expected = SparseSeries(dense_expected, fill_value=other)

    tm.assert_series_equal(result, dense_expected)
    tm.assert_sp_series_equal(result, sparse_expected)
