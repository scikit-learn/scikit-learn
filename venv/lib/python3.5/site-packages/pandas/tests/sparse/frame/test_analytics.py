import pytest
import numpy as np
from pandas import SparseDataFrame, DataFrame, SparseSeries
from pandas.util import testing as tm


@pytest.mark.xfail(reason='Wrong SparseBlock initialization '
                   '(GH 17386)')
def test_quantile():
    # GH 17386
    data = [[1, 1], [2, 10], [3, 100], [np.nan, np.nan]]
    q = 0.1

    sparse_df = SparseDataFrame(data)
    result = sparse_df.quantile(q)

    dense_df = DataFrame(data)
    dense_expected = dense_df.quantile(q)
    sparse_expected = SparseSeries(dense_expected)

    tm.assert_series_equal(result, dense_expected)
    tm.assert_sp_series_equal(result, sparse_expected)


@pytest.mark.xfail(reason='Wrong SparseBlock initialization '
                   '(GH 17386)')
def test_quantile_multi():
    # GH 17386
    data = [[1, 1], [2, 10], [3, 100], [np.nan, np.nan]]
    q = [0.1, 0.5]

    sparse_df = SparseDataFrame(data)
    result = sparse_df.quantile(q)

    dense_df = DataFrame(data)
    dense_expected = dense_df.quantile(q)
    sparse_expected = SparseDataFrame(dense_expected)

    tm.assert_frame_equal(result, dense_expected)
    tm.assert_sp_frame_equal(result, sparse_expected)
