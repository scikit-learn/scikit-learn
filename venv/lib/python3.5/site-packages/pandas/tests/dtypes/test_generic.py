# -*- coding: utf-8 -*-

from warnings import catch_warnings
import numpy as np
import pandas as pd
from pandas.core.dtypes import generic as gt
from pandas.util import testing as tm


class TestABCClasses(object):
    tuples = [[1, 2, 2], ['red', 'blue', 'red']]
    multi_index = pd.MultiIndex.from_arrays(tuples, names=('number', 'color'))
    datetime_index = pd.to_datetime(['2000/1/1', '2010/1/1'])
    timedelta_index = pd.to_timedelta(np.arange(5), unit='s')
    period_index = pd.period_range('2000/1/1', '2010/1/1/', freq='M')
    categorical = pd.Categorical([1, 2, 3], categories=[2, 3, 1])
    categorical_df = pd.DataFrame({"values": [1, 2, 3]}, index=categorical)
    df = pd.DataFrame({'names': ['a', 'b', 'c']}, index=multi_index)
    sparse_series = pd.Series([1, 2, 3]).to_sparse()
    sparse_array = pd.SparseArray(np.random.randn(10))
    sparse_frame = pd.SparseDataFrame({'a': [1, -1, None]})

    def test_abc_types(self):
        assert isinstance(pd.Index(['a', 'b', 'c']), gt.ABCIndex)
        assert isinstance(pd.Int64Index([1, 2, 3]), gt.ABCInt64Index)
        assert isinstance(pd.UInt64Index([1, 2, 3]), gt.ABCUInt64Index)
        assert isinstance(pd.Float64Index([1, 2, 3]), gt.ABCFloat64Index)
        assert isinstance(self.multi_index, gt.ABCMultiIndex)
        assert isinstance(self.datetime_index, gt.ABCDatetimeIndex)
        assert isinstance(self.timedelta_index, gt.ABCTimedeltaIndex)
        assert isinstance(self.period_index, gt.ABCPeriodIndex)
        assert isinstance(self.categorical_df.index, gt.ABCCategoricalIndex)
        assert isinstance(pd.Index(['a', 'b', 'c']), gt.ABCIndexClass)
        assert isinstance(pd.Int64Index([1, 2, 3]), gt.ABCIndexClass)
        assert isinstance(pd.Series([1, 2, 3]), gt.ABCSeries)
        assert isinstance(self.df, gt.ABCDataFrame)
        with catch_warnings(record=True):
            assert isinstance(self.df.to_panel(), gt.ABCPanel)
        assert isinstance(self.sparse_series, gt.ABCSparseSeries)
        assert isinstance(self.sparse_array, gt.ABCSparseArray)
        assert isinstance(self.sparse_frame, gt.ABCSparseDataFrame)
        assert isinstance(self.categorical, gt.ABCCategorical)
        assert isinstance(pd.Period('2012', freq='A-DEC'), gt.ABCPeriod)

        assert isinstance(pd.DateOffset(), gt.ABCDateOffset)
        assert isinstance(pd.Period('2012', freq='A-DEC').freq,
                          gt.ABCDateOffset)
        assert not isinstance(pd.Period('2012', freq='A-DEC'),
                              gt.ABCDateOffset)
        assert isinstance(pd.Interval(0, 1.5), gt.ABCInterval)
        assert not isinstance(pd.Period('2012', freq='A-DEC'), gt.ABCInterval)


def test_setattr_warnings():
    # GH7175 - GOTCHA: You can't use dot notation to add a column...
    d = {'one': pd.Series([1., 2., 3.], index=['a', 'b', 'c']),
         'two': pd.Series([1., 2., 3., 4.], index=['a', 'b', 'c', 'd'])}
    df = pd.DataFrame(d)

    with catch_warnings(record=True) as w:
        #  successfully add new column
        #  this should not raise a warning
        df['three'] = df.two + 1
        assert len(w) == 0
        assert df.three.sum() > df.two.sum()

    with catch_warnings(record=True) as w:
        #  successfully modify column in place
        #  this should not raise a warning
        df.one += 1
        assert len(w) == 0
        assert df.one.iloc[0] == 2

    with catch_warnings(record=True) as w:
        #  successfully add an attribute to a series
        #  this should not raise a warning
        df.two.not_an_index = [1, 2]
        assert len(w) == 0

    with tm.assert_produces_warning(UserWarning):
        #  warn when setting column to nonexistent name
        df.four = df.two + 2
        assert df.four.sum() > df.two.sum()
