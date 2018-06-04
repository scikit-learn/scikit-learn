# -*- coding: utf-8 -*-
# pylint: disable-msg=E1101,W0612

from operator import methodcaller
from copy import deepcopy
from distutils.version import LooseVersion

import pytest
import numpy as np
import pandas as pd

from pandas import Series, DataFrame, date_range, MultiIndex

from pandas.compat import range
from pandas.util.testing import (assert_series_equal,
                                 assert_frame_equal,
                                 assert_almost_equal)

import pandas.util.testing as tm
import pandas.util._test_decorators as td
from .test_generic import Generic

try:
    import xarray
    _XARRAY_INSTALLED = True
except ImportError:
    _XARRAY_INSTALLED = False


class TestDataFrame(Generic):
    _typ = DataFrame
    _comparator = lambda self, x, y: assert_frame_equal(x, y)

    def test_rename_mi(self):
        df = DataFrame([
            11, 21, 31
        ], index=MultiIndex.from_tuples([("A", x) for x in ["a", "B", "c"]]))
        df.rename(str.lower)

    def test_set_axis_name(self):
        df = pd.DataFrame([[1, 2], [3, 4]])
        funcs = ['_set_axis_name', 'rename_axis']
        for func in funcs:
            result = methodcaller(func, 'foo')(df)
            assert df.index.name is None
            assert result.index.name == 'foo'

            result = methodcaller(func, 'cols', axis=1)(df)
            assert df.columns.name is None
            assert result.columns.name == 'cols'

    def test_set_axis_name_mi(self):
        df = DataFrame(
            np.empty((3, 3)),
            index=MultiIndex.from_tuples([("A", x) for x in list('aBc')]),
            columns=MultiIndex.from_tuples([('C', x) for x in list('xyz')])
        )

        level_names = ['L1', 'L2']
        funcs = ['_set_axis_name', 'rename_axis']
        for func in funcs:
            result = methodcaller(func, level_names)(df)
            assert result.index.names == level_names
            assert result.columns.names == [None, None]

            result = methodcaller(func, level_names, axis=1)(df)
            assert result.columns.names == ["L1", "L2"]
            assert result.index.names == [None, None]

    def test_nonzero_single_element(self):

        # allow single item via bool method
        df = DataFrame([[True]])
        assert df.bool()

        df = DataFrame([[False]])
        assert not df.bool()

        df = DataFrame([[False, False]])
        pytest.raises(ValueError, lambda: df.bool())
        pytest.raises(ValueError, lambda: bool(df))

    def test_get_numeric_data_preserve_dtype(self):

        # get the numeric data
        o = DataFrame({'A': [1, '2', 3.]})
        result = o._get_numeric_data()
        expected = DataFrame(index=[0, 1, 2], dtype=object)
        self._compare(result, expected)

    def test_metadata_propagation_indiv(self):

        # groupby
        df = DataFrame(
            {'A': ['foo', 'bar', 'foo', 'bar', 'foo', 'bar', 'foo', 'foo'],
             'B': ['one', 'one', 'two', 'three', 'two', 'two', 'one', 'three'],
             'C': np.random.randn(8),
             'D': np.random.randn(8)})
        result = df.groupby('A').sum()
        self.check_metadata(df, result)

        # resample
        df = DataFrame(np.random.randn(1000, 2),
                       index=date_range('20130101', periods=1000, freq='s'))
        result = df.resample('1T')
        self.check_metadata(df, result)

        # merging with override
        # GH 6923
        _metadata = DataFrame._metadata
        _finalize = DataFrame.__finalize__

        np.random.seed(10)
        df1 = DataFrame(np.random.randint(0, 4, (3, 2)), columns=['a', 'b'])
        df2 = DataFrame(np.random.randint(0, 4, (3, 2)), columns=['c', 'd'])
        DataFrame._metadata = ['filename']
        df1.filename = 'fname1.csv'
        df2.filename = 'fname2.csv'

        def finalize(self, other, method=None, **kwargs):

            for name in self._metadata:
                if method == 'merge':
                    left, right = other.left, other.right
                    value = getattr(left, name, '') + '|' + getattr(right,
                                                                    name, '')
                    object.__setattr__(self, name, value)
                else:
                    object.__setattr__(self, name, getattr(other, name, ''))

            return self

        DataFrame.__finalize__ = finalize
        result = df1.merge(df2, left_on=['a'], right_on=['c'], how='inner')
        assert result.filename == 'fname1.csv|fname2.csv'

        # concat
        # GH 6927
        DataFrame._metadata = ['filename']
        df1 = DataFrame(np.random.randint(0, 4, (3, 2)), columns=list('ab'))
        df1.filename = 'foo'

        def finalize(self, other, method=None, **kwargs):
            for name in self._metadata:
                if method == 'concat':
                    value = '+'.join([getattr(
                        o, name) for o in other.objs if getattr(o, name, None)
                    ])
                    object.__setattr__(self, name, value)
                else:
                    object.__setattr__(self, name, getattr(other, name, None))

            return self

        DataFrame.__finalize__ = finalize

        result = pd.concat([df1, df1])
        assert result.filename == 'foo+foo'

        # reset
        DataFrame._metadata = _metadata
        DataFrame.__finalize__ = _finalize

    def test_set_attribute(self):
        # Test for consistent setattr behavior when an attribute and a column
        # have the same name (Issue #8994)
        df = DataFrame({'x': [1, 2, 3]})

        df.y = 2
        df['y'] = [2, 4, 6]
        df.y = 5

        assert df.y == 5
        assert_series_equal(df['y'], Series([2, 4, 6], name='y'))

    @pytest.mark.skipif(not _XARRAY_INSTALLED or _XARRAY_INSTALLED and
                        LooseVersion(xarray.__version__) <
                        LooseVersion('0.10.0'),
                        reason='xarray >= 0.10.0 required')
    @pytest.mark.parametrize(
        "index", ['FloatIndex', 'IntIndex',
                  'StringIndex', 'UnicodeIndex',
                  'DateIndex', 'PeriodIndex',
                  'CategoricalIndex', 'TimedeltaIndex'])
    def test_to_xarray_index_types(self, index):
        from xarray import Dataset

        index = getattr(tm, 'make{}'.format(index))
        df = DataFrame({'a': list('abc'),
                        'b': list(range(1, 4)),
                        'c': np.arange(3, 6).astype('u1'),
                        'd': np.arange(4.0, 7.0, dtype='float64'),
                        'e': [True, False, True],
                        'f': pd.Categorical(list('abc')),
                        'g': pd.date_range('20130101', periods=3),
                        'h': pd.date_range('20130101',
                                           periods=3,
                                           tz='US/Eastern')}
                       )

        df.index = index(3)
        df.index.name = 'foo'
        df.columns.name = 'bar'
        result = df.to_xarray()
        assert result.dims['foo'] == 3
        assert len(result.coords) == 1
        assert len(result.data_vars) == 8
        assert_almost_equal(list(result.coords.keys()), ['foo'])
        assert isinstance(result, Dataset)

        # idempotency
        # categoricals are not preserved
        # datetimes w/tz are not preserved
        # column names are lost
        expected = df.copy()
        expected['f'] = expected['f'].astype(object)
        expected['h'] = expected['h'].astype('datetime64[ns]')
        expected.columns.name = None
        assert_frame_equal(result.to_dataframe(), expected,
                           check_index_type=False, check_categorical=False)

    @td.skip_if_no('xarray', min_version='0.7.0')
    def test_to_xarray(self):
        from xarray import Dataset

        df = DataFrame({'a': list('abc'),
                        'b': list(range(1, 4)),
                        'c': np.arange(3, 6).astype('u1'),
                        'd': np.arange(4.0, 7.0, dtype='float64'),
                        'e': [True, False, True],
                        'f': pd.Categorical(list('abc')),
                        'g': pd.date_range('20130101', periods=3),
                        'h': pd.date_range('20130101',
                                           periods=3,
                                           tz='US/Eastern')}
                       )

        df.index.name = 'foo'
        result = df[0:0].to_xarray()
        assert result.dims['foo'] == 0
        assert isinstance(result, Dataset)

        # available in 0.7.1
        # MultiIndex
        df.index = pd.MultiIndex.from_product([['a'], range(3)],
                                              names=['one', 'two'])
        result = df.to_xarray()
        assert result.dims['one'] == 1
        assert result.dims['two'] == 3
        assert len(result.coords) == 2
        assert len(result.data_vars) == 8
        assert_almost_equal(list(result.coords.keys()), ['one', 'two'])
        assert isinstance(result, Dataset)

        result = result.to_dataframe()
        expected = df.copy()
        expected['f'] = expected['f'].astype(object)
        expected['h'] = expected['h'].astype('datetime64[ns]')
        expected.columns.name = None
        assert_frame_equal(result,
                           expected,
                           check_index_type=False)

    def test_deepcopy_empty(self):
        # This test covers empty frame copying with non-empty column sets
        # as reported in issue GH15370
        empty_frame = DataFrame(data=[], index=[], columns=['A'])
        empty_frame_copy = deepcopy(empty_frame)

        self._compare(empty_frame_copy, empty_frame)
