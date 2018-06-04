# -*- coding: utf-8 -*-
# pylint: disable-msg=E1101,W0612

from operator import methodcaller

import pytest
import numpy as np
import pandas as pd

from distutils.version import LooseVersion
from pandas import Series, date_range, MultiIndex

from pandas.compat import range
from pandas.util.testing import (assert_series_equal,
                                 assert_almost_equal)

import pandas.util.testing as tm
import pandas.util._test_decorators as td
from .test_generic import Generic

try:
    import xarray
    _XARRAY_INSTALLED = True
except ImportError:
    _XARRAY_INSTALLED = False


class TestSeries(Generic):
    _typ = Series
    _comparator = lambda self, x, y: assert_series_equal(x, y)

    def setup_method(self):
        self.ts = tm.makeTimeSeries()  # Was at top level in test_series
        self.ts.name = 'ts'

        self.series = tm.makeStringSeries()
        self.series.name = 'series'

    def test_rename_mi(self):
        s = Series([11, 21, 31],
                   index=MultiIndex.from_tuples(
                       [("A", x) for x in ["a", "B", "c"]]))
        s.rename(str.lower)

    def test_set_axis_name(self):
        s = Series([1, 2, 3], index=['a', 'b', 'c'])
        funcs = ['rename_axis', '_set_axis_name']
        name = 'foo'
        for func in funcs:
            result = methodcaller(func, name)(s)
            assert s.index.name is None
            assert result.index.name == name

    def test_set_axis_name_mi(self):
        s = Series([11, 21, 31], index=MultiIndex.from_tuples(
            [("A", x) for x in ["a", "B", "c"]],
            names=['l1', 'l2'])
        )
        funcs = ['rename_axis', '_set_axis_name']
        for func in funcs:
            result = methodcaller(func, ['L1', 'L2'])(s)
            assert s.index.name is None
            assert s.index.names == ['l1', 'l2']
            assert result.index.name is None
            assert result.index.names, ['L1', 'L2']

    def test_set_axis_name_raises(self):
        s = pd.Series([1])
        with pytest.raises(ValueError):
            s._set_axis_name(name='a', axis=1)

    def test_get_numeric_data_preserve_dtype(self):

        # get the numeric data
        o = Series([1, 2, 3])
        result = o._get_numeric_data()
        self._compare(result, o)

        o = Series([1, '2', 3.])
        result = o._get_numeric_data()
        expected = Series([], dtype=object, index=pd.Index([], dtype=object))
        self._compare(result, expected)

        o = Series([True, False, True])
        result = o._get_numeric_data()
        self._compare(result, o)

        o = Series([True, False, True])
        result = o._get_bool_data()
        self._compare(result, o)

        o = Series(date_range('20130101', periods=3))
        result = o._get_numeric_data()
        expected = Series([], dtype='M8[ns]', index=pd.Index([], dtype=object))
        self._compare(result, expected)

    def test_nonzero_single_element(self):

        # allow single item via bool method
        s = Series([True])
        assert s.bool()

        s = Series([False])
        assert not s.bool()

        # single item nan to raise
        for s in [Series([np.nan]), Series([pd.NaT]), Series([True]),
                  Series([False])]:
            pytest.raises(ValueError, lambda: bool(s))

        for s in [Series([np.nan]), Series([pd.NaT])]:
            pytest.raises(ValueError, lambda: s.bool())

        # multiple bool are still an error
        for s in [Series([True, True]), Series([False, False])]:
            pytest.raises(ValueError, lambda: bool(s))
            pytest.raises(ValueError, lambda: s.bool())

        # single non-bool are an error
        for s in [Series([1]), Series([0]), Series(['a']), Series([0.0])]:
            pytest.raises(ValueError, lambda: bool(s))
            pytest.raises(ValueError, lambda: s.bool())

    def test_metadata_propagation_indiv(self):
        # check that the metadata matches up on the resulting ops

        o = Series(range(3), range(3))
        o.name = 'foo'
        o2 = Series(range(3), range(3))
        o2.name = 'bar'

        result = o.T
        self.check_metadata(o, result)

        # resample
        ts = Series(np.random.rand(1000),
                    index=date_range('20130101', periods=1000, freq='s'),
                    name='foo')
        result = ts.resample('1T').mean()
        self.check_metadata(ts, result)

        result = ts.resample('1T').min()
        self.check_metadata(ts, result)

        result = ts.resample('1T').apply(lambda x: x.sum())
        self.check_metadata(ts, result)

        _metadata = Series._metadata
        _finalize = Series.__finalize__
        Series._metadata = ['name', 'filename']
        o.filename = 'foo'
        o2.filename = 'bar'

        def finalize(self, other, method=None, **kwargs):
            for name in self._metadata:
                if method == 'concat' and name == 'filename':
                    value = '+'.join([getattr(
                        o, name) for o in other.objs if getattr(o, name, None)
                    ])
                    object.__setattr__(self, name, value)
                else:
                    object.__setattr__(self, name, getattr(other, name, None))

            return self

        Series.__finalize__ = finalize

        result = pd.concat([o, o2])
        assert result.filename == 'foo+bar'
        assert result.name is None

        # reset
        Series._metadata = _metadata
        Series.__finalize__ = _finalize

    @pytest.mark.skipif(not _XARRAY_INSTALLED or _XARRAY_INSTALLED and
                        LooseVersion(xarray.__version__) <
                        LooseVersion('0.10.0'),
                        reason='xarray >= 0.10.0 required')
    @pytest.mark.parametrize(
        "index",
        ['FloatIndex', 'IntIndex',
         'StringIndex', 'UnicodeIndex',
         'DateIndex', 'PeriodIndex',
         'TimedeltaIndex', 'CategoricalIndex'])
    def test_to_xarray_index_types(self, index):
        from xarray import DataArray

        index = getattr(tm, 'make{}'.format(index))
        s = Series(range(6), index=index(6))
        s.index.name = 'foo'
        result = s.to_xarray()
        repr(result)
        assert len(result) == 6
        assert len(result.coords) == 1
        assert_almost_equal(list(result.coords.keys()), ['foo'])
        assert isinstance(result, DataArray)

        # idempotency
        assert_series_equal(result.to_series(), s,
                            check_index_type=False,
                            check_categorical=True)

    @td.skip_if_no('xarray', min_version='0.7.0')
    def test_to_xarray(self):
        from xarray import DataArray

        s = Series([])
        s.index.name = 'foo'
        result = s.to_xarray()
        assert len(result) == 0
        assert len(result.coords) == 1
        assert_almost_equal(list(result.coords.keys()), ['foo'])
        assert isinstance(result, DataArray)

        s = Series(range(6))
        s.index.name = 'foo'
        s.index = pd.MultiIndex.from_product([['a', 'b'], range(3)],
                                             names=['one', 'two'])
        result = s.to_xarray()
        assert len(result) == 2
        assert_almost_equal(list(result.coords.keys()), ['one', 'two'])
        assert isinstance(result, DataArray)
        assert_series_equal(result.to_series(), s)

    def test_valid_deprecated(self):
        # GH18800
        with tm.assert_produces_warning(FutureWarning):
            pd.Series([]).valid()
