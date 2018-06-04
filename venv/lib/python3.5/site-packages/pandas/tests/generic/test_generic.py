# -*- coding: utf-8 -*-
# pylint: disable-msg=E1101,W0612

from copy import copy, deepcopy
from warnings import catch_warnings

import pytest
import numpy as np
import pandas as pd

from pandas.core.dtypes.common import is_scalar
from pandas import (Series, DataFrame, Panel,
                    date_range, MultiIndex)

import pandas.io.formats.printing as printing

from pandas.compat import range, zip, PY3
from pandas.util.testing import (assert_raises_regex,
                                 assert_series_equal,
                                 assert_panel_equal,
                                 assert_frame_equal)

import pandas.util.testing as tm


# ----------------------------------------------------------------------
# Generic types test cases

class Generic(object):

    @property
    def _ndim(self):
        return self._typ._AXIS_LEN

    def _axes(self):
        """ return the axes for my object typ """
        return self._typ._AXIS_ORDERS

    def _construct(self, shape, value=None, dtype=None, **kwargs):
        """ construct an object for the given shape
            if value is specified use that if its a scalar
            if value is an array, repeat it as needed """

        if isinstance(shape, int):
            shape = tuple([shape] * self._ndim)
        if value is not None:
            if is_scalar(value):
                if value == 'empty':
                    arr = None

                    # remove the info axis
                    kwargs.pop(self._typ._info_axis_name, None)
                else:
                    arr = np.empty(shape, dtype=dtype)
                    arr.fill(value)
            else:
                fshape = np.prod(shape)
                arr = value.ravel()
                new_shape = fshape / arr.shape[0]
                if fshape % arr.shape[0] != 0:
                    raise Exception("invalid value passed in _construct")

                arr = np.repeat(arr, new_shape).reshape(shape)
        else:
            arr = np.random.randn(*shape)
        return self._typ(arr, dtype=dtype, **kwargs)

    def _compare(self, result, expected):
        self._comparator(result, expected)

    def test_rename(self):

        # single axis
        idx = list('ABCD')
        # relabeling values passed into self.rename
        args = [
            str.lower,
            {x: x.lower() for x in idx},
            Series({x: x.lower() for x in idx}),
        ]

        for axis in self._axes():
            kwargs = {axis: idx}
            obj = self._construct(4, **kwargs)

            for arg in args:
                # rename a single axis
                result = obj.rename(**{axis: arg})
                expected = obj.copy()
                setattr(expected, axis, list('abcd'))
                self._compare(result, expected)

        # multiple axes at once

    def test_get_numeric_data(self):

        n = 4
        kwargs = {}
        for i in range(self._ndim):
            kwargs[self._typ._AXIS_NAMES[i]] = list(range(n))

        # get the numeric data
        o = self._construct(n, **kwargs)
        result = o._get_numeric_data()
        self._compare(result, o)

        # non-inclusion
        result = o._get_bool_data()
        expected = self._construct(n, value='empty', **kwargs)
        self._compare(result, expected)

        # get the bool data
        arr = np.array([True, True, False, True])
        o = self._construct(n, value=arr, **kwargs)
        result = o._get_numeric_data()
        self._compare(result, o)

        # _get_numeric_data is includes _get_bool_data, so can't test for
        # non-inclusion

    def test_get_default(self):

        # GH 7725
        d0 = "a", "b", "c", "d"
        d1 = np.arange(4, dtype='int64')
        others = "e", 10

        for data, index in ((d0, d1), (d1, d0)):
            s = Series(data, index=index)
            for i, d in zip(index, data):
                assert s.get(i) == d
                assert s.get(i, d) == d
                assert s.get(i, "z") == d
                for other in others:
                    assert s.get(other, "z") == "z"
                    assert s.get(other, other) == other

    def test_nonzero(self):

        # GH 4633
        # look at the boolean/nonzero behavior for objects
        obj = self._construct(shape=4)
        pytest.raises(ValueError, lambda: bool(obj == 0))
        pytest.raises(ValueError, lambda: bool(obj == 1))
        pytest.raises(ValueError, lambda: bool(obj))

        obj = self._construct(shape=4, value=1)
        pytest.raises(ValueError, lambda: bool(obj == 0))
        pytest.raises(ValueError, lambda: bool(obj == 1))
        pytest.raises(ValueError, lambda: bool(obj))

        obj = self._construct(shape=4, value=np.nan)
        pytest.raises(ValueError, lambda: bool(obj == 0))
        pytest.raises(ValueError, lambda: bool(obj == 1))
        pytest.raises(ValueError, lambda: bool(obj))

        # empty
        obj = self._construct(shape=0)
        pytest.raises(ValueError, lambda: bool(obj))

        # invalid behaviors

        obj1 = self._construct(shape=4, value=1)
        obj2 = self._construct(shape=4, value=1)

        def f():
            if obj1:
                printing.pprint_thing("this works and shouldn't")

        pytest.raises(ValueError, f)
        pytest.raises(ValueError, lambda: obj1 and obj2)
        pytest.raises(ValueError, lambda: obj1 or obj2)
        pytest.raises(ValueError, lambda: not obj1)

    def test_downcast(self):
        # test close downcasting

        o = self._construct(shape=4, value=9, dtype=np.int64)
        result = o.copy()
        result._data = o._data.downcast(dtypes='infer')
        self._compare(result, o)

        o = self._construct(shape=4, value=9.)
        expected = o.astype(np.int64)
        result = o.copy()
        result._data = o._data.downcast(dtypes='infer')
        self._compare(result, expected)

        o = self._construct(shape=4, value=9.5)
        result = o.copy()
        result._data = o._data.downcast(dtypes='infer')
        self._compare(result, o)

        # are close
        o = self._construct(shape=4, value=9.000000000005)
        result = o.copy()
        result._data = o._data.downcast(dtypes='infer')
        expected = o.astype(np.int64)
        self._compare(result, expected)

    def test_constructor_compound_dtypes(self):
        # GH 5191
        # compound dtypes should raise not-implementederror

        def f(dtype):
            return self._construct(shape=3, dtype=dtype)

        pytest.raises(NotImplementedError, f, [("A", "datetime64[h]"),
                                               ("B", "str"),
                                               ("C", "int32")])

        # these work (though results may be unexpected)
        f('int64')
        f('float64')
        f('M8[ns]')

    def check_metadata(self, x, y=None):
        for m in x._metadata:
            v = getattr(x, m, None)
            if y is None:
                assert v is None
            else:
                assert v == getattr(y, m, None)

    def test_metadata_propagation(self):
        # check that the metadata matches up on the resulting ops

        o = self._construct(shape=3)
        o.name = 'foo'
        o2 = self._construct(shape=3)
        o2.name = 'bar'

        # TODO
        # Once panel can do non-trivial combine operations
        # (currently there is an a raise in the Panel arith_ops to prevent
        # this, though it actually does work)
        # can remove all of these try: except: blocks on the actual operations

        # ----------
        # preserving
        # ----------

        # simple ops with scalars
        for op in ['__add__', '__sub__', '__truediv__', '__mul__']:
            result = getattr(o, op)(1)
            self.check_metadata(o, result)

        # ops with like
        for op in ['__add__', '__sub__', '__truediv__', '__mul__']:
            try:
                result = getattr(o, op)(o)
                self.check_metadata(o, result)
            except (ValueError, AttributeError):
                pass

        # simple boolean
        for op in ['__eq__', '__le__', '__ge__']:
            v1 = getattr(o, op)(o)
            self.check_metadata(o, v1)

            try:
                self.check_metadata(o, v1 & v1)
            except (ValueError):
                pass

            try:
                self.check_metadata(o, v1 | v1)
            except (ValueError):
                pass

        # combine_first
        try:
            result = o.combine_first(o2)
            self.check_metadata(o, result)
        except (AttributeError):
            pass

        # ---------------------------
        # non-preserving (by default)
        # ---------------------------

        # add non-like
        try:
            result = o + o2
            self.check_metadata(result)
        except (ValueError, AttributeError):
            pass

        # simple boolean
        for op in ['__eq__', '__le__', '__ge__']:

            # this is a name matching op
            v1 = getattr(o, op)(o)

            v2 = getattr(o, op)(o2)
            self.check_metadata(v2)

            try:
                self.check_metadata(v1 & v2)
            except (ValueError):
                pass

            try:
                self.check_metadata(v1 | v2)
            except (ValueError):
                pass

    def test_head_tail(self):
        # GH5370

        o = self._construct(shape=10)

        # check all index types
        for index in [tm.makeFloatIndex, tm.makeIntIndex, tm.makeStringIndex,
                      tm.makeUnicodeIndex, tm.makeDateIndex,
                      tm.makePeriodIndex]:
            axis = o._get_axis_name(0)
            setattr(o, axis, index(len(getattr(o, axis))))

            # Panel + dims
            try:
                o.head()
            except (NotImplementedError):
                pytest.skip('not implemented on {0}'.format(
                    o.__class__.__name__))

            self._compare(o.head(), o.iloc[:5])
            self._compare(o.tail(), o.iloc[-5:])

            # 0-len
            self._compare(o.head(0), o.iloc[0:0])
            self._compare(o.tail(0), o.iloc[0:0])

            # bounded
            self._compare(o.head(len(o) + 1), o)
            self._compare(o.tail(len(o) + 1), o)

            # neg index
            self._compare(o.head(-3), o.head(7))
            self._compare(o.tail(-3), o.tail(7))

    def test_sample(self):
        # Fixes issue: 2419

        o = self._construct(shape=10)

        ###
        # Check behavior of random_state argument
        ###

        # Check for stability when receives seed or random state -- run 10
        # times.
        for test in range(10):
            seed = np.random.randint(0, 100)
            self._compare(
                o.sample(n=4, random_state=seed), o.sample(n=4,
                                                           random_state=seed))
            self._compare(
                o.sample(frac=0.7, random_state=seed), o.sample(
                    frac=0.7, random_state=seed))

            self._compare(
                o.sample(n=4, random_state=np.random.RandomState(test)),
                o.sample(n=4, random_state=np.random.RandomState(test)))

            self._compare(
                o.sample(frac=0.7, random_state=np.random.RandomState(test)),
                o.sample(frac=0.7, random_state=np.random.RandomState(test)))

            os1, os2 = [], []
            for _ in range(2):
                np.random.seed(test)
                os1.append(o.sample(n=4))
                os2.append(o.sample(frac=0.7))
            self._compare(*os1)
            self._compare(*os2)

        # Check for error when random_state argument invalid.
        with pytest.raises(ValueError):
            o.sample(random_state='astring!')

        ###
        # Check behavior of `frac` and `N`
        ###

        # Giving both frac and N throws error
        with pytest.raises(ValueError):
            o.sample(n=3, frac=0.3)

        # Check that raises right error for negative lengths
        with pytest.raises(ValueError):
            o.sample(n=-3)
        with pytest.raises(ValueError):
            o.sample(frac=-0.3)

        # Make sure float values of `n` give error
        with pytest.raises(ValueError):
            o.sample(n=3.2)

        # Check lengths are right
        assert len(o.sample(n=4) == 4)
        assert len(o.sample(frac=0.34) == 3)
        assert len(o.sample(frac=0.36) == 4)

        ###
        # Check weights
        ###

        # Weight length must be right
        with pytest.raises(ValueError):
            o.sample(n=3, weights=[0, 1])

        with pytest.raises(ValueError):
            bad_weights = [0.5] * 11
            o.sample(n=3, weights=bad_weights)

        with pytest.raises(ValueError):
            bad_weight_series = Series([0, 0, 0.2])
            o.sample(n=4, weights=bad_weight_series)

        # Check won't accept negative weights
        with pytest.raises(ValueError):
            bad_weights = [-0.1] * 10
            o.sample(n=3, weights=bad_weights)

        # Check inf and -inf throw errors:
        with pytest.raises(ValueError):
            weights_with_inf = [0.1] * 10
            weights_with_inf[0] = np.inf
            o.sample(n=3, weights=weights_with_inf)

        with pytest.raises(ValueError):
            weights_with_ninf = [0.1] * 10
            weights_with_ninf[0] = -np.inf
            o.sample(n=3, weights=weights_with_ninf)

        # All zeros raises errors
        zero_weights = [0] * 10
        with pytest.raises(ValueError):
            o.sample(n=3, weights=zero_weights)

        # All missing weights
        nan_weights = [np.nan] * 10
        with pytest.raises(ValueError):
            o.sample(n=3, weights=nan_weights)

        # Check np.nan are replaced by zeros.
        weights_with_nan = [np.nan] * 10
        weights_with_nan[5] = 0.5
        self._compare(
            o.sample(n=1, axis=0, weights=weights_with_nan), o.iloc[5:6])

        # Check None are also replaced by zeros.
        weights_with_None = [None] * 10
        weights_with_None[5] = 0.5
        self._compare(
            o.sample(n=1, axis=0, weights=weights_with_None), o.iloc[5:6])

    def test_size_compat(self):
        # GH8846
        # size property should be defined

        o = self._construct(shape=10)
        assert o.size == np.prod(o.shape)
        assert o.size == 10 ** len(o.axes)

    def test_split_compat(self):
        # xref GH8846
        o = self._construct(shape=10)
        assert len(np.array_split(o, 5)) == 5
        assert len(np.array_split(o, 2)) == 2

    def test_unexpected_keyword(self):  # GH8597
        df = DataFrame(np.random.randn(5, 2), columns=['jim', 'joe'])
        ca = pd.Categorical([0, 0, 2, 2, 3, np.nan])
        ts = df['joe'].copy()
        ts[2] = np.nan

        with assert_raises_regex(TypeError, 'unexpected keyword'):
            df.drop('joe', axis=1, in_place=True)

        with assert_raises_regex(TypeError, 'unexpected keyword'):
            df.reindex([1, 0], inplace=True)

        with assert_raises_regex(TypeError, 'unexpected keyword'):
            ca.fillna(0, inplace=True)

        with assert_raises_regex(TypeError, 'unexpected keyword'):
            ts.fillna(0, in_place=True)

    # See gh-12301
    def test_stat_unexpected_keyword(self):
        obj = self._construct(5)
        starwars = 'Star Wars'
        errmsg = 'unexpected keyword'

        with assert_raises_regex(TypeError, errmsg):
            obj.max(epic=starwars)  # stat_function
        with assert_raises_regex(TypeError, errmsg):
            obj.var(epic=starwars)  # stat_function_ddof
        with assert_raises_regex(TypeError, errmsg):
            obj.sum(epic=starwars)  # cum_function
        with assert_raises_regex(TypeError, errmsg):
            obj.any(epic=starwars)  # logical_function

    def test_api_compat(self):

        # GH 12021
        # compat for __name__, __qualname__

        obj = self._construct(5)
        for func in ['sum', 'cumsum', 'any', 'var']:
            f = getattr(obj, func)
            assert f.__name__ == func
            if PY3:
                assert f.__qualname__.endswith(func)

    def test_stat_non_defaults_args(self):
        obj = self._construct(5)
        out = np.array([0])
        errmsg = "the 'out' parameter is not supported"

        with assert_raises_regex(ValueError, errmsg):
            obj.max(out=out)  # stat_function
        with assert_raises_regex(ValueError, errmsg):
            obj.var(out=out)  # stat_function_ddof
        with assert_raises_regex(ValueError, errmsg):
            obj.sum(out=out)  # cum_function
        with assert_raises_regex(ValueError, errmsg):
            obj.any(out=out)  # logical_function

    def test_truncate_out_of_bounds(self):
        # GH11382

        # small
        shape = [int(2e3)] + ([1] * (self._ndim - 1))
        small = self._construct(shape, dtype='int8')
        self._compare(small.truncate(), small)
        self._compare(small.truncate(before=0, after=3e3), small)
        self._compare(small.truncate(before=-1, after=2e3), small)

        # big
        shape = [int(2e6)] + ([1] * (self._ndim - 1))
        big = self._construct(shape, dtype='int8')
        self._compare(big.truncate(), big)
        self._compare(big.truncate(before=0, after=3e6), big)
        self._compare(big.truncate(before=-1, after=2e6), big)

    def test_validate_bool_args(self):
        df = DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]})
        invalid_values = [1, "True", [1, 2, 3], 5.0]

        for value in invalid_values:
            with pytest.raises(ValueError):
                super(DataFrame, df).rename_axis(mapper={'a': 'x', 'b': 'y'},
                                                 axis=1, inplace=value)

            with pytest.raises(ValueError):
                super(DataFrame, df).drop('a', axis=1, inplace=value)

            with pytest.raises(ValueError):
                super(DataFrame, df).sort_index(inplace=value)

            with pytest.raises(ValueError):
                super(DataFrame, df)._consolidate(inplace=value)

            with pytest.raises(ValueError):
                super(DataFrame, df).fillna(value=0, inplace=value)

            with pytest.raises(ValueError):
                super(DataFrame, df).replace(to_replace=1, value=7,
                                             inplace=value)

            with pytest.raises(ValueError):
                super(DataFrame, df).interpolate(inplace=value)

            with pytest.raises(ValueError):
                super(DataFrame, df)._where(cond=df.a > 2, inplace=value)

            with pytest.raises(ValueError):
                super(DataFrame, df).mask(cond=df.a > 2, inplace=value)

    def test_copy_and_deepcopy(self):
        # GH 15444
        for shape in [0, 1, 2]:
            obj = self._construct(shape)
            for func in [copy,
                         deepcopy,
                         lambda x: x.copy(deep=False),
                         lambda x: x.copy(deep=True)]:
                obj_copy = func(obj)
                assert obj_copy is not obj
                self._compare(obj_copy, obj)

    @pytest.mark.parametrize("periods,fill_method,limit,exp", [
        (1, "ffill", None, [np.nan, np.nan, np.nan, 1, 1, 1.5, 0, 0]),
        (1, "ffill", 1, [np.nan, np.nan, np.nan, 1, 1, 1.5, 0, np.nan]),
        (1, "bfill", None, [np.nan, 0, 0, 1, 1, 1.5, np.nan, np.nan]),
        (1, "bfill", 1, [np.nan, np.nan, 0, 1, 1, 1.5, np.nan, np.nan]),
        (-1, "ffill", None, [np.nan, np.nan, -.5, -.5, -.6, 0, 0, np.nan]),
        (-1, "ffill", 1, [np.nan, np.nan, -.5, -.5, -.6, 0, np.nan, np.nan]),
        (-1, "bfill", None, [0, 0, -.5, -.5, -.6, np.nan, np.nan, np.nan]),
        (-1, "bfill", 1, [np.nan, 0, -.5, -.5, -.6, np.nan, np.nan, np.nan])
    ])
    def test_pct_change(self, periods, fill_method, limit, exp):
        vals = [np.nan, np.nan, 1, 2, 4, 10, np.nan, np.nan]
        obj = self._typ(vals)
        func = getattr(obj, 'pct_change')
        res = func(periods=periods, fill_method=fill_method, limit=limit)
        if type(obj) is DataFrame:
            tm.assert_frame_equal(res, DataFrame(exp))
        else:
            tm.assert_series_equal(res, Series(exp))


class TestNDFrame(object):
    # tests that don't fit elsewhere

    def test_sample(sel):
        # Fixes issue: 2419
        # additional specific object based tests

        # A few dataframe test with degenerate weights.
        easy_weight_list = [0] * 10
        easy_weight_list[5] = 1

        df = pd.DataFrame({'col1': range(10, 20),
                           'col2': range(20, 30),
                           'colString': ['a'] * 10,
                           'easyweights': easy_weight_list})
        sample1 = df.sample(n=1, weights='easyweights')
        assert_frame_equal(sample1, df.iloc[5:6])

        # Ensure proper error if string given as weight for Series, panel, or
        # DataFrame with axis = 1.
        s = Series(range(10))
        with pytest.raises(ValueError):
            s.sample(n=3, weights='weight_column')

        with catch_warnings(record=True):
            panel = Panel(items=[0, 1, 2], major_axis=[2, 3, 4],
                          minor_axis=[3, 4, 5])
            with pytest.raises(ValueError):
                panel.sample(n=1, weights='weight_column')

        with pytest.raises(ValueError):
            df.sample(n=1, weights='weight_column', axis=1)

        # Check weighting key error
        with pytest.raises(KeyError):
            df.sample(n=3, weights='not_a_real_column_name')

        # Check that re-normalizes weights that don't sum to one.
        weights_less_than_1 = [0] * 10
        weights_less_than_1[0] = 0.5
        tm.assert_frame_equal(
            df.sample(n=1, weights=weights_less_than_1), df.iloc[:1])

        ###
        # Test axis argument
        ###

        # Test axis argument
        df = pd.DataFrame({'col1': range(10), 'col2': ['a'] * 10})
        second_column_weight = [0, 1]
        assert_frame_equal(
            df.sample(n=1, axis=1, weights=second_column_weight), df[['col2']])

        # Different axis arg types
        assert_frame_equal(df.sample(n=1, axis='columns',
                                     weights=second_column_weight),
                           df[['col2']])

        weight = [0] * 10
        weight[5] = 0.5
        assert_frame_equal(df.sample(n=1, axis='rows', weights=weight),
                           df.iloc[5:6])
        assert_frame_equal(df.sample(n=1, axis='index', weights=weight),
                           df.iloc[5:6])

        # Check out of range axis values
        with pytest.raises(ValueError):
            df.sample(n=1, axis=2)

        with pytest.raises(ValueError):
            df.sample(n=1, axis='not_a_name')

        with pytest.raises(ValueError):
            s = pd.Series(range(10))
            s.sample(n=1, axis=1)

        # Test weight length compared to correct axis
        with pytest.raises(ValueError):
            df.sample(n=1, axis=1, weights=[0.5] * 10)

        # Check weights with axis = 1
        easy_weight_list = [0] * 3
        easy_weight_list[2] = 1

        df = pd.DataFrame({'col1': range(10, 20),
                           'col2': range(20, 30),
                           'colString': ['a'] * 10})
        sample1 = df.sample(n=1, axis=1, weights=easy_weight_list)
        assert_frame_equal(sample1, df[['colString']])

        # Test default axes
        with catch_warnings(record=True):
            p = Panel(items=['a', 'b', 'c'], major_axis=[2, 4, 6],
                      minor_axis=[1, 3, 5])
            assert_panel_equal(
                p.sample(n=3, random_state=42), p.sample(n=3, axis=1,
                                                         random_state=42))
            assert_frame_equal(
                df.sample(n=3, random_state=42), df.sample(n=3, axis=0,
                                                           random_state=42))

        # Test that function aligns weights with frame
        df = DataFrame(
            {'col1': [5, 6, 7],
             'col2': ['a', 'b', 'c'], }, index=[9, 5, 3])
        s = Series([1, 0, 0], index=[3, 5, 9])
        assert_frame_equal(df.loc[[3]], df.sample(1, weights=s))

        # Weights have index values to be dropped because not in
        # sampled DataFrame
        s2 = Series([0.001, 0, 10000], index=[3, 5, 10])
        assert_frame_equal(df.loc[[3]], df.sample(1, weights=s2))

        # Weights have empty values to be filed with zeros
        s3 = Series([0.01, 0], index=[3, 5])
        assert_frame_equal(df.loc[[3]], df.sample(1, weights=s3))

        # No overlap in weight and sampled DataFrame indices
        s4 = Series([1, 0], index=[1, 2])
        with pytest.raises(ValueError):
            df.sample(1, weights=s4)

    def test_squeeze(self):
        # noop
        for s in [tm.makeFloatSeries(), tm.makeStringSeries(),
                  tm.makeObjectSeries()]:
            tm.assert_series_equal(s.squeeze(), s)
        for df in [tm.makeTimeDataFrame()]:
            tm.assert_frame_equal(df.squeeze(), df)
        with catch_warnings(record=True):
            for p in [tm.makePanel()]:
                tm.assert_panel_equal(p.squeeze(), p)

        # squeezing
        df = tm.makeTimeDataFrame().reindex(columns=['A'])
        tm.assert_series_equal(df.squeeze(), df['A'])

        with catch_warnings(record=True):
            p = tm.makePanel().reindex(items=['ItemA'])
            tm.assert_frame_equal(p.squeeze(), p['ItemA'])

            p = tm.makePanel().reindex(items=['ItemA'], minor_axis=['A'])
            tm.assert_series_equal(p.squeeze(), p.loc['ItemA', :, 'A'])

        # don't fail with 0 length dimensions GH11229 & GH8999
        empty_series = Series([], name='five')
        empty_frame = DataFrame([empty_series])
        with catch_warnings(record=True):
            empty_panel = Panel({'six': empty_frame})

        [tm.assert_series_equal(empty_series, higher_dim.squeeze())
         for higher_dim in [empty_series, empty_frame, empty_panel]]

        # axis argument
        df = tm.makeTimeDataFrame(nper=1).iloc[:, :1]
        assert df.shape == (1, 1)
        tm.assert_series_equal(df.squeeze(axis=0), df.iloc[0])
        tm.assert_series_equal(df.squeeze(axis='index'), df.iloc[0])
        tm.assert_series_equal(df.squeeze(axis=1), df.iloc[:, 0])
        tm.assert_series_equal(df.squeeze(axis='columns'), df.iloc[:, 0])
        assert df.squeeze() == df.iloc[0, 0]
        pytest.raises(ValueError, df.squeeze, axis=2)
        pytest.raises(ValueError, df.squeeze, axis='x')

        df = tm.makeTimeDataFrame(3)
        tm.assert_frame_equal(df.squeeze(axis=0), df)

    def test_numpy_squeeze(self):
        s = tm.makeFloatSeries()
        tm.assert_series_equal(np.squeeze(s), s)

        df = tm.makeTimeDataFrame().reindex(columns=['A'])
        tm.assert_series_equal(np.squeeze(df), df['A'])

    def test_transpose(self):
        msg = (r"transpose\(\) got multiple values for "
               r"keyword argument 'axes'")
        for s in [tm.makeFloatSeries(), tm.makeStringSeries(),
                  tm.makeObjectSeries()]:
            # calls implementation in pandas/core/base.py
            tm.assert_series_equal(s.transpose(), s)
        for df in [tm.makeTimeDataFrame()]:
            tm.assert_frame_equal(df.transpose().transpose(), df)

        with catch_warnings(record=True):
            for p in [tm.makePanel()]:
                tm.assert_panel_equal(p.transpose(2, 0, 1)
                                      .transpose(1, 2, 0), p)
                tm.assert_raises_regex(TypeError, msg, p.transpose,
                                       2, 0, 1, axes=(2, 0, 1))

    def test_numpy_transpose(self):
        msg = "the 'axes' parameter is not supported"

        s = tm.makeFloatSeries()
        tm.assert_series_equal(
            np.transpose(s), s)
        tm.assert_raises_regex(ValueError, msg,
                               np.transpose, s, axes=1)

        df = tm.makeTimeDataFrame()
        tm.assert_frame_equal(np.transpose(
            np.transpose(df)), df)
        tm.assert_raises_regex(ValueError, msg,
                               np.transpose, df, axes=1)

        with catch_warnings(record=True):
            p = tm.makePanel()
            tm.assert_panel_equal(np.transpose(
                np.transpose(p, axes=(2, 0, 1)),
                axes=(1, 2, 0)), p)

    def test_take(self):
        indices = [1, 5, -2, 6, 3, -1]
        for s in [tm.makeFloatSeries(), tm.makeStringSeries(),
                  tm.makeObjectSeries()]:
            out = s.take(indices)
            expected = Series(data=s.values.take(indices),
                              index=s.index.take(indices), dtype=s.dtype)
            tm.assert_series_equal(out, expected)
        for df in [tm.makeTimeDataFrame()]:
            out = df.take(indices)
            expected = DataFrame(data=df.values.take(indices, axis=0),
                                 index=df.index.take(indices),
                                 columns=df.columns)
            tm.assert_frame_equal(out, expected)

        indices = [-3, 2, 0, 1]
        with catch_warnings(record=True):
            for p in [tm.makePanel()]:
                out = p.take(indices)
                expected = Panel(data=p.values.take(indices, axis=0),
                                 items=p.items.take(indices),
                                 major_axis=p.major_axis,
                                 minor_axis=p.minor_axis)
                tm.assert_panel_equal(out, expected)

    def test_take_invalid_kwargs(self):
        indices = [-3, 2, 0, 1]
        s = tm.makeFloatSeries()
        df = tm.makeTimeDataFrame()

        with catch_warnings(record=True):
            p = tm.makePanel()

        for obj in (s, df, p):
            msg = r"take\(\) got an unexpected keyword argument 'foo'"
            tm.assert_raises_regex(TypeError, msg, obj.take,
                                   indices, foo=2)

            msg = "the 'out' parameter is not supported"
            tm.assert_raises_regex(ValueError, msg, obj.take,
                                   indices, out=indices)

            msg = "the 'mode' parameter is not supported"
            tm.assert_raises_regex(ValueError, msg, obj.take,
                                   indices, mode='clip')

    def test_equals(self):
        s1 = pd.Series([1, 2, 3], index=[0, 2, 1])
        s2 = s1.copy()
        assert s1.equals(s2)

        s1[1] = 99
        assert not s1.equals(s2)

        # NaNs compare as equal
        s1 = pd.Series([1, np.nan, 3, np.nan], index=[0, 2, 1, 3])
        s2 = s1.copy()
        assert s1.equals(s2)

        s2[0] = 9.9
        assert not s1.equals(s2)

        idx = MultiIndex.from_tuples([(0, 'a'), (1, 'b'), (2, 'c')])
        s1 = Series([1, 2, np.nan], index=idx)
        s2 = s1.copy()
        assert s1.equals(s2)

        # Add object dtype column with nans
        index = np.random.random(10)
        df1 = DataFrame(
            np.random.random(10, ), index=index, columns=['floats'])
        df1['text'] = 'the sky is so blue. we could use more chocolate.'.split(
        )
        df1['start'] = date_range('2000-1-1', periods=10, freq='T')
        df1['end'] = date_range('2000-1-1', periods=10, freq='D')
        df1['diff'] = df1['end'] - df1['start']
        df1['bool'] = (np.arange(10) % 3 == 0)
        df1.loc[::2] = np.nan
        df2 = df1.copy()
        assert df1['text'].equals(df2['text'])
        assert df1['start'].equals(df2['start'])
        assert df1['end'].equals(df2['end'])
        assert df1['diff'].equals(df2['diff'])
        assert df1['bool'].equals(df2['bool'])
        assert df1.equals(df2)
        assert not df1.equals(object)

        # different dtype
        different = df1.copy()
        different['floats'] = different['floats'].astype('float32')
        assert not df1.equals(different)

        # different index
        different_index = -index
        different = df2.set_index(different_index)
        assert not df1.equals(different)

        # different columns
        different = df2.copy()
        different.columns = df2.columns[::-1]
        assert not df1.equals(different)

        # DatetimeIndex
        index = pd.date_range('2000-1-1', periods=10, freq='T')
        df1 = df1.set_index(index)
        df2 = df1.copy()
        assert df1.equals(df2)

        # MultiIndex
        df3 = df1.set_index(['text'], append=True)
        df2 = df1.set_index(['text'], append=True)
        assert df3.equals(df2)

        df2 = df1.set_index(['floats'], append=True)
        assert not df3.equals(df2)

        # NaN in index
        df3 = df1.set_index(['floats'], append=True)
        df2 = df1.set_index(['floats'], append=True)
        assert df3.equals(df2)

        # GH 8437
        a = pd.Series([False, np.nan])
        b = pd.Series([False, np.nan])
        c = pd.Series(index=range(2))
        d = pd.Series(index=range(2))
        e = pd.Series(index=range(2))
        f = pd.Series(index=range(2))
        c[:-1] = d[:-1] = e[0] = f[0] = False
        assert a.equals(a)
        assert a.equals(b)
        assert a.equals(c)
        assert a.equals(d)
        assert a.equals(e)
        assert e.equals(f)

    def test_describe_raises(self):
        with catch_warnings(record=True):
            with pytest.raises(NotImplementedError):
                tm.makePanel().describe()

    def test_pipe(self):
        df = DataFrame({'A': [1, 2, 3]})
        f = lambda x, y: x ** y
        result = df.pipe(f, 2)
        expected = DataFrame({'A': [1, 4, 9]})
        assert_frame_equal(result, expected)

        result = df.A.pipe(f, 2)
        assert_series_equal(result, expected.A)

    def test_pipe_tuple(self):
        df = DataFrame({'A': [1, 2, 3]})
        f = lambda x, y: y
        result = df.pipe((f, 'y'), 0)
        assert_frame_equal(result, df)

        result = df.A.pipe((f, 'y'), 0)
        assert_series_equal(result, df.A)

    def test_pipe_tuple_error(self):
        df = DataFrame({"A": [1, 2, 3]})
        f = lambda x, y: y
        with pytest.raises(ValueError):
            df.pipe((f, 'y'), x=1, y=0)

        with pytest.raises(ValueError):
            df.A.pipe((f, 'y'), x=1, y=0)

    def test_pipe_panel(self):
        with catch_warnings(record=True):
            wp = Panel({'r1': DataFrame({"A": [1, 2, 3]})})
            f = lambda x, y: x + y
            result = wp.pipe(f, 2)
            expected = wp + 2
            assert_panel_equal(result, expected)

            result = wp.pipe((f, 'y'), x=1)
            expected = wp + 1
            assert_panel_equal(result, expected)

            with pytest.raises(ValueError):
                result = wp.pipe((f, 'y'), x=1, y=1)
