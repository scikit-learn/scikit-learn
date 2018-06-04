# coding=utf-8
# pylint: disable-msg=E1101,W0612

import pytest
import pytz

from collections import Iterable
from datetime import datetime, timedelta
import operator
from itertools import product, starmap

from numpy import nan
import numpy as np
import pandas as pd

from pandas import (Index, Series, DataFrame, isna, bdate_range,
                    NaT, date_range, timedelta_range, Categorical)
from pandas.core.indexes.datetimes import Timestamp
from pandas.core.indexes.timedeltas import Timedelta
import pandas.core.nanops as nanops

from pandas.errors import PerformanceWarning
from pandas.compat import range, zip
from pandas import compat
from pandas.util.testing import (assert_series_equal, assert_almost_equal,
                                 assert_frame_equal, assert_index_equal)
import pandas.util.testing as tm

from .common import TestData


class TestSeriesComparisons(object):
    def test_series_comparison_scalars(self):
        series = Series(date_range('1/1/2000', periods=10))

        val = datetime(2000, 1, 4)
        result = series > val
        expected = Series([x > val for x in series])
        tm.assert_series_equal(result, expected)

        val = series[5]
        result = series > val
        expected = Series([x > val for x in series])
        tm.assert_series_equal(result, expected)

    def test_comparisons(self):
        left = np.random.randn(10)
        right = np.random.randn(10)
        left[:3] = np.nan

        result = nanops.nangt(left, right)
        with np.errstate(invalid='ignore'):
            expected = (left > right).astype('O')
        expected[:3] = np.nan

        assert_almost_equal(result, expected)

        s = Series(['a', 'b', 'c'])
        s2 = Series([False, True, False])

        # it works!
        exp = Series([False, False, False])
        assert_series_equal(s == s2, exp)
        assert_series_equal(s2 == s, exp)

    def test_operator_series_comparison_zerorank(self):
        # GH 13006
        result = np.float64(0) > pd.Series([1, 2, 3])
        expected = 0.0 > pd.Series([1, 2, 3])
        tm.assert_series_equal(result, expected)
        result = pd.Series([1, 2, 3]) < np.float64(0)
        expected = pd.Series([1, 2, 3]) < 0.0
        tm.assert_series_equal(result, expected)
        result = np.array([0, 1, 2])[0] > pd.Series([0, 1, 2])
        expected = 0.0 > pd.Series([1, 2, 3])
        tm.assert_series_equal(result, expected)

    def test_object_comparisons(self):
        s = Series(['a', 'b', np.nan, 'c', 'a'])

        result = s == 'a'
        expected = Series([True, False, False, False, True])
        assert_series_equal(result, expected)

        result = s < 'a'
        expected = Series([False, False, False, False, False])
        assert_series_equal(result, expected)

        result = s != 'a'
        expected = -(s == 'a')
        assert_series_equal(result, expected)

    def test_categorical_comparisons(self):
        # GH 8938
        # allow equality comparisons
        a = Series(list('abc'), dtype="category")
        b = Series(list('abc'), dtype="object")
        c = Series(['a', 'b', 'cc'], dtype="object")
        d = Series(list('acb'), dtype="object")
        e = Categorical(list('abc'))
        f = Categorical(list('acb'))

        # vs scalar
        assert not (a == 'a').all()
        assert ((a != 'a') == ~(a == 'a')).all()

        assert not ('a' == a).all()
        assert (a == 'a')[0]
        assert ('a' == a)[0]
        assert not ('a' != a)[0]

        # vs list-like
        assert (a == a).all()
        assert not (a != a).all()

        assert (a == list(a)).all()
        assert (a == b).all()
        assert (b == a).all()
        assert ((~(a == b)) == (a != b)).all()
        assert ((~(b == a)) == (b != a)).all()

        assert not (a == c).all()
        assert not (c == a).all()
        assert not (a == d).all()
        assert not (d == a).all()

        # vs a cat-like
        assert (a == e).all()
        assert (e == a).all()
        assert not (a == f).all()
        assert not (f == a).all()

        assert ((~(a == e) == (a != e)).all())
        assert ((~(e == a) == (e != a)).all())
        assert ((~(a == f) == (a != f)).all())
        assert ((~(f == a) == (f != a)).all())

        # non-equality is not comparable
        pytest.raises(TypeError, lambda: a < b)
        pytest.raises(TypeError, lambda: b < a)
        pytest.raises(TypeError, lambda: a > b)
        pytest.raises(TypeError, lambda: b > a)

    def test_comparison_tuples(self):
        # GH11339
        # comparisons vs tuple
        s = Series([(1, 1), (1, 2)])

        result = s == (1, 2)
        expected = Series([False, True])
        assert_series_equal(result, expected)

        result = s != (1, 2)
        expected = Series([True, False])
        assert_series_equal(result, expected)

        result = s == (0, 0)
        expected = Series([False, False])
        assert_series_equal(result, expected)

        result = s != (0, 0)
        expected = Series([True, True])
        assert_series_equal(result, expected)

        s = Series([(1, 1), (1, 1)])

        result = s == (1, 1)
        expected = Series([True, True])
        assert_series_equal(result, expected)

        result = s != (1, 1)
        expected = Series([False, False])
        assert_series_equal(result, expected)

        s = Series([frozenset([1]), frozenset([1, 2])])

        result = s == frozenset([1])
        expected = Series([True, False])
        assert_series_equal(result, expected)

    def test_comparison_operators_with_nas(self):
        ser = Series(bdate_range('1/1/2000', periods=10), dtype=object)
        ser[::2] = np.nan

        # test that comparisons work
        ops = ['lt', 'le', 'gt', 'ge', 'eq', 'ne']
        for op in ops:
            val = ser[5]

            f = getattr(operator, op)
            result = f(ser, val)

            expected = f(ser.dropna(), val).reindex(ser.index)

            if op == 'ne':
                expected = expected.fillna(True).astype(bool)
            else:
                expected = expected.fillna(False).astype(bool)

            assert_series_equal(result, expected)

            # fffffffuuuuuuuuuuuu
            # result = f(val, s)
            # expected = f(val, s.dropna()).reindex(s.index)
            # assert_series_equal(result, expected)

            # boolean &, |, ^ should work with object arrays and propagate NAs

        ops = ['and_', 'or_', 'xor']
        mask = ser.isna()
        for bool_op in ops:
            func = getattr(operator, bool_op)

            filled = ser.fillna(ser[0])

            result = func(ser < ser[9], ser > ser[3])

            expected = func(filled < filled[9], filled > filled[3])
            expected[mask] = False
            assert_series_equal(result, expected)

    def test_comparison_object_numeric_nas(self):
        ser = Series(np.random.randn(10), dtype=object)
        shifted = ser.shift(2)

        ops = ['lt', 'le', 'gt', 'ge', 'eq', 'ne']
        for op in ops:
            func = getattr(operator, op)

            result = func(ser, shifted)
            expected = func(ser.astype(float), shifted.astype(float))
            assert_series_equal(result, expected)

    def test_comparison_invalid(self):
        # GH4968
        # invalid date/int comparisons
        s = Series(range(5))
        s2 = Series(date_range('20010101', periods=5))

        for (x, y) in [(s, s2), (s2, s)]:
            pytest.raises(TypeError, lambda: x == y)
            pytest.raises(TypeError, lambda: x != y)
            pytest.raises(TypeError, lambda: x >= y)
            pytest.raises(TypeError, lambda: x > y)
            pytest.raises(TypeError, lambda: x < y)
            pytest.raises(TypeError, lambda: x <= y)

    def test_unequal_categorical_comparison_raises_type_error(self):
        # unequal comparison should raise for unordered cats
        cat = Series(Categorical(list("abc")))

        def f():
            cat > "b"

        pytest.raises(TypeError, f)
        cat = Series(Categorical(list("abc"), ordered=False))

        def f():
            cat > "b"

        pytest.raises(TypeError, f)

        # https://github.com/pandas-dev/pandas/issues/9836#issuecomment-92123057
        # and following comparisons with scalars not in categories should raise
        # for unequal comps, but not for equal/not equal
        cat = Series(Categorical(list("abc"), ordered=True))

        pytest.raises(TypeError, lambda: cat < "d")
        pytest.raises(TypeError, lambda: cat > "d")
        pytest.raises(TypeError, lambda: "d" < cat)
        pytest.raises(TypeError, lambda: "d" > cat)

        tm.assert_series_equal(cat == "d", Series([False, False, False]))
        tm.assert_series_equal(cat != "d", Series([True, True, True]))

    @pytest.mark.parametrize('dtype', [None, object])
    def test_more_na_comparisons(self, dtype):
        left = Series(['a', np.nan, 'c'], dtype=dtype)
        right = Series(['a', np.nan, 'd'], dtype=dtype)

        result = left == right
        expected = Series([True, False, False])
        assert_series_equal(result, expected)

        result = left != right
        expected = Series([False, True, True])
        assert_series_equal(result, expected)

        result = left == np.nan
        expected = Series([False, False, False])
        assert_series_equal(result, expected)

        result = left != np.nan
        expected = Series([True, True, True])
        assert_series_equal(result, expected)

    @pytest.mark.parametrize('pair', [
        ([pd.Timestamp('2011-01-01'), NaT, pd.Timestamp('2011-01-03')],
         [NaT, NaT, pd.Timestamp('2011-01-03')]),

        ([pd.Timedelta('1 days'), NaT, pd.Timedelta('3 days')],
         [NaT, NaT, pd.Timedelta('3 days')]),

        ([pd.Period('2011-01', freq='M'), NaT, pd.Period('2011-03', freq='M')],
         [NaT, NaT, pd.Period('2011-03', freq='M')])])
    @pytest.mark.parametrize('reverse', [True, False])
    @pytest.mark.parametrize('box', [Series, Index])
    @pytest.mark.parametrize('dtype', [None, object])
    def test_nat_comparisons(self, dtype, box, reverse, pair):
        l, r = pair
        if reverse:
            # add lhs / rhs switched data
            l, r = r, l

        left = Series(l, dtype=dtype)
        right = box(r, dtype=dtype)
        # Series, Index

        expected = Series([False, False, True])
        assert_series_equal(left == right, expected)

        expected = Series([True, True, False])
        assert_series_equal(left != right, expected)

        expected = Series([False, False, False])
        assert_series_equal(left < right, expected)

        expected = Series([False, False, False])
        assert_series_equal(left > right, expected)

        expected = Series([False, False, True])
        assert_series_equal(left >= right, expected)

        expected = Series([False, False, True])
        assert_series_equal(left <= right, expected)

    @pytest.mark.parametrize('data', [
        [pd.Timestamp('2011-01-01'), NaT, pd.Timestamp('2011-01-03')],
        [pd.Timedelta('1 days'), NaT, pd.Timedelta('3 days')],
        [pd.Period('2011-01', freq='M'), NaT, pd.Period('2011-03', freq='M')]
    ])
    @pytest.mark.parametrize('dtype', [None, object])
    def test_nat_comparisons_scalar(self, dtype, data):
        left = Series(data, dtype=dtype)

        expected = Series([False, False, False])
        assert_series_equal(left == pd.NaT, expected)
        assert_series_equal(pd.NaT == left, expected)

        expected = Series([True, True, True])
        assert_series_equal(left != pd.NaT, expected)
        assert_series_equal(pd.NaT != left, expected)

        expected = Series([False, False, False])
        assert_series_equal(left < pd.NaT, expected)
        assert_series_equal(pd.NaT > left, expected)
        assert_series_equal(left <= pd.NaT, expected)
        assert_series_equal(pd.NaT >= left, expected)

        assert_series_equal(left > pd.NaT, expected)
        assert_series_equal(pd.NaT < left, expected)
        assert_series_equal(left >= pd.NaT, expected)
        assert_series_equal(pd.NaT <= left, expected)

    def test_comparison_different_length(self):
        a = Series(['a', 'b', 'c'])
        b = Series(['b', 'a'])
        pytest.raises(ValueError, a.__lt__, b)

        a = Series([1, 2])
        b = Series([2, 3, 4])
        pytest.raises(ValueError, a.__eq__, b)

    def test_comparison_label_based(self):

        # GH 4947
        # comparisons should be label based

        a = Series([True, False, True], list('bca'))
        b = Series([False, True, False], list('abc'))

        expected = Series([False, True, False], list('abc'))
        result = a & b
        assert_series_equal(result, expected)

        expected = Series([True, True, False], list('abc'))
        result = a | b
        assert_series_equal(result, expected)

        expected = Series([True, False, False], list('abc'))
        result = a ^ b
        assert_series_equal(result, expected)

        # rhs is bigger
        a = Series([True, False, True], list('bca'))
        b = Series([False, True, False, True], list('abcd'))

        expected = Series([False, True, False, False], list('abcd'))
        result = a & b
        assert_series_equal(result, expected)

        expected = Series([True, True, False, False], list('abcd'))
        result = a | b
        assert_series_equal(result, expected)

        # filling

        # vs empty
        result = a & Series([])
        expected = Series([False, False, False], list('bca'))
        assert_series_equal(result, expected)

        result = a | Series([])
        expected = Series([True, False, True], list('bca'))
        assert_series_equal(result, expected)

        # vs non-matching
        result = a & Series([1], ['z'])
        expected = Series([False, False, False, False], list('abcz'))
        assert_series_equal(result, expected)

        result = a | Series([1], ['z'])
        expected = Series([True, True, False, False], list('abcz'))
        assert_series_equal(result, expected)

        # identity
        # we would like s[s|e] == s to hold for any e, whether empty or not
        for e in [Series([]), Series([1], ['z']),
                  Series(np.nan, b.index), Series(np.nan, a.index)]:
            result = a[a | e]
            assert_series_equal(result, a[a])

        for e in [Series(['z'])]:
            if compat.PY3:
                with tm.assert_produces_warning(RuntimeWarning):
                    result = a[a | e]
            else:
                result = a[a | e]
            assert_series_equal(result, a[a])

        # vs scalars
        index = list('bca')
        t = Series([True, False, True])

        for v in [True, 1, 2]:
            result = Series([True, False, True], index=index) | v
            expected = Series([True, True, True], index=index)
            assert_series_equal(result, expected)

        for v in [np.nan, 'foo']:
            pytest.raises(TypeError, lambda: t | v)

        for v in [False, 0]:
            result = Series([True, False, True], index=index) | v
            expected = Series([True, False, True], index=index)
            assert_series_equal(result, expected)

        for v in [True, 1]:
            result = Series([True, False, True], index=index) & v
            expected = Series([True, False, True], index=index)
            assert_series_equal(result, expected)

        for v in [False, 0]:
            result = Series([True, False, True], index=index) & v
            expected = Series([False, False, False], index=index)
            assert_series_equal(result, expected)
        for v in [np.nan]:
            pytest.raises(TypeError, lambda: t & v)

    def test_comparison_flex_basic(self):
        left = pd.Series(np.random.randn(10))
        right = pd.Series(np.random.randn(10))

        assert_series_equal(left.eq(right), left == right)
        assert_series_equal(left.ne(right), left != right)
        assert_series_equal(left.le(right), left < right)
        assert_series_equal(left.lt(right), left <= right)
        assert_series_equal(left.gt(right), left > right)
        assert_series_equal(left.ge(right), left >= right)

        # axis
        for axis in [0, None, 'index']:
            assert_series_equal(left.eq(right, axis=axis), left == right)
            assert_series_equal(left.ne(right, axis=axis), left != right)
            assert_series_equal(left.le(right, axis=axis), left < right)
            assert_series_equal(left.lt(right, axis=axis), left <= right)
            assert_series_equal(left.gt(right, axis=axis), left > right)
            assert_series_equal(left.ge(right, axis=axis), left >= right)

        #
        msg = 'No axis named 1 for object type'
        for op in ['eq', 'ne', 'le', 'le', 'gt', 'ge']:
            with tm.assert_raises_regex(ValueError, msg):
                getattr(left, op)(right, axis=1)

    def test_comparison_flex_alignment(self):
        left = Series([1, 3, 2], index=list('abc'))
        right = Series([2, 2, 2], index=list('bcd'))

        exp = pd.Series([False, False, True, False], index=list('abcd'))
        assert_series_equal(left.eq(right), exp)

        exp = pd.Series([True, True, False, True], index=list('abcd'))
        assert_series_equal(left.ne(right), exp)

        exp = pd.Series([False, False, True, False], index=list('abcd'))
        assert_series_equal(left.le(right), exp)

        exp = pd.Series([False, False, False, False], index=list('abcd'))
        assert_series_equal(left.lt(right), exp)

        exp = pd.Series([False, True, True, False], index=list('abcd'))
        assert_series_equal(left.ge(right), exp)

        exp = pd.Series([False, True, False, False], index=list('abcd'))
        assert_series_equal(left.gt(right), exp)

    def test_comparison_flex_alignment_fill(self):
        left = Series([1, 3, 2], index=list('abc'))
        right = Series([2, 2, 2], index=list('bcd'))

        exp = pd.Series([False, False, True, True], index=list('abcd'))
        assert_series_equal(left.eq(right, fill_value=2), exp)

        exp = pd.Series([True, True, False, False], index=list('abcd'))
        assert_series_equal(left.ne(right, fill_value=2), exp)

        exp = pd.Series([False, False, True, True], index=list('abcd'))
        assert_series_equal(left.le(right, fill_value=0), exp)

        exp = pd.Series([False, False, False, True], index=list('abcd'))
        assert_series_equal(left.lt(right, fill_value=0), exp)

        exp = pd.Series([True, True, True, False], index=list('abcd'))
        assert_series_equal(left.ge(right, fill_value=0), exp)

        exp = pd.Series([True, True, False, False], index=list('abcd'))
        assert_series_equal(left.gt(right, fill_value=0), exp)

    def test_ne(self):
        ts = Series([3, 4, 5, 6, 7], [3, 4, 5, 6, 7], dtype=float)
        expected = [True, True, False, True, True]
        assert tm.equalContents(ts.index != 5, expected)
        assert tm.equalContents(~(ts.index == 5), expected)

    def test_comp_ops_df_compat(self):
        # GH 1134
        s1 = pd.Series([1, 2, 3], index=list('ABC'), name='x')
        s2 = pd.Series([2, 2, 2], index=list('ABD'), name='x')

        s3 = pd.Series([1, 2, 3], index=list('ABC'), name='x')
        s4 = pd.Series([2, 2, 2, 2], index=list('ABCD'), name='x')

        for left, right in [(s1, s2), (s2, s1), (s3, s4), (s4, s3)]:

            msg = "Can only compare identically-labeled Series objects"
            with tm.assert_raises_regex(ValueError, msg):
                left == right

            with tm.assert_raises_regex(ValueError, msg):
                left != right

            with tm.assert_raises_regex(ValueError, msg):
                left < right

            msg = "Can only compare identically-labeled DataFrame objects"
            with tm.assert_raises_regex(ValueError, msg):
                left.to_frame() == right.to_frame()

            with tm.assert_raises_regex(ValueError, msg):
                left.to_frame() != right.to_frame()

            with tm.assert_raises_regex(ValueError, msg):
                left.to_frame() < right.to_frame()


class TestTimedeltaSeriesArithmetic(object):

    def test_timedelta64_operations_with_DateOffset(self):
        # GH 10699
        td = Series([timedelta(minutes=5, seconds=3)] * 3)
        result = td + pd.offsets.Minute(1)
        expected = Series([timedelta(minutes=6, seconds=3)] * 3)
        assert_series_equal(result, expected)

        result = td - pd.offsets.Minute(1)
        expected = Series([timedelta(minutes=4, seconds=3)] * 3)
        assert_series_equal(result, expected)

        with tm.assert_produces_warning(PerformanceWarning):
            result = td + Series([pd.offsets.Minute(1), pd.offsets.Second(3),
                                  pd.offsets.Hour(2)])
        expected = Series([timedelta(minutes=6, seconds=3), timedelta(
            minutes=5, seconds=6), timedelta(hours=2, minutes=5, seconds=3)])
        assert_series_equal(result, expected)

        result = td + pd.offsets.Minute(1) + pd.offsets.Second(12)
        expected = Series([timedelta(minutes=6, seconds=15)] * 3)
        assert_series_equal(result, expected)

        # valid DateOffsets
        for do in ['Hour', 'Minute', 'Second', 'Day', 'Micro', 'Milli',
                   'Nano']:
            op = getattr(pd.offsets, do)
            td + op(5)
            op(5) + td
            td - op(5)
            op(5) - td

    def test_timedelta64_operations_with_timedeltas(self):
        # td operate with td
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td2 = timedelta(minutes=5, seconds=4)
        result = td1 - td2
        expected = (Series([timedelta(seconds=0)] * 3) -
                    Series([timedelta(seconds=1)] * 3))
        assert result.dtype == 'm8[ns]'
        assert_series_equal(result, expected)

        result2 = td2 - td1
        expected = (Series([timedelta(seconds=1)] * 3) -
                    Series([timedelta(seconds=0)] * 3))
        assert_series_equal(result2, expected)

        # roundtrip
        assert_series_equal(result + td2, td1)

        # Now again, using pd.to_timedelta, which should build
        # a Series or a scalar, depending on input.
        td1 = Series(pd.to_timedelta(['00:05:03'] * 3))
        td2 = pd.to_timedelta('00:05:04')
        result = td1 - td2
        expected = (Series([timedelta(seconds=0)] * 3) -
                    Series([timedelta(seconds=1)] * 3))
        assert result.dtype == 'm8[ns]'
        assert_series_equal(result, expected)

        result2 = td2 - td1
        expected = (Series([timedelta(seconds=1)] * 3) -
                    Series([timedelta(seconds=0)] * 3))
        assert_series_equal(result2, expected)

        # roundtrip
        assert_series_equal(result + td2, td1)

    def test_operators_timedelta64(self):
        # series ops
        v1 = date_range('2012-1-1', periods=3, freq='D')
        v2 = date_range('2012-1-2', periods=3, freq='D')
        rs = Series(v2) - Series(v1)
        xp = Series(1e9 * 3600 * 24,
                    rs.index).astype('int64').astype('timedelta64[ns]')
        assert_series_equal(rs, xp)
        assert rs.dtype == 'timedelta64[ns]'

        df = DataFrame(dict(A=v1))
        td = Series([timedelta(days=i) for i in range(3)])
        assert td.dtype == 'timedelta64[ns]'

        # series on the rhs
        result = df['A'] - df['A'].shift()
        assert result.dtype == 'timedelta64[ns]'

        result = df['A'] + td
        assert result.dtype == 'M8[ns]'

        # scalar Timestamp on rhs
        maxa = df['A'].max()
        assert isinstance(maxa, Timestamp)

        resultb = df['A'] - df['A'].max()
        assert resultb.dtype == 'timedelta64[ns]'

        # timestamp on lhs
        result = resultb + df['A']
        values = [Timestamp('20111230'), Timestamp('20120101'),
                  Timestamp('20120103')]
        expected = Series(values, name='A')
        assert_series_equal(result, expected)

        # datetimes on rhs
        result = df['A'] - datetime(2001, 1, 1)
        expected = Series(
            [timedelta(days=4017 + i) for i in range(3)], name='A')
        assert_series_equal(result, expected)
        assert result.dtype == 'm8[ns]'

        d = datetime(2001, 1, 1, 3, 4)
        resulta = df['A'] - d
        assert resulta.dtype == 'm8[ns]'

        # roundtrip
        resultb = resulta + d
        assert_series_equal(df['A'], resultb)

        # timedeltas on rhs
        td = timedelta(days=1)
        resulta = df['A'] + td
        resultb = resulta - td
        assert_series_equal(resultb, df['A'])
        assert resultb.dtype == 'M8[ns]'

        # roundtrip
        td = timedelta(minutes=5, seconds=3)
        resulta = df['A'] + td
        resultb = resulta - td
        assert_series_equal(df['A'], resultb)
        assert resultb.dtype == 'M8[ns]'

        # inplace
        value = rs[2] + np.timedelta64(timedelta(minutes=5, seconds=1))
        rs[2] += np.timedelta64(timedelta(minutes=5, seconds=1))
        assert rs[2] == value

    def test_timedelta64_ops_nat(self):
        # GH 11349
        timedelta_series = Series([NaT, Timedelta('1s')])
        nat_series_dtype_timedelta = Series([NaT, NaT],
                                            dtype='timedelta64[ns]')
        single_nat_dtype_timedelta = Series([NaT], dtype='timedelta64[ns]')

        # subtraction
        assert_series_equal(timedelta_series - NaT,
                            nat_series_dtype_timedelta)
        assert_series_equal(-NaT + timedelta_series,
                            nat_series_dtype_timedelta)

        assert_series_equal(timedelta_series - single_nat_dtype_timedelta,
                            nat_series_dtype_timedelta)
        assert_series_equal(-single_nat_dtype_timedelta + timedelta_series,
                            nat_series_dtype_timedelta)

        # addition
        assert_series_equal(nat_series_dtype_timedelta + NaT,
                            nat_series_dtype_timedelta)
        assert_series_equal(NaT + nat_series_dtype_timedelta,
                            nat_series_dtype_timedelta)

        assert_series_equal(nat_series_dtype_timedelta +
                            single_nat_dtype_timedelta,
                            nat_series_dtype_timedelta)
        assert_series_equal(single_nat_dtype_timedelta +
                            nat_series_dtype_timedelta,
                            nat_series_dtype_timedelta)

        assert_series_equal(timedelta_series + NaT,
                            nat_series_dtype_timedelta)
        assert_series_equal(NaT + timedelta_series,
                            nat_series_dtype_timedelta)

        assert_series_equal(timedelta_series + single_nat_dtype_timedelta,
                            nat_series_dtype_timedelta)
        assert_series_equal(single_nat_dtype_timedelta + timedelta_series,
                            nat_series_dtype_timedelta)

        assert_series_equal(nat_series_dtype_timedelta + NaT,
                            nat_series_dtype_timedelta)
        assert_series_equal(NaT + nat_series_dtype_timedelta,
                            nat_series_dtype_timedelta)

        assert_series_equal(nat_series_dtype_timedelta +
                            single_nat_dtype_timedelta,
                            nat_series_dtype_timedelta)
        assert_series_equal(single_nat_dtype_timedelta +
                            nat_series_dtype_timedelta,
                            nat_series_dtype_timedelta)

        # multiplication
        assert_series_equal(nat_series_dtype_timedelta * 1.0,
                            nat_series_dtype_timedelta)
        assert_series_equal(1.0 * nat_series_dtype_timedelta,
                            nat_series_dtype_timedelta)

        assert_series_equal(timedelta_series * 1, timedelta_series)
        assert_series_equal(1 * timedelta_series, timedelta_series)

        assert_series_equal(timedelta_series * 1.5,
                            Series([NaT, Timedelta('1.5s')]))
        assert_series_equal(1.5 * timedelta_series,
                            Series([NaT, Timedelta('1.5s')]))

        assert_series_equal(timedelta_series * nan,
                            nat_series_dtype_timedelta)
        assert_series_equal(nan * timedelta_series,
                            nat_series_dtype_timedelta)

        # division
        assert_series_equal(timedelta_series / 2,
                            Series([NaT, Timedelta('0.5s')]))
        assert_series_equal(timedelta_series / 2.0,
                            Series([NaT, Timedelta('0.5s')]))
        assert_series_equal(timedelta_series / nan,
                            nat_series_dtype_timedelta)

    @pytest.mark.parametrize('scalar_td', [timedelta(minutes=5, seconds=4),
                                           Timedelta(minutes=5, seconds=4),
                                           Timedelta('5m4s').to_timedelta64()])
    def test_operators_timedelta64_with_timedelta(self, scalar_td):
        # smoke tests
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        td1 + scalar_td
        scalar_td + td1
        td1 - scalar_td
        scalar_td - td1
        td1 / scalar_td
        scalar_td / td1


class TestDatetimeSeriesArithmetic(object):
    @pytest.mark.parametrize(
        'box, assert_func',
        [(Series, tm.assert_series_equal),
         (pd.Index, tm.assert_index_equal)])
    def test_sub_datetime64_not_ns(self, box, assert_func):
        # GH#7996
        dt64 = np.datetime64('2013-01-01')
        assert dt64.dtype == 'datetime64[D]'

        obj = box(date_range('20130101', periods=3))
        res = obj - dt64
        expected = box([Timedelta(days=0), Timedelta(days=1),
                        Timedelta(days=2)])
        assert_func(res, expected)

        res = dt64 - obj
        assert_func(res, -expected)

    def test_operators_datetimelike_invalid(self, all_arithmetic_operators):
        # these are all TypeEror ops
        op_str = all_arithmetic_operators

        def check(get_ser, test_ser):

            # check that we are getting a TypeError
            # with 'operate' (from core/ops.py) for the ops that are not
            # defined
            op = getattr(get_ser, op_str, None)
            with tm.assert_raises_regex(TypeError, 'operate|cannot'):
                op(test_ser)

        # ## timedelta64 ###
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        # ## datetime64 ###
        dt1 = Series([Timestamp('20111230'), Timestamp('20120101'),
                      Timestamp('20120103')])
        dt1.iloc[2] = np.nan
        dt2 = Series([Timestamp('20111231'), Timestamp('20120102'),
                      Timestamp('20120104')])
        if op_str not in ['__sub__', '__rsub__']:
            check(dt1, dt2)

        # ## datetime64 with timetimedelta ###
        # TODO(jreback) __rsub__ should raise?
        if op_str not in ['__add__', '__radd__', '__sub__']:
            check(dt1, td1)

        # 8260, 10763
        # datetime64 with tz
        tz = 'US/Eastern'
        dt1 = Series(date_range('2000-01-01 09:00:00', periods=5,
                                tz=tz), name='foo')
        dt2 = dt1.copy()
        dt2.iloc[2] = np.nan
        td1 = Series(timedelta_range('1 days 1 min', periods=5, freq='H'))
        td2 = td1.copy()
        td2.iloc[1] = np.nan

        if op_str not in ['__add__', '__radd__', '__sub__', '__rsub__']:
            check(dt2, td2)

    def test_operators_datetimelike(self):

        # ## timedelta64 ###
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        # ## datetime64 ###
        dt1 = Series([Timestamp('20111230'), Timestamp('20120101'),
                      Timestamp('20120103')])
        dt1.iloc[2] = np.nan
        dt2 = Series([Timestamp('20111231'), Timestamp('20120102'),
                      Timestamp('20120104')])
        dt1 - dt2
        dt2 - dt1

        # ## datetime64 with timetimedelta ###
        dt1 + td1
        td1 + dt1
        dt1 - td1
        # TODO: Decide if this ought to work.
        # td1 - dt1

        # ## timetimedelta with datetime64 ###
        td1 + dt1
        dt1 + td1

    def test_operators_datetimelike_with_timezones(self):

        tz = 'US/Eastern'
        dt1 = Series(date_range('2000-01-01 09:00:00', periods=5,
                                tz=tz), name='foo')
        dt2 = dt1.copy()
        dt2.iloc[2] = np.nan

        td1 = Series(timedelta_range('1 days 1 min', periods=5, freq='H'))
        td2 = td1.copy()
        td2.iloc[1] = np.nan

        result = dt1 + td1[0]
        exp = (dt1.dt.tz_localize(None) + td1[0]).dt.tz_localize(tz)
        assert_series_equal(result, exp)

        result = dt2 + td2[0]
        exp = (dt2.dt.tz_localize(None) + td2[0]).dt.tz_localize(tz)
        assert_series_equal(result, exp)

        # odd numpy behavior with scalar timedeltas
        result = td1[0] + dt1
        exp = (dt1.dt.tz_localize(None) + td1[0]).dt.tz_localize(tz)
        assert_series_equal(result, exp)

        result = td2[0] + dt2
        exp = (dt2.dt.tz_localize(None) + td2[0]).dt.tz_localize(tz)
        assert_series_equal(result, exp)

        result = dt1 - td1[0]
        exp = (dt1.dt.tz_localize(None) - td1[0]).dt.tz_localize(tz)
        assert_series_equal(result, exp)
        pytest.raises(TypeError, lambda: td1[0] - dt1)

        result = dt2 - td2[0]
        exp = (dt2.dt.tz_localize(None) - td2[0]).dt.tz_localize(tz)
        assert_series_equal(result, exp)
        pytest.raises(TypeError, lambda: td2[0] - dt2)

        result = dt1 + td1
        exp = (dt1.dt.tz_localize(None) + td1).dt.tz_localize(tz)
        assert_series_equal(result, exp)

        result = dt2 + td2
        exp = (dt2.dt.tz_localize(None) + td2).dt.tz_localize(tz)
        assert_series_equal(result, exp)

        result = dt1 - td1
        exp = (dt1.dt.tz_localize(None) - td1).dt.tz_localize(tz)
        assert_series_equal(result, exp)

        result = dt2 - td2
        exp = (dt2.dt.tz_localize(None) - td2).dt.tz_localize(tz)
        assert_series_equal(result, exp)

        pytest.raises(TypeError, lambda: td1 - dt1)
        pytest.raises(TypeError, lambda: td2 - dt2)

    def test_sub_single_tz(self):
        # GH12290
        s1 = Series([pd.Timestamp('2016-02-10', tz='America/Sao_Paulo')])
        s2 = Series([pd.Timestamp('2016-02-08', tz='America/Sao_Paulo')])
        result = s1 - s2
        expected = Series([Timedelta('2days')])
        assert_series_equal(result, expected)
        result = s2 - s1
        expected = Series([Timedelta('-2days')])
        assert_series_equal(result, expected)

    def test_dt64tz_series_sub_dtitz(self):
        # GH#19071 subtracting tzaware DatetimeIndex from tzaware Series
        # (with same tz) raises, fixed by #19024
        dti = pd.date_range('1999-09-30', periods=10, tz='US/Pacific')
        ser = pd.Series(dti)
        expected = pd.Series(pd.TimedeltaIndex(['0days'] * 10))

        res = dti - ser
        tm.assert_series_equal(res, expected)
        res = ser - dti
        tm.assert_series_equal(res, expected)

    def test_sub_datetime_compat(self):
        # see gh-14088
        s = Series([datetime(2016, 8, 23, 12, tzinfo=pytz.utc), pd.NaT])
        dt = datetime(2016, 8, 22, 12, tzinfo=pytz.utc)
        exp = Series([Timedelta('1 days'), pd.NaT])
        assert_series_equal(s - dt, exp)
        assert_series_equal(s - Timestamp(dt), exp)

    def test_dt64_series_with_timedelta(self):
        # scalar timedeltas/np.timedelta64 objects
        # operate with np.timedelta64 correctly
        s = Series([Timestamp('20130101 9:01'), Timestamp('20130101 9:02')])

        result = s + np.timedelta64(1, 's')
        result2 = np.timedelta64(1, 's') + s
        expected = Series([Timestamp('20130101 9:01:01'),
                           Timestamp('20130101 9:02:01')])
        assert_series_equal(result, expected)
        assert_series_equal(result2, expected)

        result = s + np.timedelta64(5, 'ms')
        result2 = np.timedelta64(5, 'ms') + s
        expected = Series([Timestamp('20130101 9:01:00.005'),
                           Timestamp('20130101 9:02:00.005')])
        assert_series_equal(result, expected)
        assert_series_equal(result2, expected)

    def test_dt64_series_add_tick_DateOffset(self):
        # GH 4532
        # operate with pd.offsets
        ser = Series([Timestamp('20130101 9:01'), Timestamp('20130101 9:02')])
        expected = Series([Timestamp('20130101 9:01:05'),
                           Timestamp('20130101 9:02:05')])

        result = ser + pd.offsets.Second(5)
        assert_series_equal(result, expected)

        result2 = pd.offsets.Second(5) + ser
        assert_series_equal(result2, expected)

    def test_dt64_series_sub_tick_DateOffset(self):
        # GH 4532
        # operate with pd.offsets
        ser = Series([Timestamp('20130101 9:01'), Timestamp('20130101 9:02')])
        expected = Series([Timestamp('20130101 9:00:55'),
                           Timestamp('20130101 9:01:55')])

        result = ser - pd.offsets.Second(5)
        assert_series_equal(result, expected)

        result2 = -pd.offsets.Second(5) + ser
        assert_series_equal(result2, expected)

        with pytest.raises(TypeError):
            pd.offsets.Second(5) - ser

    @pytest.mark.parametrize('cls_name', ['Day', 'Hour', 'Minute', 'Second',
                                          'Milli', 'Micro', 'Nano'])
    def test_dt64_series_with_tick_DateOffset_smoke(self, cls_name):
        # GH 4532
        # smoke tests for valid DateOffsets
        ser = Series([Timestamp('20130101 9:01'), Timestamp('20130101 9:02')])

        offset_cls = getattr(pd.offsets, cls_name)
        ser + offset_cls(5)
        offset_cls(5) + ser

    def test_dt64_series_add_mixed_tick_DateOffset(self):
        # GH 4532
        # operate with pd.offsets
        s = Series([Timestamp('20130101 9:01'), Timestamp('20130101 9:02')])

        result = s + pd.offsets.Milli(5)
        result2 = pd.offsets.Milli(5) + s
        expected = Series([Timestamp('20130101 9:01:00.005'),
                           Timestamp('20130101 9:02:00.005')])
        assert_series_equal(result, expected)
        assert_series_equal(result2, expected)

        result = s + pd.offsets.Minute(5) + pd.offsets.Milli(5)
        expected = Series([Timestamp('20130101 9:06:00.005'),
                           Timestamp('20130101 9:07:00.005')])
        assert_series_equal(result, expected)

    def test_dt64_series_sub_NaT(self):
        # GH#18808
        dti = pd.DatetimeIndex([pd.NaT, pd.Timestamp('19900315')])
        ser = pd.Series(dti)
        res = ser - pd.NaT
        expected = pd.Series([pd.NaT, pd.NaT], dtype='timedelta64[ns]')
        tm.assert_series_equal(res, expected)

        dti_tz = dti.tz_localize('Asia/Tokyo')
        ser_tz = pd.Series(dti_tz)
        res = ser_tz - pd.NaT
        expected = pd.Series([pd.NaT, pd.NaT], dtype='timedelta64[ns]')
        tm.assert_series_equal(res, expected)

    def test_datetime64_ops_nat(self):
        # GH 11349
        datetime_series = Series([NaT, Timestamp('19900315')])
        nat_series_dtype_timestamp = Series([NaT, NaT], dtype='datetime64[ns]')
        single_nat_dtype_datetime = Series([NaT], dtype='datetime64[ns]')

        # subtraction
        assert_series_equal(-NaT + datetime_series, nat_series_dtype_timestamp)
        with pytest.raises(TypeError):
            -single_nat_dtype_datetime + datetime_series

        assert_series_equal(-NaT + nat_series_dtype_timestamp,
                            nat_series_dtype_timestamp)
        with pytest.raises(TypeError):
            -single_nat_dtype_datetime + nat_series_dtype_timestamp

        # addition
        assert_series_equal(nat_series_dtype_timestamp + NaT,
                            nat_series_dtype_timestamp)
        assert_series_equal(NaT + nat_series_dtype_timestamp,
                            nat_series_dtype_timestamp)

        assert_series_equal(nat_series_dtype_timestamp + NaT,
                            nat_series_dtype_timestamp)
        assert_series_equal(NaT + nat_series_dtype_timestamp,
                            nat_series_dtype_timestamp)

    @pytest.mark.parametrize('dt64_series', [
        Series([Timestamp('19900315'), Timestamp('19900315')]),
        Series([NaT, Timestamp('19900315')]),
        Series([NaT, NaT], dtype='datetime64[ns]')])
    @pytest.mark.parametrize('one', [1, 1.0, np.array(1)])
    def test_dt64_mul_div_numeric_invalid(self, one, dt64_series):
        # multiplication
        with pytest.raises(TypeError):
            dt64_series * one
        with pytest.raises(TypeError):
            one * dt64_series

        # division
        with pytest.raises(TypeError):
            dt64_series / one
        with pytest.raises(TypeError):
            one / dt64_series

    def test_dt64_series_arith_overflow(self):
        # GH#12534, fixed by #19024
        dt = pd.Timestamp('1700-01-31')
        td = pd.Timedelta('20000 Days')
        dti = pd.date_range('1949-09-30', freq='100Y', periods=4)
        ser = pd.Series(dti)
        with pytest.raises(OverflowError):
            ser - dt
        with pytest.raises(OverflowError):
            dt - ser
        with pytest.raises(OverflowError):
            ser + td
        with pytest.raises(OverflowError):
            td + ser

        ser.iloc[-1] = pd.NaT
        expected = pd.Series(['2004-10-03', '2104-10-04', '2204-10-04', 'NaT'],
                             dtype='datetime64[ns]')
        res = ser + td
        tm.assert_series_equal(res, expected)
        res = td + ser
        tm.assert_series_equal(res, expected)

        ser.iloc[1:] = pd.NaT
        expected = pd.Series(['91279 Days', 'NaT', 'NaT', 'NaT'],
                             dtype='timedelta64[ns]')
        res = ser - dt
        tm.assert_series_equal(res, expected)
        res = dt - ser
        tm.assert_series_equal(res, -expected)

    @pytest.mark.parametrize('op', ['__add__', '__radd__',
                                    '__sub__', '__rsub__'])
    @pytest.mark.parametrize('tz', [None, 'Asia/Tokyo'])
    def test_dt64_series_add_intlike(self, tz, op):
        # GH#19123
        dti = pd.DatetimeIndex(['2016-01-02', '2016-02-03', 'NaT'], tz=tz)
        ser = Series(dti)

        other = Series([20, 30, 40], dtype='uint8')

        pytest.raises(TypeError, getattr(ser, op), 1)

        pytest.raises(TypeError, getattr(ser, op), other)

        pytest.raises(TypeError, getattr(ser, op), other.values)

        pytest.raises(TypeError, getattr(ser, op), pd.Index(other))


class TestSeriesOperators(TestData):
    @pytest.mark.parametrize(
        'ts',
        [
            (lambda x: x, lambda x: x * 2, False),
            (lambda x: x, lambda x: x[::2], False),
            (lambda x: x, lambda x: 5, True),
            (lambda x: tm.makeFloatSeries(),
             lambda x: tm.makeFloatSeries(),
             True)
        ])
    @pytest.mark.parametrize('opname', ['add', 'sub', 'mul', 'floordiv',
                                        'truediv', 'div', 'pow'])
    def test_op_method(self, opname, ts):
        # check that Series.{opname} behaves like Series.__{opname}__,
        series = ts[0](self.ts)
        other = ts[1](self.ts)
        check_reverse = ts[2]

        if opname == 'div' and compat.PY3:
            pytest.skip('div test only for Py3')

        op = getattr(Series, opname)

        if op == 'div':
            alt = operator.truediv
        else:
            alt = getattr(operator, opname)

        result = op(series, other)
        expected = alt(series, other)
        assert_almost_equal(result, expected)
        if check_reverse:
            rop = getattr(Series, "r" + opname)
            result = rop(series, other)
            expected = alt(other, series)
            assert_almost_equal(result, expected)

    def test_neg(self):
        assert_series_equal(-self.series, -1 * self.series)

    def test_invert(self):
        assert_series_equal(-(self.series < 0), ~(self.series < 0))

    def test_operators(self):
        def _check_op(series, other, op, pos_only=False,
                      check_dtype=True):
            left = np.abs(series) if pos_only else series
            right = np.abs(other) if pos_only else other

            cython_or_numpy = op(left, right)
            python = left.combine(right, op)
            assert_series_equal(cython_or_numpy, python,
                                check_dtype=check_dtype)

        def check(series, other):
            simple_ops = ['add', 'sub', 'mul', 'truediv', 'floordiv', 'mod']

            for opname in simple_ops:
                _check_op(series, other, getattr(operator, opname))

            _check_op(series, other, operator.pow, pos_only=True)

            _check_op(series, other, lambda x, y: operator.add(y, x))
            _check_op(series, other, lambda x, y: operator.sub(y, x))
            _check_op(series, other, lambda x, y: operator.truediv(y, x))
            _check_op(series, other, lambda x, y: operator.floordiv(y, x))
            _check_op(series, other, lambda x, y: operator.mul(y, x))
            _check_op(series, other, lambda x, y: operator.pow(y, x),
                      pos_only=True)
            _check_op(series, other, lambda x, y: operator.mod(y, x))

        check(self.ts, self.ts * 2)
        check(self.ts, self.ts * 0)
        check(self.ts, self.ts[::2])
        check(self.ts, 5)

        def check_comparators(series, other, check_dtype=True):
            _check_op(series, other, operator.gt, check_dtype=check_dtype)
            _check_op(series, other, operator.ge, check_dtype=check_dtype)
            _check_op(series, other, operator.eq, check_dtype=check_dtype)
            _check_op(series, other, operator.lt, check_dtype=check_dtype)
            _check_op(series, other, operator.le, check_dtype=check_dtype)

        check_comparators(self.ts, 5)
        check_comparators(self.ts, self.ts + 1, check_dtype=False)

    def test_divmod(self):
        def check(series, other):
            results = divmod(series, other)
            if isinstance(other, Iterable) and len(series) != len(other):
                # if the lengths don't match, this is the test where we use
                # `self.ts[::2]`. Pad every other value in `other_np` with nan.
                other_np = []
                for n in other:
                    other_np.append(n)
                    other_np.append(np.nan)
            else:
                other_np = other
            other_np = np.asarray(other_np)
            with np.errstate(all='ignore'):
                expecteds = divmod(series.values, np.asarray(other_np))

            for result, expected in zip(results, expecteds):
                # check the values, name, and index separately
                assert_almost_equal(np.asarray(result), expected)

                assert result.name == series.name
                assert_index_equal(result.index, series.index)

        check(self.ts, self.ts * 2)
        check(self.ts, self.ts * 0)
        check(self.ts, self.ts[::2])
        check(self.ts, 5)

    def test_operators_empty_int_corner(self):
        s1 = Series([], [], dtype=np.int32)
        s2 = Series({'x': 0.})
        assert_series_equal(s1 * s2, Series([np.nan], index=['x']))

    @pytest.mark.parametrize("m", [1, 3, 10])
    @pytest.mark.parametrize("unit", ['D', 'h', 'm', 's', 'ms', 'us', 'ns'])
    def test_timedelta64_conversions(self, m, unit):

        startdate = Series(date_range('2013-01-01', '2013-01-03'))
        enddate = Series(date_range('2013-03-01', '2013-03-03'))

        s1 = enddate - startdate
        s1[2] = np.nan

        # op
        expected = s1.apply(lambda x: x / np.timedelta64(m, unit))
        result = s1 / np.timedelta64(m, unit)
        assert_series_equal(result, expected)

        # reverse op
        expected = s1.apply(
            lambda x: Timedelta(np.timedelta64(m, unit)) / x)
        result = np.timedelta64(m, unit) / s1
        assert_series_equal(result, expected)

    @pytest.mark.parametrize('op', [operator.add, operator.sub])
    def test_timedelta64_equal_timedelta_supported_ops(self, op):
        ser = Series([Timestamp('20130301'), Timestamp('20130228 23:00:00'),
                      Timestamp('20130228 22:00:00'),
                      Timestamp('20130228 21:00:00')])

        intervals = 'D', 'h', 'm', 's', 'us'

        # TODO: unused
        # npy16_mappings = {'D': 24 * 60 * 60 * 1000000,
        #                   'h': 60 * 60 * 1000000,
        #                   'm': 60 * 1000000,
        #                   's': 1000000,
        #                   'us': 1}

        def timedelta64(*args):
            return sum(starmap(np.timedelta64, zip(args, intervals)))

        for d, h, m, s, us in product(*([range(2)] * 5)):
            nptd = timedelta64(d, h, m, s, us)
            pytd = timedelta(days=d, hours=h, minutes=m, seconds=s,
                             microseconds=us)
            lhs = op(ser, nptd)
            rhs = op(ser, pytd)

            assert_series_equal(lhs, rhs)

    def test_ops_nat_mixed_datetime64_timedelta64(self):
        # GH 11349
        timedelta_series = Series([NaT, Timedelta('1s')])
        datetime_series = Series([NaT, Timestamp('19900315')])
        nat_series_dtype_timedelta = Series([NaT, NaT],
                                            dtype='timedelta64[ns]')
        nat_series_dtype_timestamp = Series([NaT, NaT], dtype='datetime64[ns]')
        single_nat_dtype_datetime = Series([NaT], dtype='datetime64[ns]')
        single_nat_dtype_timedelta = Series([NaT], dtype='timedelta64[ns]')

        # subtraction
        assert_series_equal(datetime_series - single_nat_dtype_datetime,
                            nat_series_dtype_timedelta)

        assert_series_equal(datetime_series - single_nat_dtype_timedelta,
                            nat_series_dtype_timestamp)
        assert_series_equal(-single_nat_dtype_timedelta + datetime_series,
                            nat_series_dtype_timestamp)

        # without a Series wrapping the NaT, it is ambiguous
        # whether it is a datetime64 or timedelta64
        # defaults to interpreting it as timedelta64
        assert_series_equal(nat_series_dtype_timestamp -
                            single_nat_dtype_datetime,
                            nat_series_dtype_timedelta)

        assert_series_equal(nat_series_dtype_timestamp -
                            single_nat_dtype_timedelta,
                            nat_series_dtype_timestamp)
        assert_series_equal(-single_nat_dtype_timedelta +
                            nat_series_dtype_timestamp,
                            nat_series_dtype_timestamp)

        with pytest.raises(TypeError):
            timedelta_series - single_nat_dtype_datetime

        # addition
        assert_series_equal(nat_series_dtype_timestamp +
                            single_nat_dtype_timedelta,
                            nat_series_dtype_timestamp)
        assert_series_equal(single_nat_dtype_timedelta +
                            nat_series_dtype_timestamp,
                            nat_series_dtype_timestamp)

        assert_series_equal(nat_series_dtype_timestamp +
                            single_nat_dtype_timedelta,
                            nat_series_dtype_timestamp)
        assert_series_equal(single_nat_dtype_timedelta +
                            nat_series_dtype_timestamp,
                            nat_series_dtype_timestamp)

        assert_series_equal(nat_series_dtype_timedelta +
                            single_nat_dtype_datetime,
                            nat_series_dtype_timestamp)
        assert_series_equal(single_nat_dtype_datetime +
                            nat_series_dtype_timedelta,
                            nat_series_dtype_timestamp)

    def test_ops_datetimelike_align(self):
        # GH 7500
        # datetimelike ops need to align
        dt = Series(date_range('2012-1-1', periods=3, freq='D'))
        dt.iloc[2] = np.nan
        dt2 = dt[::-1]

        expected = Series([timedelta(0), timedelta(0), pd.NaT])
        # name is reset
        result = dt2 - dt
        assert_series_equal(result, expected)

        expected = Series(expected, name=0)
        result = (dt2.to_frame() - dt.to_frame())[0]
        assert_series_equal(result, expected)

    def test_operators_bitwise(self):
        # GH 9016: support bitwise op for integer types
        index = list('bca')

        s_tft = Series([True, False, True], index=index)
        s_fff = Series([False, False, False], index=index)
        s_tff = Series([True, False, False], index=index)
        s_empty = Series([])

        # TODO: unused
        # s_0101 = Series([0, 1, 0, 1])

        s_0123 = Series(range(4), dtype='int64')
        s_3333 = Series([3] * 4)
        s_4444 = Series([4] * 4)

        res = s_tft & s_empty
        expected = s_fff
        assert_series_equal(res, expected)

        res = s_tft | s_empty
        expected = s_tft
        assert_series_equal(res, expected)

        res = s_0123 & s_3333
        expected = Series(range(4), dtype='int64')
        assert_series_equal(res, expected)

        res = s_0123 | s_4444
        expected = Series(range(4, 8), dtype='int64')
        assert_series_equal(res, expected)

        s_a0b1c0 = Series([1], list('b'))

        res = s_tft & s_a0b1c0
        expected = s_tff.reindex(list('abc'))
        assert_series_equal(res, expected)

        res = s_tft | s_a0b1c0
        expected = s_tft.reindex(list('abc'))
        assert_series_equal(res, expected)

        n0 = 0
        res = s_tft & n0
        expected = s_fff
        assert_series_equal(res, expected)

        res = s_0123 & n0
        expected = Series([0] * 4)
        assert_series_equal(res, expected)

        n1 = 1
        res = s_tft & n1
        expected = s_tft
        assert_series_equal(res, expected)

        res = s_0123 & n1
        expected = Series([0, 1, 0, 1])
        assert_series_equal(res, expected)

        s_1111 = Series([1] * 4, dtype='int8')
        res = s_0123 & s_1111
        expected = Series([0, 1, 0, 1], dtype='int64')
        assert_series_equal(res, expected)

        res = s_0123.astype(np.int16) | s_1111.astype(np.int32)
        expected = Series([1, 1, 3, 3], dtype='int32')
        assert_series_equal(res, expected)

        pytest.raises(TypeError, lambda: s_1111 & 'a')
        pytest.raises(TypeError, lambda: s_1111 & ['a', 'b', 'c', 'd'])
        pytest.raises(TypeError, lambda: s_0123 & np.NaN)
        pytest.raises(TypeError, lambda: s_0123 & 3.14)
        pytest.raises(TypeError, lambda: s_0123 & [0.1, 4, 3.14, 2])

        # s_0123 will be all false now because of reindexing like s_tft
        if compat.PY3:
            # unable to sort incompatible object via .union.
            exp = Series([False] * 7, index=['b', 'c', 'a', 0, 1, 2, 3])
            with tm.assert_produces_warning(RuntimeWarning):
                assert_series_equal(s_tft & s_0123, exp)
        else:
            exp = Series([False] * 7, index=[0, 1, 2, 3, 'a', 'b', 'c'])
            assert_series_equal(s_tft & s_0123, exp)

        # s_tft will be all false now because of reindexing like s_0123
        if compat.PY3:
            # unable to sort incompatible object via .union.
            exp = Series([False] * 7, index=[0, 1, 2, 3, 'b', 'c', 'a'])
            with tm.assert_produces_warning(RuntimeWarning):
                assert_series_equal(s_0123 & s_tft, exp)
        else:
            exp = Series([False] * 7, index=[0, 1, 2, 3, 'a', 'b', 'c'])
            assert_series_equal(s_0123 & s_tft, exp)

        assert_series_equal(s_0123 & False, Series([False] * 4))
        assert_series_equal(s_0123 ^ False, Series([False, True, True, True]))
        assert_series_equal(s_0123 & [False], Series([False] * 4))
        assert_series_equal(s_0123 & (False), Series([False] * 4))
        assert_series_equal(s_0123 & Series([False, np.NaN, False, False]),
                            Series([False] * 4))

        s_ftft = Series([False, True, False, True])
        assert_series_equal(s_0123 & Series([0.1, 4, -3.14, 2]), s_ftft)

        s_abNd = Series(['a', 'b', np.NaN, 'd'])
        res = s_0123 & s_abNd
        expected = s_ftft
        assert_series_equal(res, expected)

    def test_scalar_na_cmp_corners(self):
        s = Series([2, 3, 4, 5, 6, 7, 8, 9, 10])

        def tester(a, b):
            return a & b

        pytest.raises(TypeError, tester, s, datetime(2005, 1, 1))

        s = Series([2, 3, 4, 5, 6, 7, 8, 9, datetime(2005, 1, 1)])
        s[::2] = np.nan

        expected = Series(True, index=s.index)
        expected[::2] = False
        assert_series_equal(tester(s, list(s)), expected)

        d = DataFrame({'A': s})
        # TODO: Fix this exception - needs to be fixed! (see GH5035)
        # (previously this was a TypeError because series returned
        # NotImplemented

        # this is an alignment issue; these are equivalent
        # https://github.com/pandas-dev/pandas/issues/5284

        pytest.raises(ValueError, lambda: d.__and__(s, axis='columns'))
        pytest.raises(ValueError, tester, s, d)

        # this is wrong as its not a boolean result
        # result = d.__and__(s,axis='index')

    def test_operators_corner(self):
        series = self.ts

        empty = Series([], index=Index([]))

        result = series + empty
        assert np.isnan(result).all()

        result = empty + Series([], index=Index([]))
        assert len(result) == 0

        # TODO: this returned NotImplemented earlier, what to do?
        # deltas = Series([timedelta(1)] * 5, index=np.arange(5))
        # sub_deltas = deltas[::2]
        # deltas5 = deltas * 5
        # deltas = deltas + sub_deltas

        # float + int
        int_ts = self.ts.astype(int)[:-5]
        added = self.ts + int_ts
        expected = Series(self.ts.values[:-5] + int_ts.values,
                          index=self.ts.index[:-5], name='ts')
        tm.assert_series_equal(added[:-5], expected)

    @pytest.mark.parametrize('op', [operator.add, operator.sub, operator.mul,
                                    operator.truediv, operator.floordiv])
    def test_operators_reverse_object(self, op):
        # GH 56
        arr = Series(np.random.randn(10), index=np.arange(10), dtype=object)

        result = op(1., arr)
        expected = op(1., arr.astype(float))
        assert_series_equal(result.astype(float), expected)

    def test_operators_combine(self):
        def _check_fill(meth, op, a, b, fill_value=0):
            exp_index = a.index.union(b.index)
            a = a.reindex(exp_index)
            b = b.reindex(exp_index)

            amask = isna(a)
            bmask = isna(b)

            exp_values = []
            for i in range(len(exp_index)):
                with np.errstate(all='ignore'):
                    if amask[i]:
                        if bmask[i]:
                            exp_values.append(nan)
                            continue
                        exp_values.append(op(fill_value, b[i]))
                    elif bmask[i]:
                        if amask[i]:
                            exp_values.append(nan)
                            continue
                        exp_values.append(op(a[i], fill_value))
                    else:
                        exp_values.append(op(a[i], b[i]))

            result = meth(a, b, fill_value=fill_value)
            expected = Series(exp_values, exp_index)
            assert_series_equal(result, expected)

        a = Series([nan, 1., 2., 3., nan], index=np.arange(5))
        b = Series([nan, 1, nan, 3, nan, 4.], index=np.arange(6))

        pairings = []
        for op in ['add', 'sub', 'mul', 'pow', 'truediv', 'floordiv']:
            fv = 0
            lop = getattr(Series, op)
            lequiv = getattr(operator, op)
            rop = getattr(Series, 'r' + op)
            # bind op at definition time...
            requiv = lambda x, y, op=op: getattr(operator, op)(y, x)
            pairings.append((lop, lequiv, fv))
            pairings.append((rop, requiv, fv))

        if compat.PY3:
            pairings.append((Series.div, operator.truediv, 1))
            pairings.append((Series.rdiv, lambda x, y: operator.truediv(y, x),
                             1))
        else:
            pairings.append((Series.div, operator.div, 1))
            pairings.append((Series.rdiv, lambda x, y: operator.div(y, x), 1))

        for op, equiv_op, fv in pairings:
            result = op(a, b)
            exp = equiv_op(a, b)
            assert_series_equal(result, exp)
            _check_fill(op, equiv_op, a, b, fill_value=fv)
            # should accept axis=0 or axis='rows'
            op(a, b, axis=0)

    def test_operators_na_handling(self):
        from decimal import Decimal
        from datetime import date
        s = Series([Decimal('1.3'), Decimal('2.3')],
                   index=[date(2012, 1, 1), date(2012, 1, 2)])

        result = s + s.shift(1)
        result2 = s.shift(1) + s
        assert isna(result[0])
        assert isna(result2[0])

        s = Series(['foo', 'bar', 'baz', np.nan])
        result = 'prefix_' + s
        expected = Series(['prefix_foo', 'prefix_bar', 'prefix_baz', np.nan])
        assert_series_equal(result, expected)

        result = s + '_suffix'
        expected = Series(['foo_suffix', 'bar_suffix', 'baz_suffix', np.nan])
        assert_series_equal(result, expected)

    def test_datetime64_with_index(self):
        # arithmetic integer ops with an index
        ser = Series(np.random.randn(5))
        expected = ser - ser.index.to_series()
        result = ser - ser.index
        assert_series_equal(result, expected)

        # GH 4629
        # arithmetic datetime64 ops with an index
        ser = Series(date_range('20130101', periods=5),
                     index=date_range('20130101', periods=5))
        expected = ser - ser.index.to_series()
        result = ser - ser.index
        assert_series_equal(result, expected)

        with pytest.raises(TypeError):
            # GH#18850
            result = ser - ser.index.to_period()

        df = DataFrame(np.random.randn(5, 2),
                       index=date_range('20130101', periods=5))
        df['date'] = Timestamp('20130102')
        df['expected'] = df['date'] - df.index.to_series()
        df['result'] = df['date'] - df.index
        assert_series_equal(df['result'], df['expected'], check_names=False)

    def test_op_duplicate_index(self):
        # GH14227
        s1 = Series([1, 2], index=[1, 1])
        s2 = Series([10, 10], index=[1, 2])
        result = s1 + s2
        expected = pd.Series([11, 12, np.nan], index=[1, 1, 2])
        assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "test_input,error_type",
        [
            (pd.Series([]), ValueError),

            # For strings, or any Series with dtype 'O'
            (pd.Series(['foo', 'bar', 'baz']), TypeError),
            (pd.Series([(1,), (2,)]), TypeError),

            # For mixed data types
            (
                pd.Series(['foo', 'foo', 'bar', 'bar', None, np.nan, 'baz']),
                TypeError
            ),
        ]
    )
    def test_assert_idxminmax_raises(self, test_input, error_type):
        """
        Cases where ``Series.argmax`` and related should raise an exception
        """
        with pytest.raises(error_type):
            test_input.idxmin()
        with pytest.raises(error_type):
            test_input.idxmin(skipna=False)
        with pytest.raises(error_type):
            test_input.idxmax()
        with pytest.raises(error_type):
            test_input.idxmax(skipna=False)

    def test_idxminmax_with_inf(self):
        # For numeric data with NA and Inf (GH #13595)
        s = pd.Series([0, -np.inf, np.inf, np.nan])

        assert s.idxmin() == 1
        assert np.isnan(s.idxmin(skipna=False))

        assert s.idxmax() == 2
        assert np.isnan(s.idxmax(skipna=False))

        # Using old-style behavior that treats floating point nan, -inf, and
        # +inf as missing
        with pd.option_context('mode.use_inf_as_na', True):
            assert s.idxmin() == 0
            assert np.isnan(s.idxmin(skipna=False))
            assert s.idxmax() == 0
            np.isnan(s.idxmax(skipna=False))


class TestSeriesOperationsDataFrameCompat(object):
    def test_operators_frame(self):
        # rpow does not work with DataFrame
        ts = tm.makeTimeSeries()
        ts.name = 'ts'

        df = DataFrame({'A': ts})

        assert_series_equal(ts + ts, ts + df['A'],
                            check_names=False)
        assert_series_equal(ts ** ts, ts ** df['A'],
                            check_names=False)
        assert_series_equal(ts < ts, ts < df['A'],
                            check_names=False)
        assert_series_equal(ts / ts, ts / df['A'],
                            check_names=False)

    def test_series_frame_radd_bug(self):
        # GH#353
        vals = Series(tm.rands_array(5, 10))
        result = 'foo_' + vals
        expected = vals.map(lambda x: 'foo_' + x)
        assert_series_equal(result, expected)

        frame = DataFrame({'vals': vals})
        result = 'foo_' + frame
        expected = DataFrame({'vals': vals.map(lambda x: 'foo_' + x)})
        assert_frame_equal(result, expected)

        ts = tm.makeTimeSeries()
        ts.name = 'ts'

        # really raise this time
        with pytest.raises(TypeError):
            datetime.now() + ts

        with pytest.raises(TypeError):
            ts + datetime.now()

    def test_bool_ops_df_compat(self):
        # GH 1134
        s1 = pd.Series([True, False, True], index=list('ABC'), name='x')
        s2 = pd.Series([True, True, False], index=list('ABD'), name='x')

        exp = pd.Series([True, False, False, False],
                        index=list('ABCD'), name='x')
        assert_series_equal(s1 & s2, exp)
        assert_series_equal(s2 & s1, exp)

        # True | np.nan => True
        exp = pd.Series([True, True, True, False],
                        index=list('ABCD'), name='x')
        assert_series_equal(s1 | s2, exp)
        # np.nan | True => np.nan, filled with False
        exp = pd.Series([True, True, False, False],
                        index=list('ABCD'), name='x')
        assert_series_equal(s2 | s1, exp)

        # DataFrame doesn't fill nan with False
        exp = pd.DataFrame({'x': [True, False, np.nan, np.nan]},
                           index=list('ABCD'))
        assert_frame_equal(s1.to_frame() & s2.to_frame(), exp)
        assert_frame_equal(s2.to_frame() & s1.to_frame(), exp)

        exp = pd.DataFrame({'x': [True, True, np.nan, np.nan]},
                           index=list('ABCD'))
        assert_frame_equal(s1.to_frame() | s2.to_frame(), exp)
        assert_frame_equal(s2.to_frame() | s1.to_frame(), exp)

        # different length
        s3 = pd.Series([True, False, True], index=list('ABC'), name='x')
        s4 = pd.Series([True, True, True, True], index=list('ABCD'), name='x')

        exp = pd.Series([True, False, True, False],
                        index=list('ABCD'), name='x')
        assert_series_equal(s3 & s4, exp)
        assert_series_equal(s4 & s3, exp)

        # np.nan | True => np.nan, filled with False
        exp = pd.Series([True, True, True, False],
                        index=list('ABCD'), name='x')
        assert_series_equal(s3 | s4, exp)
        # True | np.nan => True
        exp = pd.Series([True, True, True, True],
                        index=list('ABCD'), name='x')
        assert_series_equal(s4 | s3, exp)

        exp = pd.DataFrame({'x': [True, False, True, np.nan]},
                           index=list('ABCD'))
        assert_frame_equal(s3.to_frame() & s4.to_frame(), exp)
        assert_frame_equal(s4.to_frame() & s3.to_frame(), exp)

        exp = pd.DataFrame({'x': [True, True, True, np.nan]},
                           index=list('ABCD'))
        assert_frame_equal(s3.to_frame() | s4.to_frame(), exp)
        assert_frame_equal(s4.to_frame() | s3.to_frame(), exp)

    def test_arith_ops_df_compat(self):
        # GH#1134
        s1 = pd.Series([1, 2, 3], index=list('ABC'), name='x')
        s2 = pd.Series([2, 2, 2], index=list('ABD'), name='x')

        exp = pd.Series([3.0, 4.0, np.nan, np.nan],
                        index=list('ABCD'), name='x')
        assert_series_equal(s1 + s2, exp)
        assert_series_equal(s2 + s1, exp)

        exp = pd.DataFrame({'x': [3.0, 4.0, np.nan, np.nan]},
                           index=list('ABCD'))
        assert_frame_equal(s1.to_frame() + s2.to_frame(), exp)
        assert_frame_equal(s2.to_frame() + s1.to_frame(), exp)

        # different length
        s3 = pd.Series([1, 2, 3], index=list('ABC'), name='x')
        s4 = pd.Series([2, 2, 2, 2], index=list('ABCD'), name='x')

        exp = pd.Series([3, 4, 5, np.nan],
                        index=list('ABCD'), name='x')
        assert_series_equal(s3 + s4, exp)
        assert_series_equal(s4 + s3, exp)

        exp = pd.DataFrame({'x': [3, 4, 5, np.nan]},
                           index=list('ABCD'))
        assert_frame_equal(s3.to_frame() + s4.to_frame(), exp)
        assert_frame_equal(s4.to_frame() + s3.to_frame(), exp)
