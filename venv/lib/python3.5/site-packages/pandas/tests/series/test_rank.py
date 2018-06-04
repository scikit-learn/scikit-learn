# -*- coding: utf-8 -*-
from pandas import compat, Timestamp

import pytest

from distutils.version import LooseVersion
from numpy import nan
import numpy as np

from pandas import Series, date_range, NaT
from pandas.api.types import CategoricalDtype

from pandas.compat import product
from pandas.util.testing import assert_series_equal
import pandas.util.testing as tm
from pandas.tests.series.common import TestData
from pandas._libs.tslib import iNaT
from pandas._libs.algos import Infinity, NegInfinity
from itertools import chain
import pandas.util._test_decorators as td


class TestSeriesRank(TestData):
    s = Series([1, 3, 4, 2, nan, 2, 1, 5, nan, 3])

    results = {
        'average': np.array([1.5, 5.5, 7.0, 3.5, nan,
                             3.5, 1.5, 8.0, nan, 5.5]),
        'min': np.array([1, 5, 7, 3, nan, 3, 1, 8, nan, 5]),
        'max': np.array([2, 6, 7, 4, nan, 4, 2, 8, nan, 6]),
        'first': np.array([1, 5, 7, 3, nan, 4, 2, 8, nan, 6]),
        'dense': np.array([1, 3, 4, 2, nan, 2, 1, 5, nan, 3]),
    }

    def test_rank(self):
        pytest.importorskip('scipy.stats.special')
        rankdata = pytest.importorskip('scipy.stats.rankdata')

        self.ts[::2] = np.nan
        self.ts[:10][::3] = 4.

        ranks = self.ts.rank()
        oranks = self.ts.astype('O').rank()

        assert_series_equal(ranks, oranks)

        mask = np.isnan(self.ts)
        filled = self.ts.fillna(np.inf)

        # rankdata returns a ndarray
        exp = Series(rankdata(filled), index=filled.index, name='ts')
        exp[mask] = np.nan

        tm.assert_series_equal(ranks, exp)

        iseries = Series(np.arange(5).repeat(2))

        iranks = iseries.rank()
        exp = iseries.astype(float).rank()
        assert_series_equal(iranks, exp)
        iseries = Series(np.arange(5)) + 1.0
        exp = iseries / 5.0
        iranks = iseries.rank(pct=True)

        assert_series_equal(iranks, exp)

        iseries = Series(np.repeat(1, 100))
        exp = Series(np.repeat(0.505, 100))
        iranks = iseries.rank(pct=True)
        assert_series_equal(iranks, exp)

        iseries[1] = np.nan
        exp = Series(np.repeat(50.0 / 99.0, 100))
        exp[1] = np.nan
        iranks = iseries.rank(pct=True)
        assert_series_equal(iranks, exp)

        iseries = Series(np.arange(5)) + 1.0
        iseries[4] = np.nan
        exp = iseries / 4.0
        iranks = iseries.rank(pct=True)
        assert_series_equal(iranks, exp)

        iseries = Series(np.repeat(np.nan, 100))
        exp = iseries.copy()
        iranks = iseries.rank(pct=True)
        assert_series_equal(iranks, exp)

        iseries = Series(np.arange(5)) + 1
        iseries[4] = np.nan
        exp = iseries / 4.0
        iranks = iseries.rank(pct=True)
        assert_series_equal(iranks, exp)

        rng = date_range('1/1/1990', periods=5)
        iseries = Series(np.arange(5), rng) + 1
        iseries.iloc[4] = np.nan
        exp = iseries / 4.0
        iranks = iseries.rank(pct=True)
        assert_series_equal(iranks, exp)

        iseries = Series([1e-50, 1e-100, 1e-20, 1e-2, 1e-20 + 1e-30, 1e-1])
        exp = Series([2, 1, 3, 5, 4, 6.0])
        iranks = iseries.rank()
        assert_series_equal(iranks, exp)

        # GH 5968
        iseries = Series(['3 day', '1 day 10m', '-2 day', NaT],
                         dtype='m8[ns]')
        exp = Series([3, 2, 1, np.nan])
        iranks = iseries.rank()
        assert_series_equal(iranks, exp)

        values = np.array(
            [-50, -1, -1e-20, -1e-25, -1e-50, 0, 1e-40, 1e-20, 1e-10, 2, 40
             ], dtype='float64')
        random_order = np.random.permutation(len(values))
        iseries = Series(values[random_order])
        exp = Series(random_order + 1.0, dtype='float64')
        iranks = iseries.rank()
        assert_series_equal(iranks, exp)

    def test_rank_categorical(self):
        # GH issue #15420 rank incorrectly orders ordered categories

        # Test ascending/descending ranking for ordered categoricals
        exp = Series([1., 2., 3., 4., 5., 6.])
        exp_desc = Series([6., 5., 4., 3., 2., 1.])
        ordered = Series(
            ['first', 'second', 'third', 'fourth', 'fifth', 'sixth']
        ).astype(CategoricalDtype(categories=['first', 'second', 'third',
                                              'fourth', 'fifth', 'sixth'],
                                  ordered=True))
        assert_series_equal(ordered.rank(), exp)
        assert_series_equal(ordered.rank(ascending=False), exp_desc)

        # Unordered categoricals should be ranked as objects
        unordered = Series(['first', 'second', 'third', 'fourth',
                            'fifth', 'sixth']).astype(
            CategoricalDtype(categories=['first', 'second', 'third',
                                         'fourth', 'fifth', 'sixth'],
                             ordered=False))
        exp_unordered = Series([2., 4., 6., 3., 1., 5.])
        res = unordered.rank()
        assert_series_equal(res, exp_unordered)

        unordered1 = Series(
            [1, 2, 3, 4, 5, 6],
        ).astype(CategoricalDtype([1, 2, 3, 4, 5, 6], False))
        exp_unordered1 = Series([1., 2., 3., 4., 5., 6.])
        res1 = unordered1.rank()
        assert_series_equal(res1, exp_unordered1)

        # Test na_option for rank data
        na_ser = Series(
            ['first', 'second', 'third', 'fourth', 'fifth', 'sixth', np.NaN]
        ).astype(CategoricalDtype(['first', 'second', 'third', 'fourth',
                                   'fifth', 'sixth', 'seventh'], True))

        exp_top = Series([2., 3., 4., 5., 6., 7., 1.])
        exp_bot = Series([1., 2., 3., 4., 5., 6., 7.])
        exp_keep = Series([1., 2., 3., 4., 5., 6., np.NaN])

        assert_series_equal(na_ser.rank(na_option='top'), exp_top)
        assert_series_equal(na_ser.rank(na_option='bottom'), exp_bot)
        assert_series_equal(na_ser.rank(na_option='keep'), exp_keep)

        # Test na_option for rank data with ascending False
        exp_top = Series([7., 6., 5., 4., 3., 2., 1.])
        exp_bot = Series([6., 5., 4., 3., 2., 1., 7.])
        exp_keep = Series([6., 5., 4., 3., 2., 1., np.NaN])

        assert_series_equal(
            na_ser.rank(na_option='top', ascending=False),
            exp_top
        )
        assert_series_equal(
            na_ser.rank(na_option='bottom', ascending=False),
            exp_bot
        )
        assert_series_equal(
            na_ser.rank(na_option='keep', ascending=False),
            exp_keep
        )

        # Test with pct=True
        na_ser = Series(['first', 'second', 'third', 'fourth', np.NaN]).astype(
            CategoricalDtype(['first', 'second', 'third', 'fourth'], True))
        exp_top = Series([0.4, 0.6, 0.8, 1., 0.2])
        exp_bot = Series([0.2, 0.4, 0.6, 0.8, 1.])
        exp_keep = Series([0.25, 0.5, 0.75, 1., np.NaN])

        assert_series_equal(na_ser.rank(na_option='top', pct=True), exp_top)
        assert_series_equal(na_ser.rank(na_option='bottom', pct=True), exp_bot)
        assert_series_equal(na_ser.rank(na_option='keep', pct=True), exp_keep)

    def test_rank_signature(self):
        s = Series([0, 1])
        s.rank(method='average')
        pytest.raises(ValueError, s.rank, 'average')

    @pytest.mark.parametrize('contents,dtype', [
        ([-np.inf, -50, -1, -1e-20, -1e-25, -1e-50, 0, 1e-40, 1e-20, 1e-10,
          2, 40, np.inf],
         'float64'),
        ([-np.inf, -50, -1, -1e-20, -1e-25, -1e-45, 0, 1e-40, 1e-20, 1e-10,
          2, 40, np.inf],
         'float32'),
        ([np.iinfo(np.uint8).min, 1, 2, 100, np.iinfo(np.uint8).max],
         'uint8'),
        pytest.param([np.iinfo(np.int64).min, -100, 0, 1, 9999, 100000,
                      1e10, np.iinfo(np.int64).max],
                     'int64',
                     marks=pytest.mark.xfail(
                         reason="iNaT is equivalent to minimum value of dtype"
                         "int64 pending issue #16674")),
        ([NegInfinity(), '1', 'A', 'BA', 'Ba', 'C', Infinity()],
         'object')
    ])
    def test_rank_inf(self, contents, dtype):
        dtype_na_map = {
            'float64': np.nan,
            'float32': np.nan,
            'int64': iNaT,
            'object': None
        }
        # Insert nans at random positions if underlying dtype has missing
        # value. Then adjust the expected order by adding nans accordingly
        # This is for testing whether rank calculation is affected
        # when values are interwined with nan values.
        values = np.array(contents, dtype=dtype)
        exp_order = np.array(range(len(values)), dtype='float64') + 1.0
        if dtype in dtype_na_map:
            na_value = dtype_na_map[dtype]
            nan_indices = np.random.choice(range(len(values)), 5)
            values = np.insert(values, nan_indices, na_value)
            exp_order = np.insert(exp_order, nan_indices, np.nan)
        # shuffle the testing array and expected results in the same way
        random_order = np.random.permutation(len(values))
        iseries = Series(values[random_order])
        exp = Series(exp_order[random_order], dtype='float64')
        iranks = iseries.rank()
        assert_series_equal(iranks, exp)

    def test_rank_tie_methods(self):
        s = self.s

        def _check(s, expected, method='average'):
            result = s.rank(method=method)
            tm.assert_series_equal(result, Series(expected))

        dtypes = [None, object]
        disabled = set([(object, 'first')])
        results = self.results

        for method, dtype in product(results, dtypes):
            if (dtype, method) in disabled:
                continue
            series = s if dtype is None else s.astype(dtype)
            _check(series, results[method], method=method)

    @td.skip_if_no_scipy
    @pytest.mark.parametrize('ascending', [True, False])
    @pytest.mark.parametrize('method', ['average', 'min', 'max', 'first',
                                        'dense'])
    @pytest.mark.parametrize('na_option', ['top', 'bottom', 'keep'])
    def test_rank_tie_methods_on_infs_nans(self, method, na_option, ascending):
        dtypes = [('object', None, Infinity(), NegInfinity()),
                  ('float64', np.nan, np.inf, -np.inf)]
        chunk = 3
        disabled = set([('object', 'first')])

        def _check(s, method, na_option, ascending):
            exp_ranks = {
                'average': ([2, 2, 2], [5, 5, 5], [8, 8, 8]),
                'min': ([1, 1, 1], [4, 4, 4], [7, 7, 7]),
                'max': ([3, 3, 3], [6, 6, 6], [9, 9, 9]),
                'first': ([1, 2, 3], [4, 5, 6], [7, 8, 9]),
                'dense': ([1, 1, 1], [2, 2, 2], [3, 3, 3])
            }
            ranks = exp_ranks[method]
            if na_option == 'top':
                order = [ranks[1], ranks[0], ranks[2]]
            elif na_option == 'bottom':
                order = [ranks[0], ranks[2], ranks[1]]
            else:
                order = [ranks[0], [np.nan] * chunk, ranks[1]]
            expected = order if ascending else order[::-1]
            expected = list(chain.from_iterable(expected))
            result = s.rank(method=method, na_option=na_option,
                            ascending=ascending)
            tm.assert_series_equal(result, Series(expected, dtype='float64'))

        for dtype, na_value, pos_inf, neg_inf in dtypes:
            in_arr = [neg_inf] * chunk + [na_value] * chunk + [pos_inf] * chunk
            iseries = Series(in_arr, dtype=dtype)
            if (dtype, method) in disabled:
                continue
            _check(iseries, method, na_option, ascending)

    def test_rank_desc_mix_nans_infs(self):
        # GH 19538
        # check descending ranking when mix nans and infs
        iseries = Series([1, np.nan, np.inf, -np.inf, 25])
        result = iseries.rank(ascending=False)
        exp = Series([3, np.nan, 1, 4, 2], dtype='float64')
        tm.assert_series_equal(result, exp)

    def test_rank_methods_series(self):
        pytest.importorskip('scipy.stats.special')
        rankdata = pytest.importorskip('scipy.stats.rankdata')
        import scipy

        xs = np.random.randn(9)
        xs = np.concatenate([xs[i:] for i in range(0, 9, 2)])  # add duplicates
        np.random.shuffle(xs)

        index = [chr(ord('a') + i) for i in range(len(xs))]

        for vals in [xs, xs + 1e6, xs * 1e-6]:
            ts = Series(vals, index=index)

            for m in ['average', 'min', 'max', 'first', 'dense']:
                result = ts.rank(method=m)
                sprank = rankdata(vals, m if m != 'first' else 'ordinal')
                expected = Series(sprank, index=index)

                if LooseVersion(scipy.__version__) >= LooseVersion('0.17.0'):
                    expected = expected.astype('float64')
                tm.assert_series_equal(result, expected)

    def test_rank_dense_method(self):
        dtypes = ['O', 'f8', 'i8']
        in_out = [([1], [1]),
                  ([2], [1]),
                  ([0], [1]),
                  ([2, 2], [1, 1]),
                  ([1, 2, 3], [1, 2, 3]),
                  ([4, 2, 1], [3, 2, 1],),
                  ([1, 1, 5, 5, 3], [1, 1, 3, 3, 2]),
                  ([-5, -4, -3, -2, -1], [1, 2, 3, 4, 5])]

        for ser, exp in in_out:
            for dtype in dtypes:
                s = Series(ser).astype(dtype)
                result = s.rank(method='dense')
                expected = Series(exp).astype(result.dtype)
                assert_series_equal(result, expected)

    def test_rank_descending(self):
        dtypes = ['O', 'f8', 'i8']

        for dtype, method in product(dtypes, self.results):
            if 'i' in dtype:
                s = self.s.dropna()
            else:
                s = self.s.astype(dtype)

            res = s.rank(ascending=False)
            expected = (s.max() - s).rank()
            assert_series_equal(res, expected)

            if method == 'first' and dtype == 'O':
                continue

            expected = (s.max() - s).rank(method=method)
            res2 = s.rank(method=method, ascending=False)
            assert_series_equal(res2, expected)

    def test_rank_int(self):
        s = self.s.dropna().astype('i8')

        for method, res in compat.iteritems(self.results):
            result = s.rank(method=method)
            expected = Series(res).dropna()
            expected.index = result.index
            assert_series_equal(result, expected)

    def test_rank_object_bug(self):
        # GH 13445

        # smoke tests
        Series([np.nan] * 32).astype(object).rank(ascending=True)
        Series([np.nan] * 32).astype(object).rank(ascending=False)

    def test_rank_modify_inplace(self):
        # GH 18521
        # Check rank does not mutate series
        s = Series([Timestamp('2017-01-05 10:20:27.569000'), NaT])
        expected = s.copy()

        s.rank()
        result = s
        assert_series_equal(result, expected)


# GH15630, pct should be on 100% basis when method='dense'

@pytest.mark.parametrize('dtype', ['O', 'f8', 'i8'])
@pytest.mark.parametrize('ser, exp', [
    ([1], [1.]),
    ([1, 2], [1. / 2, 2. / 2]),
    ([2, 2], [1., 1.]),
    ([1, 2, 3], [1. / 3, 2. / 3, 3. / 3]),
    ([1, 2, 2], [1. / 2, 2. / 2, 2. / 2]),
    ([4, 2, 1], [3. / 3, 2. / 3, 1. / 3],),
    ([1, 1, 5, 5, 3], [1. / 3, 1. / 3, 3. / 3, 3. / 3, 2. / 3]),
    ([1, 1, 3, 3, 5, 5], [1. / 3, 1. / 3, 2. / 3, 2. / 3, 3. / 3, 3. / 3]),
    ([-5, -4, -3, -2, -1], [1. / 5, 2. / 5, 3. / 5, 4. / 5, 5. / 5])])
def test_rank_dense_pct(dtype, ser, exp):
        s = Series(ser).astype(dtype)
        result = s.rank(method='dense', pct=True)
        expected = Series(exp).astype(result.dtype)
        assert_series_equal(result, expected)


@pytest.mark.parametrize('dtype', ['O', 'f8', 'i8'])
@pytest.mark.parametrize('ser, exp', [
    ([1], [1.]),
    ([1, 2], [1. / 2, 2. / 2]),
    ([2, 2], [1. / 2, 1. / 2]),
    ([1, 2, 3], [1. / 3, 2. / 3, 3. / 3]),
    ([1, 2, 2], [1. / 3, 2. / 3, 2. / 3]),
    ([4, 2, 1], [3. / 3, 2. / 3, 1. / 3],),
    ([1, 1, 5, 5, 3], [1. / 5, 1. / 5, 4. / 5, 4. / 5, 3. / 5]),
    ([1, 1, 3, 3, 5, 5], [1. / 6, 1. / 6, 3. / 6, 3. / 6, 5. / 6, 5. / 6]),
    ([-5, -4, -3, -2, -1], [1. / 5, 2. / 5, 3. / 5, 4. / 5, 5. / 5])])
def test_rank_min_pct(dtype, ser, exp):
        s = Series(ser).astype(dtype)
        result = s.rank(method='min', pct=True)
        expected = Series(exp).astype(result.dtype)
        assert_series_equal(result, expected)


@pytest.mark.parametrize('dtype', ['O', 'f8', 'i8'])
@pytest.mark.parametrize('ser, exp', [
    ([1], [1.]),
    ([1, 2], [1. / 2, 2. / 2]),
    ([2, 2], [1., 1.]),
    ([1, 2, 3], [1. / 3, 2. / 3, 3. / 3]),
    ([1, 2, 2], [1. / 3, 3. / 3, 3. / 3]),
    ([4, 2, 1], [3. / 3, 2. / 3, 1. / 3],),
    ([1, 1, 5, 5, 3], [2. / 5, 2. / 5, 5. / 5, 5. / 5, 3. / 5]),
    ([1, 1, 3, 3, 5, 5], [2. / 6, 2. / 6, 4. / 6, 4. / 6, 6. / 6, 6. / 6]),
    ([-5, -4, -3, -2, -1], [1. / 5, 2. / 5, 3. / 5, 4. / 5, 5. / 5])])
def test_rank_max_pct(dtype, ser, exp):
        s = Series(ser).astype(dtype)
        result = s.rank(method='max', pct=True)
        expected = Series(exp).astype(result.dtype)
        assert_series_equal(result, expected)


@pytest.mark.parametrize('dtype', ['O', 'f8', 'i8'])
@pytest.mark.parametrize('ser, exp', [
    ([1], [1.]),
    ([1, 2], [1. / 2, 2. / 2]),
    ([2, 2], [1.5 / 2, 1.5 / 2]),
    ([1, 2, 3], [1. / 3, 2. / 3, 3. / 3]),
    ([1, 2, 2], [1. / 3, 2.5 / 3, 2.5 / 3]),
    ([4, 2, 1], [3. / 3, 2. / 3, 1. / 3],),
    ([1, 1, 5, 5, 3], [1.5 / 5, 1.5 / 5, 4.5 / 5, 4.5 / 5, 3. / 5]),
    ([1, 1, 3, 3, 5, 5],
     [1.5 / 6, 1.5 / 6, 3.5 / 6, 3.5 / 6, 5.5 / 6, 5.5 / 6]),
    ([-5, -4, -3, -2, -1], [1. / 5, 2. / 5, 3. / 5, 4. / 5, 5. / 5])])
def test_rank_average_pct(dtype, ser, exp):
        s = Series(ser).astype(dtype)
        result = s.rank(method='average', pct=True)
        expected = Series(exp).astype(result.dtype)
        assert_series_equal(result, expected)


@pytest.mark.parametrize('dtype', ['f8', 'i8'])
@pytest.mark.parametrize('ser, exp', [
    ([1], [1.]),
    ([1, 2], [1. / 2, 2. / 2]),
    ([2, 2], [1. / 2, 2. / 2.]),
    ([1, 2, 3], [1. / 3, 2. / 3, 3. / 3]),
    ([1, 2, 2], [1. / 3, 2. / 3, 3. / 3]),
    ([4, 2, 1], [3. / 3, 2. / 3, 1. / 3],),
    ([1, 1, 5, 5, 3], [1. / 5, 2. / 5, 4. / 5, 5. / 5, 3. / 5]),
    ([1, 1, 3, 3, 5, 5], [1. / 6, 2. / 6, 3. / 6, 4. / 6, 5. / 6, 6. / 6]),
    ([-5, -4, -3, -2, -1], [1. / 5, 2. / 5, 3. / 5, 4. / 5, 5. / 5])])
def test_rank_first_pct(dtype, ser, exp):
        s = Series(ser).astype(dtype)
        result = s.rank(method='first', pct=True)
        expected = Series(exp).astype(result.dtype)
        assert_series_equal(result, expected)
