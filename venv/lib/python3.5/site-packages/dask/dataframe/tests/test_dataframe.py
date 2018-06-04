import sys
import textwrap
from distutils.version import LooseVersion
from itertools import product
from operator import add

import pandas as pd
import pandas.util.testing as tm
import numpy as np
import pytest

import dask
import dask.array as da
import dask.dataframe as dd
from dask.base import compute_as_if_collection
from dask.utils import put_lines, M

from dask.dataframe.core import repartition_divisions, aca, _concat, Scalar
from dask.dataframe import methods
from dask.dataframe.utils import (assert_eq, make_meta, assert_max_deps,
                                  PANDAS_VERSION)

if PANDAS_VERSION >= '0.20.0':
    from pandas.io.formats import format as pandas_format
else:
    from pandas.formats import format as pandas_format


dsk = {('x', 0): pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]},
                              index=[0, 1, 3]),
       ('x', 1): pd.DataFrame({'a': [4, 5, 6], 'b': [3, 2, 1]},
                              index=[5, 6, 8]),
       ('x', 2): pd.DataFrame({'a': [7, 8, 9], 'b': [0, 0, 0]},
                              index=[9, 9, 9])}
meta = make_meta({'a': 'i8', 'b': 'i8'}, index=pd.Index([], 'i8'))
d = dd.DataFrame(dsk, 'x', meta, [0, 5, 9, 9])
full = d.compute()


def test_Dataframe():
    expected = pd.Series([2, 3, 4, 5, 6, 7, 8, 9, 10],
                         index=[0, 1, 3, 5, 6, 8, 9, 9, 9],
                         name='a')

    assert_eq(d['a'] + 1, expected)

    tm.assert_index_equal(d.columns, pd.Index(['a', 'b']))

    assert_eq(d[d['b'] > 2], full[full['b'] > 2])
    assert_eq(d[['a', 'b']], full[['a', 'b']])
    assert_eq(d.a, full.a)
    assert d.b.mean().compute() == full.b.mean()
    assert np.allclose(d.b.var().compute(), full.b.var())
    assert np.allclose(d.b.std().compute(), full.b.std())

    assert d.index._name == d.index._name  # this is deterministic

    assert repr(d)


def test_head_tail():
    assert_eq(d.head(2), full.head(2))
    assert_eq(d.head(3), full.head(3))
    assert_eq(d.head(2), dsk[('x', 0)].head(2))
    assert_eq(d['a'].head(2), full['a'].head(2))
    assert_eq(d['a'].head(3), full['a'].head(3))
    assert_eq(d['a'].head(2), dsk[('x', 0)]['a'].head(2))
    assert (sorted(d.head(2, compute=False).dask) ==
            sorted(d.head(2, compute=False).dask))
    assert (sorted(d.head(2, compute=False).dask) !=
            sorted(d.head(3, compute=False).dask))

    assert_eq(d.tail(2), full.tail(2))
    assert_eq(d.tail(3), full.tail(3))
    assert_eq(d.tail(2), dsk[('x', 2)].tail(2))
    assert_eq(d['a'].tail(2), full['a'].tail(2))
    assert_eq(d['a'].tail(3), full['a'].tail(3))
    assert_eq(d['a'].tail(2), dsk[('x', 2)]['a'].tail(2))
    assert (sorted(d.tail(2, compute=False).dask) ==
            sorted(d.tail(2, compute=False).dask))
    assert (sorted(d.tail(2, compute=False).dask) !=
            sorted(d.tail(3, compute=False).dask))


def test_head_npartitions():
    assert_eq(d.head(5, npartitions=2), full.head(5))
    assert_eq(d.head(5, npartitions=2, compute=False), full.head(5))
    assert_eq(d.head(5, npartitions=-1), full.head(5))
    assert_eq(d.head(7, npartitions=-1), full.head(7))
    assert_eq(d.head(2, npartitions=-1), full.head(2))
    with pytest.raises(ValueError):
        d.head(2, npartitions=5)


@pytest.mark.skipif(sys.version_info[:2] == (3, 3),
                    reason="Python3.3 uses pytest2.7.2, w/o warns method")
def test_head_npartitions_warn():
    with pytest.warns(None):
        d.head(100)

    with pytest.warns(None):
        d.head(7)

    with pytest.warns(None):
        d.head(7, npartitions=2)


def test_index_head():
    assert_eq(d.index.head(2), full.index[:2])
    assert_eq(d.index.head(3), full.index[:3])


def test_Series():
    assert isinstance(d.a, dd.Series)
    assert isinstance(d.a + 1, dd.Series)
    assert_eq((d + 1), full + 1)


def test_Index():
    for case in [pd.DataFrame(np.random.randn(10, 5), index=list('abcdefghij')),
                 pd.DataFrame(np.random.randn(10, 5),
                              index=pd.date_range('2011-01-01', freq='D',
                                                  periods=10))]:
        ddf = dd.from_pandas(case, 3)
        assert_eq(ddf.index, case.index)
        pytest.raises(AttributeError, lambda: ddf.index.index)


def test_Scalar():
    val = np.int64(1)
    s = Scalar({('a', 0): val}, 'a', 'i8')
    assert hasattr(s, 'dtype')
    assert 'dtype' in dir(s)
    assert_eq(s, val)
    assert repr(s) == "dd.Scalar<a, dtype=int64>"

    val = pd.Timestamp('2001-01-01')
    s = Scalar({('a', 0): val}, 'a', val)
    assert not hasattr(s, 'dtype')
    assert 'dtype' not in dir(s)
    assert_eq(s, val)
    assert repr(s) == "dd.Scalar<a, type=Timestamp>"


def test_attributes():
    assert 'a' in dir(d)
    assert 'foo' not in dir(d)
    pytest.raises(AttributeError, lambda: d.foo)

    df = dd.from_pandas(pd.DataFrame({'a b c': [1, 2, 3]}), npartitions=2)
    assert 'a b c' not in dir(df)
    df = dd.from_pandas(pd.DataFrame({'a': [1, 2], 5: [1, 2]}), npartitions=2)
    assert 'a' in dir(df)
    assert 5 not in dir(df)

    df = dd.from_pandas(tm.makeTimeDataFrame(), npartitions=3)
    pytest.raises(AttributeError, lambda: df.foo)


def test_column_names():
    tm.assert_index_equal(d.columns, pd.Index(['a', 'b']))
    tm.assert_index_equal(d[['b', 'a']].columns, pd.Index(['b', 'a']))
    assert d['a'].name == 'a'
    assert (d['a'] + 1).name == 'a'
    assert (d['a'] + d['b']).name is None


def test_index_names():
    assert d.index.name is None

    idx = pd.Index([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], name='x')
    df = pd.DataFrame(np.random.randn(10, 5), idx)
    ddf = dd.from_pandas(df, 3)
    assert ddf.index.name == 'x'
    assert ddf.index.compute().name == 'x'


@pytest.mark.parametrize(
    'npartitions',
    [1, pytest.mark.xfail(2, reason='pandas join removes freq')]
)
def test_timezone_freq(npartitions):
    s_naive = pd.Series(pd.date_range('20130101', periods=10))
    s_aware = pd.Series(pd.date_range('20130101', periods=10, tz='US/Eastern'))
    pdf = pd.DataFrame({'tz': s_aware, 'notz': s_naive})
    ddf = dd.from_pandas(pdf, npartitions=npartitions)

    assert pdf.tz[0].freq == ddf.compute().tz[0].freq == ddf.tz.compute()[0].freq


def test_rename_columns():
    # GH 819
    df = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7],
                       'b': [7, 6, 5, 4, 3, 2, 1]})
    ddf = dd.from_pandas(df, 2)

    ddf.columns = ['x', 'y']
    df.columns = ['x', 'y']
    tm.assert_index_equal(ddf.columns, pd.Index(['x', 'y']))
    tm.assert_index_equal(ddf._meta.columns, pd.Index(['x', 'y']))
    assert_eq(ddf, df)

    msg = r"Length mismatch: Expected axis has 2 elements, new values have 4 elements"
    with pytest.raises(ValueError) as err:
        ddf.columns = [1, 2, 3, 4]
    assert msg in str(err.value)

    # Multi-index columns
    df = pd.DataFrame({('A', '0') : [1, 2, 2, 3], ('B', 1) : [1, 2, 3, 4]})
    ddf = dd.from_pandas(df, npartitions=2)

    df.columns = ['x', 'y']
    ddf.columns = ['x', 'y']
    tm.assert_index_equal(ddf.columns, pd.Index(['x', 'y']))
    tm.assert_index_equal(ddf._meta.columns, pd.Index(['x', 'y']))
    assert_eq(ddf, df)


def test_rename_series():
    # GH 819
    s = pd.Series([1, 2, 3, 4, 5, 6, 7], name='x')
    ds = dd.from_pandas(s, 2)

    s.name = 'renamed'
    ds.name = 'renamed'
    assert s.name == 'renamed'
    assert_eq(ds, s)

    ind = s.index
    dind = ds.index
    ind.name = 'renamed'
    dind.name = 'renamed'
    assert ind.name == 'renamed'
    assert_eq(dind, ind)


def test_rename_series_method():
    # Series name
    s = pd.Series([1, 2, 3, 4, 5, 6, 7], name='x')
    ds = dd.from_pandas(s, 2)

    assert_eq(ds.rename('y'), s.rename('y'))
    assert ds.name == 'x'  # no mutation
    assert_eq(ds.rename(), s.rename())

    ds.rename('z', inplace=True)
    s.rename('z', inplace=True)
    assert ds.name == 'z'
    assert_eq(ds, s)

    # Series index
    s = pd.Series(['a', 'b', 'c', 'd', 'e', 'f', 'g'], name='x')
    ds = dd.from_pandas(s, 2)

    for is_sorted in [True, False]:
        res = ds.rename(lambda x: x ** 2, sorted_index=is_sorted)
        assert_eq(res, s.rename(lambda x: x ** 2))
        assert res.known_divisions == is_sorted

        res = ds.rename(s, sorted_index=is_sorted)
        assert_eq(res, s.rename(s))
        assert res.known_divisions == is_sorted

    with pytest.raises(ValueError):
        ds.rename(lambda x: -x, sorted_index=True)
    assert_eq(ds.rename(lambda x: -x), s.rename(lambda x: -x))

    res = ds.rename(ds)
    assert_eq(res, s.rename(s))
    assert not res.known_divisions

    ds2 = ds.clear_divisions()
    res = ds2.rename(lambda x: x**2, sorted_index=True)
    assert_eq(res, s.rename(lambda x: x**2))
    assert not res.known_divisions

    res = ds.rename(lambda x: x**2, inplace=True, sorted_index=True)
    assert res is ds
    s.rename(lambda x: x**2, inplace=True)
    assert_eq(ds, s)


def test_describe():
    # prepare test case which approx quantiles will be the same as actuals
    s = pd.Series(list(range(20)) * 4)
    df = pd.DataFrame({'a': list(range(20)) * 4, 'b': list(range(4)) * 20})

    ds = dd.from_pandas(s, 4)
    ddf = dd.from_pandas(df, 4)

    assert_eq(s.describe(), ds.describe())
    assert_eq(df.describe(), ddf.describe())
    assert_eq(s.describe(), ds.describe(split_every=2))
    assert_eq(df.describe(), ddf.describe(split_every=2))

    assert ds.describe(split_every=2)._name != ds.describe()._name
    assert ddf.describe(split_every=2)._name != ddf.describe()._name

    # remove string columns
    df = pd.DataFrame({'a': list(range(20)) * 4, 'b': list(range(4)) * 20,
                       'c': list('abcd') * 20})
    ddf = dd.from_pandas(df, 4)
    assert_eq(df.describe(), ddf.describe())
    assert_eq(df.describe(), ddf.describe(split_every=2))


def test_describe_empty():
    # https://github.com/dask/dask/issues/2326
    ddf = dd.from_pandas(pd.DataFrame({"A": ['a', 'b']}), 2)
    with pytest.raises(ValueError) as rec:
        ddf.describe()
    assert 'DataFrame contains only non-numeric data.' in str(rec)

    with pytest.raises(ValueError) as rec:
        ddf.A.describe()
    assert 'Cannot compute ``describe`` on object dtype.' in str(rec)


def test_cumulative():
    df = pd.DataFrame(np.random.randn(100, 5), columns=list('abcde'))
    df_out = pd.DataFrame(np.random.randn(100, 5), columns=list('abcde'))

    ddf = dd.from_pandas(df, 5)
    ddf_out = dd.from_pandas(df_out, 5)

    assert_eq(ddf.cumsum(), df.cumsum())
    assert_eq(ddf.cumprod(), df.cumprod())
    assert_eq(ddf.cummin(), df.cummin())
    assert_eq(ddf.cummax(), df.cummax())

    assert_eq(ddf.cumsum(axis=1), df.cumsum(axis=1))
    assert_eq(ddf.cumprod(axis=1), df.cumprod(axis=1))
    assert_eq(ddf.cummin(axis=1), df.cummin(axis=1))
    assert_eq(ddf.cummax(axis=1), df.cummax(axis=1))

    # testing out parameter if out parameter supported
    if LooseVersion(np.__version__) >= '1.13.0':
        np.cumsum(ddf, out=ddf_out)
        assert_eq(ddf_out, df.cumsum())
        np.cumprod(ddf, out=ddf_out)
        assert_eq(ddf_out, df.cumprod())
        ddf.cummin(out=ddf_out)
        assert_eq(ddf_out, df.cummin())
        ddf.cummax(out=ddf_out)
        assert_eq(ddf_out, df.cummax())

        np.cumsum(ddf, out=ddf_out, axis=1)
        assert_eq(ddf_out, df.cumsum(axis=1))
        np.cumprod(ddf, out=ddf_out, axis=1)
        assert_eq(ddf_out, df.cumprod(axis=1))
        ddf.cummin(out=ddf_out, axis=1)
        assert_eq(ddf_out, df.cummin(axis=1))
        ddf.cummax(out=ddf_out, axis=1)
        assert_eq(ddf_out, df.cummax(axis=1))

    assert_eq(ddf.a.cumsum(), df.a.cumsum())
    assert_eq(ddf.a.cumprod(), df.a.cumprod())
    assert_eq(ddf.a.cummin(), df.a.cummin())
    assert_eq(ddf.a.cummax(), df.a.cummax())

    # With NaNs
    df = pd.DataFrame({'a': [1, 2, np.nan, 4, 5, 6, 7, 8],
                       'b': [1, 2, np.nan, np.nan, np.nan, 5, np.nan, np.nan],
                       'c': [np.nan] * 8})
    ddf = dd.from_pandas(df, 3)

    assert_eq(df.cumsum(), ddf.cumsum())
    assert_eq(df.cummin(), ddf.cummin())
    assert_eq(df.cummax(), ddf.cummax())
    assert_eq(df.cumprod(), ddf.cumprod())

    assert_eq(df.cumsum(skipna=False), ddf.cumsum(skipna=False))
    assert_eq(df.cummin(skipna=False), ddf.cummin(skipna=False))
    assert_eq(df.cummax(skipna=False), ddf.cummax(skipna=False))
    assert_eq(df.cumprod(skipna=False), ddf.cumprod(skipna=False))

    assert_eq(df.cumsum(axis=1), ddf.cumsum(axis=1))
    assert_eq(df.cummin(axis=1), ddf.cummin(axis=1))
    assert_eq(df.cummax(axis=1), ddf.cummax(axis=1))
    assert_eq(df.cumprod(axis=1), ddf.cumprod(axis=1))

    assert_eq(df.cumsum(axis=1, skipna=False), ddf.cumsum(axis=1, skipna=False))
    assert_eq(df.cummin(axis=1, skipna=False), ddf.cummin(axis=1, skipna=False))
    assert_eq(df.cummax(axis=1, skipna=False), ddf.cummax(axis=1, skipna=False))
    assert_eq(df.cumprod(axis=1, skipna=False), ddf.cumprod(axis=1, skipna=False))


def test_dropna():
    df = pd.DataFrame({'x': [np.nan, 2, 3, 4, np.nan, 6],
                       'y': [1, 2, np.nan, 4, np.nan, np.nan],
                       'z': [1, 2, 3, 4, np.nan, np.nan]},
                      index=[10, 20, 30, 40, 50, 60])
    ddf = dd.from_pandas(df, 3)

    assert_eq(ddf.x.dropna(), df.x.dropna())
    assert_eq(ddf.y.dropna(), df.y.dropna())
    assert_eq(ddf.z.dropna(), df.z.dropna())

    assert_eq(ddf.dropna(), df.dropna())
    assert_eq(ddf.dropna(how='all'), df.dropna(how='all'))
    assert_eq(ddf.dropna(subset=['x']), df.dropna(subset=['x']))
    assert_eq(ddf.dropna(subset=['y', 'z']), df.dropna(subset=['y', 'z']))
    assert_eq(ddf.dropna(subset=['y', 'z'], how='all'),
              df.dropna(subset=['y', 'z'], how='all'))


@pytest.mark.parametrize('lower, upper', [(2, 5), (2.5, 3.5)])
def test_clip(lower, upper):

    df = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7, 8, 9],
                       'b': [3, 5, 2, 5, 7, 2, 4, 2, 4]})
    ddf = dd.from_pandas(df, 3)

    s = pd.Series([1, 2, 3, 4, 5, 6, 7, 8, 9])
    ds = dd.from_pandas(s, 3)

    assert_eq(ddf.clip(lower=lower, upper=upper),
              df.clip(lower=lower, upper=upper))
    assert_eq(ddf.clip(lower=lower), df.clip(lower=lower))
    assert_eq(ddf.clip(upper=upper), df.clip(upper=upper))

    assert_eq(ds.clip(lower=lower, upper=upper),
              s.clip(lower=lower, upper=upper))
    assert_eq(ds.clip(lower=lower), s.clip(lower=lower))
    assert_eq(ds.clip(upper=upper), s.clip(upper=upper))

    assert_eq(ddf.clip_lower(lower), df.clip_lower(lower))
    assert_eq(ddf.clip_lower(upper), df.clip_lower(upper))
    assert_eq(ddf.clip_upper(lower), df.clip_upper(lower))
    assert_eq(ddf.clip_upper(upper), df.clip_upper(upper))

    assert_eq(ds.clip_lower(lower), s.clip_lower(lower))
    assert_eq(ds.clip_lower(upper), s.clip_lower(upper))
    assert_eq(ds.clip_upper(lower), s.clip_upper(lower))
    assert_eq(ds.clip_upper(upper), s.clip_upper(upper))


def test_squeeze():
    df = pd.DataFrame({'x': [1, 3, 6]})
    df2 = pd.DataFrame({'x':[0]})
    s = pd.Series({'test': 0, 'b': 100})

    ddf = dd.from_pandas(df, 3)
    ddf2 = dd.from_pandas(df2, 3)
    ds = dd.from_pandas(s, 2)

    assert_eq(df.squeeze(), ddf.squeeze())
    assert_eq(pd.Series([0], name='x'), ddf2.squeeze())
    assert_eq(ds.squeeze(), s.squeeze())

    with pytest.raises(NotImplementedError) as info:
        ddf.squeeze(axis=0)
    msg = "{0} does not support squeeze along axis 0".format(type(ddf))
    assert msg in str(info.value)

    with pytest.raises(ValueError) as info:
        ddf.squeeze(axis=2)
    msg = 'No axis {0} for object type {1}'.format(2, type(ddf))
    assert msg in str(info.value)

    with pytest.raises(ValueError) as info:
        ddf.squeeze(axis='test')
    msg = 'No axis test for object type {0}'.format(type(ddf))
    assert msg in str(info.value)


def test_where_mask():
    pdf1 = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7, 8, 9],
                         'b': [3, 5, 2, 5, 7, 2, 4, 2, 4]})
    ddf1 = dd.from_pandas(pdf1, 2)
    pdf2 = pd.DataFrame({'a': [True, False, True] * 3,
                         'b': [False, False, True] * 3})
    ddf2 = dd.from_pandas(pdf2, 2)

    # different index
    pdf3 = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7, 8, 9],
                         'b': [3, 5, 2, 5, 7, 2, 4, 2, 4]},
                        index=[0, 1, 2, 3, 4, 5, 6, 7, 8])
    ddf3 = dd.from_pandas(pdf3, 2)
    pdf4 = pd.DataFrame({'a': [True, False, True] * 3,
                         'b': [False, False, True] * 3},
                        index=[5, 6, 7, 8, 9, 10, 11, 12, 13])
    ddf4 = dd.from_pandas(pdf4, 2)

    # different columns
    pdf5 = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7, 8, 9],
                         'b': [9, 4, 2, 6, 2, 3, 1, 6, 2],
                         'c': [5, 6, 7, 8, 9, 10, 11, 12, 13]},
                        index=[0, 1, 2, 3, 4, 5, 6, 7, 8])
    ddf5 = dd.from_pandas(pdf5, 2)
    pdf6 = pd.DataFrame({'a': [True, False, True] * 3,
                         'b': [False, False, True] * 3,
                         'c': [False] * 9,
                         'd': [True] * 9},
                        index=[5, 6, 7, 8, 9, 10, 11, 12, 13])
    ddf6 = dd.from_pandas(pdf6, 2)

    cases = [(ddf1, ddf2, pdf1, pdf2),
             (ddf1.repartition([0, 3, 6, 8]), ddf2, pdf1, pdf2),
             (ddf1, ddf4, pdf3, pdf4),
             (ddf3.repartition([0, 4, 6, 8]), ddf4.repartition([5, 9, 10, 13]),
              pdf3, pdf4),
             (ddf5, ddf6, pdf5, pdf6),
             (ddf5.repartition([0, 4, 7, 8]), ddf6, pdf5, pdf6),

             # use pd.DataFrame as cond
             (ddf1, pdf2, pdf1, pdf2),
             (ddf1, pdf4, pdf3, pdf4),
             (ddf5, pdf6, pdf5, pdf6)]

    for ddf, ddcond, pdf, pdcond in cases:
        assert isinstance(ddf, dd.DataFrame)
        assert isinstance(ddcond, (dd.DataFrame, pd.DataFrame))
        assert isinstance(pdf, pd.DataFrame)
        assert isinstance(pdcond, pd.DataFrame)

        assert_eq(ddf.where(ddcond), pdf.where(pdcond))
        assert_eq(ddf.mask(ddcond), pdf.mask(pdcond))
        assert_eq(ddf.where(ddcond, -ddf), pdf.where(pdcond, -pdf))
        assert_eq(ddf.mask(ddcond, -ddf), pdf.mask(pdcond, -pdf))

        assert_eq(ddf.where(ddcond.a, -ddf), pdf.where(pdcond.a, -pdf))
        assert_eq(ddf.mask(ddcond.a, -ddf), pdf.mask(pdcond.a, -pdf))
        assert_eq(ddf.a.where(ddcond.a), pdf.a.where(pdcond.a))
        assert_eq(ddf.a.mask(ddcond.a), pdf.a.mask(pdcond.a))
        assert_eq(ddf.a.where(ddcond.a, -ddf.a), pdf.a.where(pdcond.a, -pdf.a))
        assert_eq(ddf.a.mask(ddcond.a, -ddf.a), pdf.a.mask(pdcond.a, -pdf.a))


def test_map_partitions_multi_argument():
    assert_eq(dd.map_partitions(lambda a, b: a + b, d.a, d.b),
              full.a + full.b)
    assert_eq(dd.map_partitions(lambda a, b, c: a + b + c, d.a, d.b, 1),
              full.a + full.b + 1)


def test_map_partitions():
    assert_eq(d.map_partitions(lambda df: df, meta=d), full)
    assert_eq(d.map_partitions(lambda df: df), full)
    result = d.map_partitions(lambda df: df.sum(axis=1))
    assert_eq(result, full.sum(axis=1))

    assert_eq(d.map_partitions(lambda df: 1), pd.Series([1, 1, 1], dtype=np.int64),
              check_divisions=False)
    x = Scalar({('x', 0): 1}, 'x', int)
    result = dd.map_partitions(lambda x: 2, x)
    assert result.dtype in (np.int32, np.int64) and result.compute() == 2
    result = dd.map_partitions(lambda x: 4.0, x)
    assert result.dtype == np.float64 and result.compute() == 4.0


def test_map_partitions_names():
    func = lambda x: x
    assert (sorted(dd.map_partitions(func, d, meta=d).dask) ==
            sorted(dd.map_partitions(func, d, meta=d).dask))
    assert (sorted(dd.map_partitions(lambda x: x, d, meta=d, token=1).dask) ==
            sorted(dd.map_partitions(lambda x: x, d, meta=d, token=1).dask))

    func = lambda x, y: x
    assert (sorted(dd.map_partitions(func, d, d, meta=d).dask) ==
            sorted(dd.map_partitions(func, d, d, meta=d).dask))


def test_map_partitions_column_info():
    df = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [5, 6, 7, 8]})
    a = dd.from_pandas(df, npartitions=2)

    b = dd.map_partitions(lambda x: x, a, meta=a)
    tm.assert_index_equal(b.columns, a.columns)
    assert_eq(df, b)

    b = dd.map_partitions(lambda x: x, a.x, meta=a.x)
    assert b.name == a.x.name
    assert_eq(df.x, b)

    b = dd.map_partitions(lambda x: x, a.x, meta=a.x)
    assert b.name == a.x.name
    assert_eq(df.x, b)

    b = dd.map_partitions(lambda df: df.x + df.y, a)
    assert isinstance(b, dd.Series)
    assert b.dtype == 'i8'

    b = dd.map_partitions(lambda df: df.x + 1, a, meta=('x', 'i8'))
    assert isinstance(b, dd.Series)
    assert b.name == 'x'
    assert b.dtype == 'i8'


def test_map_partitions_method_names():
    df = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [5, 6, 7, 8]})
    a = dd.from_pandas(df, npartitions=2)

    b = a.map_partitions(lambda x: x)
    assert isinstance(b, dd.DataFrame)
    tm.assert_index_equal(b.columns, a.columns)

    b = a.map_partitions(lambda df: df.x + 1)
    assert isinstance(b, dd.Series)
    assert b.dtype == 'i8'

    b = a.map_partitions(lambda df: df.x + 1, meta=('x', 'i8'))
    assert isinstance(b, dd.Series)
    assert b.name == 'x'
    assert b.dtype == 'i8'


def test_map_partitions_keeps_kwargs_readable():
    df = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [5, 6, 7, 8]})
    a = dd.from_pandas(df, npartitions=2)

    def f(s, x=1):
        return s + x

    b = a.x.map_partitions(f, x=5)

    # NOTE: we'd like to ensure that we keep the keyword arguments readable
    # in the dask graph
    assert "['x', 5]" in str(b.dask)
    assert_eq(df.x + 5, b)

    assert a.x.map_partitions(f, x=5)._name != a.x.map_partitions(f, x=6)._name


def test_metadata_inference_single_partition_aligned_args():
    # https://github.com/dask/dask/issues/3034
    # Previously broadcastable series functionality broke this

    df = pd.DataFrame({'x': [1, 2, 3, 4, 5]})
    ddf = dd.from_pandas(df, npartitions=1)

    def check(df, df_x):
        assert len(df) == len(df_x)
        assert len(df) > 0
        return df

    res = dd.map_partitions(check, ddf, ddf.x)
    assert_eq(res, ddf)


def test_drop_duplicates():
    res = d.drop_duplicates()
    res2 = d.drop_duplicates(split_every=2)
    sol = full.drop_duplicates()
    assert_eq(res, sol)
    assert_eq(res2, sol)
    assert res._name != res2._name

    res = d.a.drop_duplicates()
    res2 = d.a.drop_duplicates(split_every=2)
    sol = full.a.drop_duplicates()
    assert_eq(res, sol)
    assert_eq(res2, sol)
    assert res._name != res2._name

    res = d.index.drop_duplicates()
    res2 = d.index.drop_duplicates(split_every=2)
    sol = full.index.drop_duplicates()
    assert_eq(res, sol)
    assert_eq(res2, sol)
    assert res._name != res2._name

    with pytest.raises(NotImplementedError):
        d.drop_duplicates(keep=False)


def test_drop_duplicates_subset():
    df = pd.DataFrame({'x': [1, 2, 3, 1, 2, 3],
                       'y': ['a', 'a', 'b', 'b', 'c', 'c']})
    ddf = dd.from_pandas(df, npartitions=2)

    for kwarg in [{'keep': 'first'}, {'keep': 'last'}]:
        assert_eq(df.x.drop_duplicates(**kwarg),
                  ddf.x.drop_duplicates(**kwarg))
        for ss in [['x'], 'y', ['x', 'y']]:
            assert_eq(df.drop_duplicates(subset=ss, **kwarg),
                      ddf.drop_duplicates(subset=ss, **kwarg))


def test_get_partition():
    pdf = pd.DataFrame(np.random.randn(10, 5), columns=list('abcde'))
    ddf = dd.from_pandas(pdf, 3)
    assert ddf.divisions == (0, 4, 8, 9)

    # DataFrame
    div1 = ddf.get_partition(0)
    assert isinstance(div1, dd.DataFrame)
    assert_eq(div1, pdf.loc[0:3])
    div2 = ddf.get_partition(1)
    assert_eq(div2, pdf.loc[4:7])
    div3 = ddf.get_partition(2)
    assert_eq(div3, pdf.loc[8:9])
    assert len(div1) + len(div2) + len(div3) == len(pdf)

    # Series
    div1 = ddf.a.get_partition(0)
    assert isinstance(div1, dd.Series)
    assert_eq(div1, pdf.a.loc[0:3])
    div2 = ddf.a.get_partition(1)
    assert_eq(div2, pdf.a.loc[4:7])
    div3 = ddf.a.get_partition(2)
    assert_eq(div3, pdf.a.loc[8:9])
    assert len(div1) + len(div2) + len(div3) == len(pdf.a)

    with pytest.raises(ValueError):
        ddf.get_partition(-1)

    with pytest.raises(ValueError):
        ddf.get_partition(3)


def test_ndim():
    assert (d.ndim == 2)
    assert (d.a.ndim == 1)
    assert (d.index.ndim == 1)


def test_dtype():
    assert (d.dtypes == full.dtypes).all()


def test_value_counts():
    df = pd.DataFrame({'x': [1, 2, 1, 3, 3, 1, 4]})
    ddf = dd.from_pandas(df, npartitions=3)
    result = ddf.x.value_counts()
    expected = df.x.value_counts()
    assert_eq(result, expected)
    result2 = ddf.x.value_counts(split_every=2)
    assert_eq(result2, expected)
    assert result._name != result2._name


def test_unique():
    pdf = pd.DataFrame({'x': [1, 2, 1, 3, 3, 1, 4, 2, 3, 1],
                        'y': ['a', 'c', 'b', np.nan, 'c',
                              'b', 'a', 'd', np.nan, 'a']})
    ddf = dd.from_pandas(pdf, npartitions=3)
    assert_eq(ddf.x.unique(), pd.Series(pdf.x.unique(), name='x'))
    assert_eq(ddf.y.unique(), pd.Series(pdf.y.unique(), name='y'))

    assert_eq(ddf.x.unique(split_every=2),
              pd.Series(pdf.x.unique(), name='x'))
    assert_eq(ddf.y.unique(split_every=2),
              pd.Series(pdf.y.unique(), name='y'))
    assert ddf.x.unique(split_every=2)._name != ddf.x.unique()._name


def test_isin():
    # Series test
    assert_eq(d.a.isin([0, 1, 2]), full.a.isin([0, 1, 2]))
    assert_eq(d.a.isin(pd.Series([0, 1, 2])),
              full.a.isin(pd.Series([0, 1, 2])))

    # DataFrame test
    assert_eq(d.isin([0, 1, 2]), full.isin([0, 1, 2]))


def test_len():
    assert len(d) == len(full)
    assert len(d.a) == len(full.a)


def test_size():
    assert_eq(d.size, full.size)
    assert_eq(d.a.size, full.a.size)
    assert_eq(d.index.size, full.index.size)


def test_nbytes():
    assert_eq(d.a.nbytes, full.a.nbytes)
    assert_eq(d.index.nbytes, full.index.nbytes)


def test_quantile():
    # series / multiple
    result = d.b.quantile([.3, .7])
    exp = full.b.quantile([.3, .7])  # result may different
    assert len(result) == 2
    assert result.divisions == (.3, .7)
    assert_eq(result.index, exp.index)
    assert isinstance(result, dd.Series)

    result = result.compute()
    assert isinstance(result, pd.Series)
    assert result.iloc[0] == 0
    assert 5 < result.iloc[1] < 6

    # index
    s = pd.Series(np.arange(10), index=np.arange(10))
    ds = dd.from_pandas(s, 2)

    result = ds.index.quantile([.3, .7])
    exp = s.quantile([.3, .7])
    assert len(result) == 2
    assert result.divisions == (.3, .7)
    assert_eq(result.index, exp.index)
    assert isinstance(result, dd.Series)

    result = result.compute()
    assert isinstance(result, pd.Series)
    assert 1 < result.iloc[0] < 2
    assert 7 < result.iloc[1] < 8

    # series / single
    result = d.b.quantile(.5)
    exp = full.b.quantile(.5)  # result may different
    assert isinstance(result, dd.core.Scalar)
    result = result.compute()
    assert 4 < result < 6


def test_quantile_missing():
    df = pd.DataFrame({"A": [0, np.nan, 2]})
    ddf = dd.from_pandas(df, 2)
    expected = df.quantile()
    result = ddf.quantile()
    assert_eq(result, expected)

    expected = df.A.quantile()
    result = ddf.A.quantile()
    assert_eq(result, expected)


def test_empty_quantile():
    result = d.b.quantile([])
    exp = full.b.quantile([])
    assert result.divisions == (None, None)

    assert result.name == 'b'
    assert result.compute().name == 'b'
    assert_eq(result, exp)


def test_dataframe_quantile():
    # column X is for test column order and result division
    df = pd.DataFrame({'A': np.arange(20),
                       'X': np.arange(20, 40),
                       'B': np.arange(10, 30),
                       'C': ['a', 'b', 'c', 'd'] * 5},
                      columns=['A', 'X', 'B', 'C'])
    ddf = dd.from_pandas(df, 3)

    result = ddf.quantile()
    assert result.npartitions == 1
    assert result.divisions == ('A', 'X')

    result = result.compute()
    assert isinstance(result, pd.Series)
    assert result.name == 0.5
    tm.assert_index_equal(result.index, pd.Index(['A', 'X', 'B']))
    assert (result > pd.Series([16, 36, 26], index=['A', 'X', 'B'])).all()
    assert (result < pd.Series([17, 37, 27], index=['A', 'X', 'B'])).all()

    result = ddf.quantile([0.25, 0.75])
    assert result.npartitions == 1
    assert result.divisions == (0.25, 0.75)

    result = result.compute()
    assert isinstance(result, pd.DataFrame)
    tm.assert_index_equal(result.index, pd.Index([0.25, 0.75]))
    tm.assert_index_equal(result.columns, pd.Index(['A', 'X', 'B']))
    minexp = pd.DataFrame([[1, 21, 11], [17, 37, 27]],
                          index=[0.25, 0.75], columns=['A', 'X', 'B'])
    assert (result > minexp).all().all()
    maxexp = pd.DataFrame([[2, 22, 12], [18, 38, 28]],
                          index=[0.25, 0.75], columns=['A', 'X', 'B'])
    assert (result < maxexp).all().all()

    assert_eq(ddf.quantile(axis=1), df.quantile(axis=1))
    pytest.raises(ValueError, lambda: ddf.quantile([0.25, 0.75], axis=1))


def test_index():
    assert_eq(d.index, full.index)


def test_assign():
    d_unknown = dd.from_pandas(full, npartitions=3, sort=False)
    assert not d_unknown.known_divisions
    res = d.assign(c=1,
                   d='string',
                   e=d.a.sum(),
                   f=d.a + d.b,
                   g=lambda x: x.a + x.b,
                   dt=pd.Timestamp(2018, 2, 13))
    res_unknown = d_unknown.assign(c=1,
                                   d='string',
                                   e=d_unknown.a.sum(),
                                   f=d_unknown.a + d_unknown.b,
                                   g=lambda x: x.a + x.b,
                                   dt=pd.Timestamp(2018, 2, 13))
    sol = full.assign(c=1,
                      d='string',
                      e=full.a.sum(),
                      f=full.a + full.b,
                      g=lambda x: x.a + x.b,
                      dt=pd.Timestamp(2018, 2, 13))
    assert_eq(res, sol)
    assert_eq(res_unknown, sol)

    res = d.assign(c=full.a + 1)
    assert_eq(res, full.assign(c=full.a + 1))

    # divisions unknown won't work with pandas
    with pytest.raises(ValueError):
        d_unknown.assign(c=full.a + 1)

    # unsupported type
    with pytest.raises(TypeError):
        d.assign(c=list(range(9)))

    # Fails when assigning known divisions to unknown divisions
    with pytest.raises(ValueError):
        d_unknown.assign(foo=d.a)
    # Fails when assigning unknown divisions to known divisions
    with pytest.raises(ValueError):
        d.assign(foo=d_unknown.a)


def test_map():
    assert_eq(d.a.map(lambda x: x + 1), full.a.map(lambda x: x + 1))
    lk = dict((v, v + 1) for v in full.a.values)
    assert_eq(d.a.map(lk), full.a.map(lk))
    assert_eq(d.b.map(lk), full.b.map(lk))
    lk = pd.Series(lk)
    assert_eq(d.a.map(lk), full.a.map(lk))
    assert_eq(d.b.map(lk), full.b.map(lk))
    assert_eq(d.b.map(lk, meta=d.b), full.b.map(lk))
    assert_eq(d.b.map(lk, meta=('b', 'i8')), full.b.map(lk))
    pytest.raises(TypeError, lambda: d.a.map(d.b))


def test_concat():
    x = _concat([pd.DataFrame(columns=['a', 'b']),
                 pd.DataFrame(columns=['a', 'b'])])
    assert list(x.columns) == ['a', 'b']
    assert len(x) == 0


def test_args():
    e = d.assign(c=d.a + 1)
    f = type(e)(*e._args)
    assert_eq(e, f)
    assert_eq(d.a, type(d.a)(*d.a._args))
    assert_eq(d.a.sum(), type(d.a.sum())(*d.a.sum()._args))


def test_known_divisions():
    assert d.known_divisions
    df = dd.DataFrame(dsk, 'x', meta, divisions=[None, None, None])
    assert not df.known_divisions


def test_unknown_divisions():
    dsk = {('x', 0): pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]}),
           ('x', 1): pd.DataFrame({'a': [4, 5, 6], 'b': [3, 2, 1]}),
           ('x', 2): pd.DataFrame({'a': [7, 8, 9], 'b': [0, 0, 0]})}
    meta = make_meta({'a': 'i8', 'b': 'i8'})
    d = dd.DataFrame(dsk, 'x', meta, [None, None, None, None])
    full = d.compute(get=dask.get)

    assert_eq(d.a.sum(), full.a.sum())
    assert_eq(d.a + d.b + 1, full.a + full.b + 1)


@pytest.mark.parametrize('join', ['inner', 'outer', 'left', 'right'])
def test_align(join):
    df1a = pd.DataFrame({'A': np.random.randn(10),
                         'B': np.random.randn(10)},
                        index=[1, 12, 5, 6, 3, 9, 10, 4, 13, 11])

    df1b = pd.DataFrame({'A': np.random.randn(10),
                         'B': np.random.randn(10)},
                        index=[0, 3, 2, 10, 5, 6, 7, 8, 12, 13])
    ddf1a = dd.from_pandas(df1a, 3)
    ddf1b = dd.from_pandas(df1b, 3)

    # DataFrame
    res1, res2 = ddf1a.align(ddf1b, join=join)
    exp1, exp2 = df1a.align(df1b, join=join)
    assert assert_eq(res1, exp1)
    assert assert_eq(res2, exp2)

    # Series
    res1, res2 = ddf1a['A'].align(ddf1b['B'], join=join)
    exp1, exp2 = df1a['A'].align(df1b['B'], join=join)
    assert assert_eq(res1, exp1)
    assert assert_eq(res2, exp2)

    # DataFrame with fill_value
    res1, res2 = ddf1a.align(ddf1b, join=join, fill_value=1)
    exp1, exp2 = df1a.align(df1b, join=join, fill_value=1)
    assert assert_eq(res1, exp1)
    assert assert_eq(res2, exp2)

    # Series
    res1, res2 = ddf1a['A'].align(ddf1b['B'], join=join, fill_value=1)
    exp1, exp2 = df1a['A'].align(df1b['B'], join=join, fill_value=1)
    assert assert_eq(res1, exp1)
    assert assert_eq(res2, exp2)


@pytest.mark.parametrize('join', ['inner', 'outer', 'left', 'right'])
def test_align_axis(join):
    df1a = pd.DataFrame({'A': np.random.randn(10),
                         'B': np.random.randn(10),
                         'C': np.random.randn(10)},
                        index=[1, 12, 5, 6, 3, 9, 10, 4, 13, 11])

    df1b = pd.DataFrame({'B': np.random.randn(10),
                         'C': np.random.randn(10),
                         'D': np.random.randn(10)},
                        index=[0, 3, 2, 10, 5, 6, 7, 8, 12, 13])
    ddf1a = dd.from_pandas(df1a, 3)
    ddf1b = dd.from_pandas(df1b, 3)

    res1, res2 = ddf1a.align(ddf1b, join=join, axis=0)
    exp1, exp2 = df1a.align(df1b, join=join, axis=0)
    assert assert_eq(res1, exp1)
    assert assert_eq(res2, exp2)

    res1, res2 = ddf1a.align(ddf1b, join=join, axis=1)
    exp1, exp2 = df1a.align(df1b, join=join, axis=1)
    assert assert_eq(res1, exp1)
    assert assert_eq(res2, exp2)

    res1, res2 = ddf1a.align(ddf1b, join=join, axis='index')
    exp1, exp2 = df1a.align(df1b, join=join, axis='index')
    assert assert_eq(res1, exp1)
    assert assert_eq(res2, exp2)

    res1, res2 = ddf1a.align(ddf1b, join=join, axis='columns')
    exp1, exp2 = df1a.align(df1b, join=join, axis='columns')
    assert assert_eq(res1, exp1)
    assert assert_eq(res2, exp2)

    # invalid
    with pytest.raises(ValueError):
        ddf1a.align(ddf1b, join=join, axis='XXX')

    with pytest.raises(ValueError):
        ddf1a['A'].align(ddf1b['B'], join=join, axis=1)


def test_combine():
    df1 = pd.DataFrame({'A': np.random.choice([1, 2, np.nan], 100),
                        'B': np.random.choice(['a', 'b', np.nan], 100)})

    df2 = pd.DataFrame({'A': np.random.choice([1, 2, 3], 100),
                        'B': np.random.choice(['a', 'b', 'c'], 100)})
    ddf1 = dd.from_pandas(df1, 4)
    ddf2 = dd.from_pandas(df2, 5)

    first = lambda a, b: a

    # DataFrame
    for dda, ddb, a, b in [(ddf1, ddf2, df1, df2),
                           (ddf1.A, ddf2.A, df1.A, df2.A),
                           (ddf1.B, ddf2.B, df1.B, df2.B)]:
        for func, fill_value in [(add, None), (add, 100), (first, None)]:
            sol = a.combine(b, func, fill_value=fill_value)
            assert_eq(dda.combine(ddb, func, fill_value=fill_value), sol)
            assert_eq(dda.combine(b, func, fill_value=fill_value), sol)

    assert_eq(ddf1.combine(ddf2, add, overwrite=False),
              df1.combine(df2, add, overwrite=False))
    assert dda.combine(ddb, add)._name == dda.combine(ddb, add)._name


def test_combine_first():
    df1 = pd.DataFrame({'A': np.random.choice([1, 2, np.nan], 100),
                        'B': np.random.choice(['a', 'b', np.nan], 100)})

    df2 = pd.DataFrame({'A': np.random.choice([1, 2, 3], 100),
                        'B': np.random.choice(['a', 'b', 'c'], 100)})
    ddf1 = dd.from_pandas(df1, 4)
    ddf2 = dd.from_pandas(df2, 5)

    # DataFrame
    assert_eq(ddf1.combine_first(ddf2), df1.combine_first(df2))
    assert_eq(ddf1.combine_first(df2), df1.combine_first(df2))

    # Series
    assert_eq(ddf1.A.combine_first(ddf2.A), df1.A.combine_first(df2.A))
    assert_eq(ddf1.A.combine_first(df2.A), df1.A.combine_first(df2.A))

    assert_eq(ddf1.B.combine_first(ddf2.B), df1.B.combine_first(df2.B))
    assert_eq(ddf1.B.combine_first(df2.B), df1.B.combine_first(df2.B))


def test_dataframe_picklable():
    from pickle import loads, dumps
    cloudpickle = pytest.importorskip('cloudpickle')
    cp_dumps = cloudpickle.dumps

    d = tm.makeTimeDataFrame()
    df = dd.from_pandas(d, npartitions=3)
    df = df + 2

    # dataframe
    df2 = loads(dumps(df))
    assert_eq(df, df2)
    df2 = loads(cp_dumps(df))
    assert_eq(df, df2)

    # series
    a2 = loads(dumps(df.A))
    assert_eq(df.A, a2)
    a2 = loads(cp_dumps(df.A))
    assert_eq(df.A, a2)

    # index
    i2 = loads(dumps(df.index))
    assert_eq(df.index, i2)
    i2 = loads(cp_dumps(df.index))
    assert_eq(df.index, i2)

    # scalar
    # lambdas are present, so only test cloudpickle
    s = df.A.sum()
    s2 = loads(cp_dumps(s))
    assert_eq(s, s2)


def test_random_partitions():
    a, b = d.random_split([0.5, 0.5], 42)
    assert isinstance(a, dd.DataFrame)
    assert isinstance(b, dd.DataFrame)
    assert a._name != b._name

    assert len(a.compute()) + len(b.compute()) == len(full)
    a2, b2 = d.random_split([0.5, 0.5], 42)
    assert a2._name == a._name
    assert b2._name == b._name

    parts = d.random_split([0.4, 0.5, 0.1], 42)
    names = set([p._name for p in parts])
    names.update([a._name, b._name])
    assert len(names) == 5

    with pytest.raises(ValueError):
        d.random_split([0.4, 0.5], 42)


def test_series_round():
    ps = pd.Series([1.123, 2.123, 3.123, 1.234, 2.234, 3.234], name='a')
    s = dd.from_pandas(ps, npartitions=3)
    assert_eq(s.round(), ps.round())


@pytest.mark.slow
def test_repartition():
    def _check_split_data(orig, d):
        """Check data is split properly"""
        keys = [k for k in d.dask if k[0].startswith('repartition-split')]
        keys = sorted(keys)
        sp = pd.concat([compute_as_if_collection(dd.DataFrame, d.dask, k)
                        for k in keys])
        assert_eq(orig, sp)
        assert_eq(orig, d)

    df = pd.DataFrame({'x': [1, 2, 3, 4, 5, 6], 'y': list('abdabd')},
                      index=[10, 20, 30, 40, 50, 60])
    a = dd.from_pandas(df, 2)

    b = a.repartition(divisions=[10, 20, 50, 60])
    assert b.divisions == (10, 20, 50, 60)
    assert_eq(a, b)
    assert_eq(compute_as_if_collection(dd.DataFrame, b.dask, (b._name, 0)),
              df.iloc[:1])

    for div in [[20, 60], [10, 50], [1],   # first / last element mismatch
                [0, 60], [10, 70],   # do not allow to expand divisions by default
                [10, 50, 20, 60],    # not sorted
                [10, 10, 20, 60]]:   # not unique (last element can be duplicated)

        pytest.raises(ValueError, lambda: a.repartition(divisions=div))

    pdf = pd.DataFrame(np.random.randn(7, 5), columns=list('abxyz'))
    for p in range(1, 7):
        ddf = dd.from_pandas(pdf, p)
        assert_eq(ddf, pdf)
        for div in [[0, 6], [0, 6, 6], [0, 5, 6], [0, 4, 6, 6],
                    [0, 2, 6], [0, 2, 6, 6],
                    [0, 2, 3, 6, 6], [0, 1, 2, 3, 4, 5, 6, 6]]:
            rddf = ddf.repartition(divisions=div)
            _check_split_data(ddf, rddf)
            assert rddf.divisions == tuple(div)
            assert_eq(pdf, rddf)

            rds = ddf.x.repartition(divisions=div)
            _check_split_data(ddf.x, rds)
            assert rds.divisions == tuple(div)
            assert_eq(pdf.x, rds)

        # expand divisions
        for div in [[-5, 10], [-2, 3, 5, 6], [0, 4, 5, 9, 10]]:
            rddf = ddf.repartition(divisions=div, force=True)
            _check_split_data(ddf, rddf)
            assert rddf.divisions == tuple(div)
            assert_eq(pdf, rddf)

            rds = ddf.x.repartition(divisions=div, force=True)
            _check_split_data(ddf.x, rds)
            assert rds.divisions == tuple(div)
            assert_eq(pdf.x, rds)

    pdf = pd.DataFrame({'x': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                        'y': [9, 8, 7, 6, 5, 4, 3, 2, 1, 0]},
                       index=list('abcdefghij'))
    for p in range(1, 7):
        ddf = dd.from_pandas(pdf, p)
        assert_eq(ddf, pdf)
        for div in [list('aj'), list('ajj'), list('adj'),
                    list('abfj'), list('ahjj'), list('acdj'), list('adfij'),
                    list('abdefgij'), list('abcdefghij')]:
            rddf = ddf.repartition(divisions=div)
            _check_split_data(ddf, rddf)
            assert rddf.divisions == tuple(div)
            assert_eq(pdf, rddf)

            rds = ddf.x.repartition(divisions=div)
            _check_split_data(ddf.x, rds)
            assert rds.divisions == tuple(div)
            assert_eq(pdf.x, rds)

        # expand divisions
        for div in [list('Yadijm'), list('acmrxz'), list('Yajz')]:
            rddf = ddf.repartition(divisions=div, force=True)
            _check_split_data(ddf, rddf)
            assert rddf.divisions == tuple(div)
            assert_eq(pdf, rddf)

            rds = ddf.x.repartition(divisions=div, force=True)
            _check_split_data(ddf.x, rds)
            assert rds.divisions == tuple(div)
            assert_eq(pdf.x, rds)


def test_repartition_divisions():
    result = repartition_divisions([0, 6], [0, 6, 6], 'a', 'b', 'c')
    assert result == {('b', 0): (methods.boundary_slice, ('a', 0), 0, 6, False),
                      ('b', 1): (methods.boundary_slice, ('a', 0), 6, 6, True),
                      ('c', 0): ('b', 0),
                      ('c', 1): ('b', 1)}

    result = repartition_divisions([1, 3, 7], [1, 4, 6, 7], 'a', 'b', 'c')
    assert result == {('b', 0): (methods.boundary_slice, ('a', 0), 1, 3, False),
                      ('b', 1): (methods.boundary_slice, ('a', 1), 3, 4, False),
                      ('b', 2): (methods.boundary_slice, ('a', 1), 4, 6, False),
                      ('b', 3): (methods.boundary_slice, ('a', 1), 6, 7, True),
                      ('c', 0): (methods.concat, [('b', 0), ('b', 1)]),
                      ('c', 1): ('b', 2),
                      ('c', 2): ('b', 3)}


def test_repartition_on_pandas_dataframe():
    df = pd.DataFrame({'x': [1, 2, 3, 4, 5, 6], 'y': list('abdabd')},
                      index=[10, 20, 30, 40, 50, 60])
    ddf = dd.repartition(df, divisions=[10, 20, 50, 60])
    assert isinstance(ddf, dd.DataFrame)
    assert ddf.divisions == (10, 20, 50, 60)
    assert_eq(ddf, df)

    ddf = dd.repartition(df.y, divisions=[10, 20, 50, 60])
    assert isinstance(ddf, dd.Series)
    assert ddf.divisions == (10, 20, 50, 60)
    assert_eq(ddf, df.y)


@pytest.mark.parametrize('use_index', [True, False])
@pytest.mark.parametrize('n', [1, 2, 4, 5])
@pytest.mark.parametrize('k', [1, 2, 4, 5])
@pytest.mark.parametrize('dtype', [int, float, 'M8[ns]'])
@pytest.mark.parametrize('transform', [lambda df: df, lambda df: df.x])
def test_repartition_npartitions(use_index, n, k, dtype, transform):
    df = pd.DataFrame({'x': [1, 2, 3, 4, 5, 6] * 10,
                       'y': list('abdabd') * 10},
                      index=pd.Series([10, 20, 30, 40, 50, 60] * 10, dtype=dtype))
    df = transform(df)
    a = dd.from_pandas(df, npartitions=n, sort=use_index)
    b = a.repartition(npartitions=k)
    assert_eq(a, b)
    assert b.npartitions == k
    parts = dask.get(b.dask, b.__dask_keys__())
    assert all(map(len, parts))


def test_repartition_npartitions_same_limits():
    df = pd.DataFrame({'x': [1, 2, 3]},
                      index=[pd.Timestamp('2017-05-09 00:00:00.006000'),
                             pd.Timestamp('2017-05-09 02:45:00.017999'),
                             pd.Timestamp('2017-05-09 05:59:58.938999')])

    ddf = dd.from_pandas(df, npartitions=2)

    ddf.repartition(npartitions=10)


def test_repartition_object_index():
    df = pd.DataFrame({'x': [1, 2, 3, 4, 5, 6] * 10},
                      index=list('abdabd') * 10)
    a = dd.from_pandas(df, npartitions=5)
    b = a.repartition(npartitions=2)
    assert b.npartitions == 2
    assert_eq(b, df)

    b = a.repartition(npartitions=10)
    assert b.npartitions == 10
    assert_eq(b, df)
    assert not b.known_divisions


@pytest.mark.slow
@pytest.mark.parametrize('npartitions', [1, 20, 243])
@pytest.mark.parametrize('freq', ['1D', '7D', '28h', '1h'])
@pytest.mark.parametrize('end', ['2000-04-15', '2000-04-15 12:37:01', '2000-01-01 12:37:00'])
@pytest.mark.parametrize('start', ['2000-01-01', '2000-01-01 12:30:00', '2000-01-01 12:30:00'])
def test_repartition_freq(npartitions, freq, start, end):
    start = pd.Timestamp(start)
    end = pd.Timestamp(end)
    ind = pd.DatetimeIndex(start=start, end=end, freq='60s')
    df = pd.DataFrame({'x': np.arange(len(ind))}, index=ind)
    ddf = dd.from_pandas(df, npartitions=npartitions, name='x')

    ddf2 = ddf.repartition(freq=freq)
    assert_eq(ddf2, df)


def test_repartition_freq_divisions():
    df = pd.DataFrame({'x': np.random.random(10)},
                      index=pd.DatetimeIndex(np.random.random(10) * 100e9))
    ddf = dd.from_pandas(df, npartitions=3)

    ddf2 = ddf.repartition(freq='15s')
    for div in ddf2.divisions[1:-1]:
        assert div == div.round('15s')
    assert ddf2.divisions[0] == df.index.min()
    assert ddf2.divisions[-1] == df.index.max()
    assert_eq(ddf2, ddf2)


def test_repartition_freq_errors():
    df = pd.DataFrame({'x': [1, 2, 3]})
    ddf = dd.from_pandas(df, npartitions=1)
    with pytest.raises(TypeError) as info:
        ddf.repartition(freq='1s')

    assert 'only' in str(info.value)
    assert 'timeseries' in str(info.value)


def test_repartition_freq_month():
    ts = pd.date_range("2015-01-01 00:00", " 2015-05-01 23:50", freq="10min")
    df = pd.DataFrame(np.random.randint(0,100,size=(len(ts),4)),
                      columns=list('ABCD'), index=ts)
    ddf = dd.from_pandas(df,npartitions=1).repartition(freq='1M')

    assert_eq(df, ddf)
    assert 2 < ddf.npartitions <= 6


def test_embarrassingly_parallel_operations():
    df = pd.DataFrame({'x': [1, 2, 3, 4, None, 6], 'y': list('abdabd')},
                      index=[10, 20, 30, 40, 50, 60])
    a = dd.from_pandas(df, 2)

    assert_eq(a.x.astype('float32'), df.x.astype('float32'))
    assert a.x.astype('float32').compute().dtype == 'float32'

    assert_eq(a.x.dropna(), df.x.dropna())

    assert_eq(a.x.between(2, 4), df.x.between(2, 4))

    assert_eq(a.x.clip(2, 4), df.x.clip(2, 4))

    assert_eq(a.x.notnull(), df.x.notnull())
    assert_eq(a.x.isnull(), df.x.isnull())
    assert_eq(a.notnull(), df.notnull())
    assert_eq(a.isnull(), df.isnull())

    assert len(a.sample(0.5).compute()) < len(df)


def test_fillna():
    df = tm.makeMissingDataframe(0.8, 42)
    ddf = dd.from_pandas(df, npartitions=5, sort=False)

    assert_eq(ddf.fillna(100), df.fillna(100))
    assert_eq(ddf.A.fillna(100), df.A.fillna(100))

    assert_eq(ddf.fillna(method='pad'), df.fillna(method='pad'))
    assert_eq(ddf.A.fillna(method='pad'), df.A.fillna(method='pad'))

    assert_eq(ddf.fillna(method='bfill'), df.fillna(method='bfill'))
    assert_eq(ddf.A.fillna(method='bfill'), df.A.fillna(method='bfill'))

    assert_eq(ddf.fillna(method='pad', limit=2),
              df.fillna(method='pad', limit=2))
    assert_eq(ddf.A.fillna(method='pad', limit=2),
              df.A.fillna(method='pad', limit=2))

    assert_eq(ddf.fillna(method='bfill', limit=2),
              df.fillna(method='bfill', limit=2))
    assert_eq(ddf.A.fillna(method='bfill', limit=2),
              df.A.fillna(method='bfill', limit=2))

    assert_eq(ddf.fillna(100, axis=1), df.fillna(100, axis=1))
    assert_eq(ddf.fillna(method='pad', axis=1), df.fillna(method='pad', axis=1))
    assert_eq(ddf.fillna(method='pad', limit=2, axis=1),
              df.fillna(method='pad', limit=2, axis=1))

    pytest.raises(ValueError, lambda: ddf.A.fillna(0, axis=1))
    pytest.raises(NotImplementedError, lambda: ddf.fillna(0, limit=10))
    pytest.raises(NotImplementedError, lambda: ddf.fillna(0, limit=10, axis=1))

    df = tm.makeMissingDataframe(0.2, 42)
    ddf = dd.from_pandas(df, npartitions=5, sort=False)
    pytest.raises(ValueError, lambda: ddf.fillna(method='pad').compute())
    assert_eq(df.fillna(method='pad', limit=3),
              ddf.fillna(method='pad', limit=3))


def test_fillna_multi_dataframe():
    df = tm.makeMissingDataframe(0.8, 42)
    ddf = dd.from_pandas(df, npartitions=5, sort=False)

    assert_eq(ddf.A.fillna(ddf.B), df.A.fillna(df.B))
    assert_eq(ddf.B.fillna(ddf.A), df.B.fillna(df.A))


def test_ffill_bfill():
    df = tm.makeMissingDataframe(0.8, 42)
    ddf = dd.from_pandas(df, npartitions=5, sort=False)

    assert_eq(ddf.ffill(), df.ffill())
    assert_eq(ddf.bfill(), df.bfill())
    assert_eq(ddf.ffill(axis=1), df.ffill(axis=1))
    assert_eq(ddf.bfill(axis=1), df.bfill(axis=1))


def test_fillna_series_types():
    # https://github.com/dask/dask/issues/2809
    df = pd.DataFrame({"A": [1, np.nan, 3], "B": [1, np.nan, 3]})
    ddf = dd.from_pandas(df, npartitions=2)
    fill_value = pd.Series([1, 10], index=['A', 'C'])
    assert_eq(ddf.fillna(fill_value), df.fillna(fill_value))


def test_sample():
    df = pd.DataFrame({'x': [1, 2, 3, 4, None, 6], 'y': list('abdabd')},
                      index=[10, 20, 30, 40, 50, 60])
    a = dd.from_pandas(df, 2)

    b = a.sample(0.5)

    assert_eq(b, b)

    c = a.sample(0.5, random_state=1234)
    d = a.sample(0.5, random_state=1234)
    assert_eq(c, d)

    assert a.sample(0.5)._name != a.sample(0.5)._name


def test_sample_without_replacement():
    df = pd.DataFrame({'x': [1, 2, 3, 4, None, 6], 'y': list('abdabd')},
                      index=[10, 20, 30, 40, 50, 60])
    a = dd.from_pandas(df, 2)
    b = a.sample(0.7, replace=False)
    bb = b.index.compute()
    assert len(bb) == len(set(bb))


def test_datetime_accessor():
    df = pd.DataFrame({'x': [1, 2, 3, 4]})
    df['x'] = df.x.astype('M8[us]')

    a = dd.from_pandas(df, 2)

    assert 'date' in dir(a.x.dt)

    # pandas loses Series.name via datetime accessor
    # see https://github.com/pydata/pandas/issues/10712
    assert_eq(a.x.dt.date, df.x.dt.date, check_names=False)

    # to_pydatetime returns a numpy array in pandas, but a Series in dask
    assert_eq(a.x.dt.to_pydatetime(),
              pd.Series(df.x.dt.to_pydatetime(), index=df.index, dtype=object))

    assert set(a.x.dt.date.dask) == set(a.x.dt.date.dask)
    assert set(a.x.dt.to_pydatetime().dask) == set(a.x.dt.to_pydatetime().dask)


def test_str_accessor():
    df = pd.DataFrame({'x': ['abc', 'bcd', 'cdef', 'DEFG'], 'y': [1, 2, 3, 4]},
                      index=['E', 'f', 'g', 'h'])

    ddf = dd.from_pandas(df, 2)

    # Check that str not in dir/hasattr for non-object columns
    assert 'str' not in dir(ddf.y)
    assert not hasattr(ddf.y, 'str')

    # not implemented methods don't show up
    assert 'get_dummies' not in dir(ddf.x.str)
    assert not hasattr(ddf.x.str, 'get_dummies')

    assert 'upper' in dir(ddf.x.str)
    assert_eq(ddf.x.str.upper(), df.x.str.upper())
    assert set(ddf.x.str.upper().dask) == set(ddf.x.str.upper().dask)

    assert 'upper' in dir(ddf.index.str)
    assert_eq(ddf.index.str.upper(), df.index.str.upper())
    assert set(ddf.index.str.upper().dask) == set(ddf.index.str.upper().dask)

    # make sure to pass thru args & kwargs
    assert 'contains' in dir(ddf.x.str)
    assert_eq(ddf.x.str.contains('a'), df.x.str.contains('a'))
    assert set(ddf.x.str.contains('a').dask) == set(ddf.x.str.contains('a').dask)

    assert_eq(ddf.x.str.contains('d', case=False), df.x.str.contains('d', case=False))
    assert (set(ddf.x.str.contains('d', case=False).dask) ==
            set(ddf.x.str.contains('d', case=False).dask))

    for na in [True, False]:
        assert_eq(ddf.x.str.contains('a', na=na), df.x.str.contains('a', na=na))
        assert (set(ddf.x.str.contains('a', na=na).dask) ==
                set(ddf.x.str.contains('a', na=na).dask))

    for regex in [True, False]:
        assert_eq(ddf.x.str.contains('a', regex=regex), df.x.str.contains('a', regex=regex))
        assert (set(ddf.x.str.contains('a', regex=regex).dask) ==
                set(ddf.x.str.contains('a', regex=regex).dask))

    assert_eq(ddf.x.str[:2], df.x.str[:2])
    assert_eq(ddf.x.str[1], df.x.str[1])

    # str.extractall
    assert_eq(ddf.x.str.extractall('(.*)b(.*)'),
              df.x.str.extractall('(.*)b(.*)'))

    # str.cat
    sol = df.x.str.cat(df.x.str.upper(), sep=':')
    assert_eq(ddf.x.str.cat(ddf.x.str.upper(), sep=':'), sol)
    assert_eq(ddf.x.str.cat(df.x.str.upper(), sep=':'), sol)
    assert_eq(ddf.x.str.cat([ddf.x.str.upper(), df.x.str.lower()], sep=':'),
              df.x.str.cat([df.x.str.upper(), df.x.str.lower()], sep=':'))

    for o in ['foo', ['foo']]:
        with pytest.raises(TypeError):
            ddf.x.str.cat(o)

    with pytest.raises(NotImplementedError):
        ddf.x.str.cat(sep=':')


def test_empty_max():
    meta = make_meta({'x': 'i8'})
    a = dd.DataFrame({('x', 0): pd.DataFrame({'x': [1]}),
                      ('x', 1): pd.DataFrame({'x': []})}, 'x',
                     meta, [None, None, None])
    assert_eq(a.x.max(), 1)


def test_query():
    pytest.importorskip('numexpr')

    df = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [5, 6, 7, 8]})
    ddf = dd.from_pandas(df, npartitions=2)
    assert_eq(ddf.query('x**2 > y'),
              df.query('x**2 > y'))
    assert_eq(ddf.query('x**2 > @value', local_dict={'value': 4}),
              df.query('x**2 > @value', local_dict={'value': 4}))


def test_eval():
    pytest.importorskip('numexpr')

    p = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [5, 6, 7, 8]})
    d = dd.from_pandas(p, npartitions=2)

    assert_eq(p.eval('x + y'), d.eval('x + y'))
    assert_eq(p.eval('z = x + y', inplace=False),
              d.eval('z = x + y', inplace=False))
    with pytest.raises(NotImplementedError):
        d.eval('z = x + y', inplace=True)

    # catch FutureWarning from pandas about assignment in eval
    with pytest.warns(None):
        if PANDAS_VERSION < '0.21.0':
            if p.eval('z = x + y', inplace=None) is None:
                with pytest.raises(NotImplementedError):
                    d.eval('z = x + y', inplace=None)


@pytest.mark.parametrize('include, exclude', [
    ([int], None),
    (None, [int]),
    ([np.number, object], [float]),
    (['datetime'], None)
])
def test_select_dtypes(include, exclude):
    n = 10
    df = pd.DataFrame({'cint': [1] * n,
                       'cstr': ['a'] * n,
                       'clfoat': [1.] * n,
                       'cdt': pd.date_range('2016-01-01', periods=n)})
    a = dd.from_pandas(df, npartitions=2)
    result = a.select_dtypes(include=include, exclude=exclude)
    expected = df.select_dtypes(include=include, exclude=exclude)
    assert_eq(result, expected)

    # count dtypes
    tm.assert_series_equal(a.get_dtype_counts(), df.get_dtype_counts())

    tm.assert_series_equal(result.get_dtype_counts(),
                           expected.get_dtype_counts())

    if PANDAS_VERSION >= '0.23.0':
        ctx = pytest.warns(FutureWarning)
    else:
        ctx = pytest.warns(None)

    with ctx:
        tm.assert_series_equal(a.get_ftype_counts(), df.get_ftype_counts())
        tm.assert_series_equal(result.get_ftype_counts(),
                               expected.get_ftype_counts())


def test_deterministic_apply_concat_apply_names():
    df = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [5, 6, 7, 8]})
    a = dd.from_pandas(df, npartitions=2)

    assert sorted(a.x.nlargest(2).dask) == sorted(a.x.nlargest(2).dask)
    assert sorted(a.x.nlargest(2).dask) != sorted(a.x.nlargest(3).dask)
    assert (sorted(a.x.drop_duplicates().dask) ==
            sorted(a.x.drop_duplicates().dask))
    assert (sorted(a.groupby('x').y.mean().dask) ==
            sorted(a.groupby('x').y.mean().dask))

    # Test aca without passing in token string
    f = lambda a: a.nlargest(5)
    f2 = lambda a: a.nlargest(3)
    assert (sorted(aca(a.x, f, f, a.x._meta).dask) !=
            sorted(aca(a.x, f2, f2, a.x._meta).dask))
    assert (sorted(aca(a.x, f, f, a.x._meta).dask) ==
            sorted(aca(a.x, f, f, a.x._meta).dask))

    # Test aca with keywords
    def chunk(x, c_key=0, both_key=0):
        return x.sum() + c_key + both_key

    def agg(x, a_key=0, both_key=0):
        return pd.Series(x).sum() + a_key + both_key

    c_key = 2
    a_key = 3
    both_key = 4

    res = aca(a.x, chunk=chunk, aggregate=agg, chunk_kwargs={'c_key': c_key},
              aggregate_kwargs={'a_key': a_key}, both_key=both_key)
    assert (sorted(res.dask) ==
            sorted(aca(a.x, chunk=chunk, aggregate=agg,
                       chunk_kwargs={'c_key': c_key},
                       aggregate_kwargs={'a_key': a_key},
                       both_key=both_key).dask))
    assert (sorted(res.dask) !=
            sorted(aca(a.x, chunk=chunk, aggregate=agg,
                       chunk_kwargs={'c_key': c_key},
                       aggregate_kwargs={'a_key': a_key},
                       both_key=0).dask))

    assert_eq(res, df.x.sum() + 2 * (c_key + both_key) + a_key + both_key)


def test_aca_meta_infer():
    df = pd.DataFrame({'x': [1, 2, 3, 4],
                       'y': [5, 6, 7, 8]})
    ddf = dd.from_pandas(df, npartitions=2)

    def chunk(x, y, constant=1.0):
        return (x + y + constant).head()

    def agg(x):
        return x.head()

    res = aca([ddf, 2.0], chunk=chunk, aggregate=agg,
              chunk_kwargs=dict(constant=2.0))
    sol = (df + 2.0 + 2.0).head()
    assert_eq(res, sol)

    # Should infer as a scalar
    res = aca([ddf.x], chunk=lambda x: pd.Series([x.sum()]),
              aggregate=lambda x: x.sum())
    assert isinstance(res, Scalar)
    assert res.compute() == df.x.sum()


def test_aca_split_every():
    df = pd.DataFrame({'x': [1] * 60})
    ddf = dd.from_pandas(df, npartitions=15)

    def chunk(x, y, constant=0):
        return x.sum() + y + constant

    def combine(x, constant=0):
        return x.sum() + constant + 1

    def agg(x, constant=0):
        return x.sum() + constant + 2

    f = lambda n: aca([ddf, 2.0], chunk=chunk, aggregate=agg, combine=combine,
                      chunk_kwargs=dict(constant=1.0),
                      combine_kwargs=dict(constant=2.0),
                      aggregate_kwargs=dict(constant=3.0),
                      split_every=n)

    assert_max_deps(f(3), 3)
    assert_max_deps(f(4), 4, False)
    assert_max_deps(f(5), 5)
    assert set(f(15).dask.keys()) == set(f(ddf.npartitions).dask.keys())

    r3 = f(3)
    r4 = f(4)
    assert r3._name != r4._name
    # Only intersect on reading operations
    assert len(set(r3.dask.keys()) & set(r4.dask.keys())) == len(ddf.dask.keys())

    # Keywords are different for each step
    assert f(3).compute() == 60 + 15 * (2 + 1) + 7 * (2 + 1) + (3 + 2)
    # Keywords are same for each step
    res = aca([ddf, 2.0], chunk=chunk, aggregate=agg, combine=combine,
              constant=3.0, split_every=3)
    assert res.compute() == 60 + 15 * (2 + 3) + 7 * (3 + 1) + (3 + 2)
    # No combine provided, combine is agg
    res = aca([ddf, 2.0], chunk=chunk, aggregate=agg, constant=3, split_every=3)
    assert res.compute() == 60 + 15 * (2 + 3) + 8 * (3 + 2)

    # split_every must be >= 2
    with pytest.raises(ValueError):
        f(1)

    # combine_kwargs with no combine provided
    with pytest.raises(ValueError):
        aca([ddf, 2.0], chunk=chunk, aggregate=agg, split_every=3,
            chunk_kwargs=dict(constant=1.0),
            combine_kwargs=dict(constant=2.0),
            aggregate_kwargs=dict(constant=3.0))


def test_reduction_method():
    df = pd.DataFrame({'x': range(50), 'y': range(50, 100)})
    ddf = dd.from_pandas(df, npartitions=4)

    chunk = lambda x, val=0: (x >= val).sum()
    agg = lambda x: x.sum()

    # Output of chunk is a scalar
    res = ddf.x.reduction(chunk, aggregate=agg)
    assert_eq(res, df.x.count())

    # Output of chunk is a series
    res = ddf.reduction(chunk, aggregate=agg)
    assert res._name == ddf.reduction(chunk, aggregate=agg)._name
    assert_eq(res, df.count())

    # Test with keywords
    res2 = ddf.reduction(chunk, aggregate=agg, chunk_kwargs={'val': 25})
    res2._name == ddf.reduction(chunk, aggregate=agg,
                                chunk_kwargs={'val': 25})._name
    assert res2._name != res._name
    assert_eq(res2, (df >= 25).sum())

    # Output of chunk is a dataframe
    def sum_and_count(x):
        return pd.DataFrame({'sum': x.sum(), 'count': x.count()})
    res = ddf.reduction(sum_and_count,
                        aggregate=lambda x: x.groupby(level=0).sum())

    assert_eq(res, pd.DataFrame({'sum': df.sum(), 'count': df.count()}))


def test_reduction_method_split_every():
    df = pd.Series([1] * 60)
    ddf = dd.from_pandas(df, npartitions=15)

    def chunk(x, constant=0):
        return x.sum() + constant

    def combine(x, constant=0):
        return x.sum() + constant + 1

    def agg(x, constant=0):
        return x.sum() + constant + 2

    f = lambda n: ddf.reduction(chunk, aggregate=agg, combine=combine,
                                chunk_kwargs=dict(constant=1.0),
                                combine_kwargs=dict(constant=2.0),
                                aggregate_kwargs=dict(constant=3.0),
                                split_every=n)

    assert_max_deps(f(3), 3)
    assert_max_deps(f(4), 4, False)
    assert_max_deps(f(5), 5)
    assert set(f(15).dask.keys()) == set(f(ddf.npartitions).dask.keys())

    r3 = f(3)
    r4 = f(4)
    assert r3._name != r4._name
    # Only intersect on reading operations
    assert len(set(r3.dask.keys()) & set(r4.dask.keys())) == len(ddf.dask.keys())

    # Keywords are different for each step
    assert f(3).compute() == 60 + 15 + 7 * (2 + 1) + (3 + 2)
    # Keywords are same for each step
    res = ddf.reduction(chunk, aggregate=agg, combine=combine, constant=3.0,
                        split_every=3)
    assert res.compute() == 60 + 15 * 3 + 7 * (3 + 1) + (3 + 2)
    # No combine provided, combine is agg
    res = ddf.reduction(chunk, aggregate=agg, constant=3.0, split_every=3)
    assert res.compute() == 60 + 15 * 3 + 8 * (3 + 2)

    # split_every must be >= 2
    with pytest.raises(ValueError):
        f(1)

    # combine_kwargs with no combine provided
    with pytest.raises(ValueError):
        ddf.reduction(chunk, aggregate=agg, split_every=3,
                      chunk_kwargs=dict(constant=1.0),
                      combine_kwargs=dict(constant=2.0),
                      aggregate_kwargs=dict(constant=3.0))


def test_pipe():
    df = pd.DataFrame({'x': range(50), 'y': range(50, 100)})
    ddf = dd.from_pandas(df, npartitions=4)

    def f(x, y, z=0):
        return x + y + z

    assert_eq(ddf.pipe(f, 1, z=2), f(ddf, 1, z=2))
    assert_eq(ddf.x.pipe(f, 1, z=2), f(ddf.x, 1, z=2))


def test_gh_517():
    arr = np.random.randn(100, 2)
    df = pd.DataFrame(arr, columns=['a', 'b'])
    ddf = dd.from_pandas(df, 2)
    assert ddf.index.nunique().compute() == 100

    ddf2 = dd.from_pandas(pd.concat([df, df]), 5)
    assert ddf2.index.nunique().compute() == 100


def test_drop_axis_1():
    df = pd.DataFrame({'x': [1, 2, 3, 4],
                       'y': [5, 6, 7, 8],
                       'z': [9, 10, 11, 12]})
    ddf = dd.from_pandas(df, npartitions=2)

    assert_eq(ddf.drop('y', axis=1), df.drop('y', axis=1))
    assert_eq(ddf.drop(['y', 'z'], axis=1), df.drop(['y', 'z'], axis=1))
    with pytest.raises(ValueError):
        ddf.drop(['a', 'x'], axis=1)
    assert_eq(ddf.drop(['a', 'x'], axis=1, errors='ignore'),
              df.drop(['a', 'x'], axis=1, errors='ignore'))


def test_gh580():
    df = pd.DataFrame({'x': np.arange(10, dtype=float)})
    ddf = dd.from_pandas(df, 2)
    assert_eq(np.cos(df['x']), np.cos(ddf['x']))
    assert_eq(np.cos(df['x']), np.cos(ddf['x']))


def test_rename_dict():
    renamer = {'a': 'A', 'b': 'B'}
    assert_eq(d.rename(columns=renamer),
              full.rename(columns=renamer))


def test_rename_function():
    renamer = lambda x: x.upper()
    assert_eq(d.rename(columns=renamer),
              full.rename(columns=renamer))


def test_rename_index():
    renamer = {0: 1}
    pytest.raises(ValueError, lambda: d.rename(index=renamer))


def test_to_timestamp():
    index = pd.PeriodIndex(freq='A', start='1/1/2001', end='12/1/2004')
    df = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [10, 20, 30, 40]}, index=index)
    ddf = dd.from_pandas(df, npartitions=3)
    assert_eq(ddf.to_timestamp(), df.to_timestamp())
    assert_eq(ddf.to_timestamp(freq='M', how='s').compute(),
              df.to_timestamp(freq='M', how='s'))
    assert_eq(ddf.x.to_timestamp(), df.x.to_timestamp())
    assert_eq(ddf.x.to_timestamp(freq='M', how='s').compute(),
              df.x.to_timestamp(freq='M', how='s'))


def test_to_frame():
    s = pd.Series([1, 2, 3], name='foo')
    a = dd.from_pandas(s, npartitions=2)

    assert_eq(s.to_frame(), a.to_frame())
    assert_eq(s.to_frame('bar'), a.to_frame('bar'))


def test_apply():
    df = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)

    func = lambda row: row['x'] + row['y']
    assert_eq(ddf.x.apply(lambda x: x + 1, meta=("x", int)),
              df.x.apply(lambda x: x + 1))

    # specify meta
    assert_eq(ddf.apply(lambda xy: xy[0] + xy[1], axis=1, meta=(None, int)),
              df.apply(lambda xy: xy[0] + xy[1], axis=1))
    assert_eq(ddf.apply(lambda xy: xy[0] + xy[1], axis='columns', meta=(None, int)),
              df.apply(lambda xy: xy[0] + xy[1], axis='columns'))

    # inference
    with pytest.warns(None):
        assert_eq(ddf.apply(lambda xy: xy[0] + xy[1], axis=1),
                  df.apply(lambda xy: xy[0] + xy[1], axis=1))
    with pytest.warns(None):
        assert_eq(ddf.apply(lambda xy: xy, axis=1),
                  df.apply(lambda xy: xy, axis=1))

    # specify meta
    func = lambda x: pd.Series([x, x])
    assert_eq(ddf.x.apply(func, meta=[(0, int), (1, int)]), df.x.apply(func))
    # inference
    with pytest.warns(None):
        assert_eq(ddf.x.apply(func), df.x.apply(func))

    # axis=0
    with pytest.raises(NotImplementedError):
        ddf.apply(lambda xy: xy, axis=0)

    with pytest.raises(NotImplementedError):
        ddf.apply(lambda xy: xy, axis='index')


@pytest.mark.skipif(sys.version_info <= (3, 0),
                    reason="Global filter is applied by another library, and "
                           "not reset properly.")
def test_apply_warns():
    df = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)

    func = lambda row: row['x'] + row['y']

    with pytest.warns(UserWarning) as w:
        ddf.apply(func, axis=1)
    assert len(w) == 1

    with pytest.warns(None) as w:
        ddf.apply(func, axis=1, meta=(None, int))
    assert len(w) == 0


def test_applymap():
    df = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)
    assert_eq(ddf.applymap(lambda x: x + 1), df.applymap(lambda x: x + 1))

    assert_eq(ddf.applymap(lambda x: (x, x)), df.applymap(lambda x: (x, x)))


def test_abs():
    df = pd.DataFrame({'A': [1, -2, 3, -4, 5],
                       'B': [-6., -7, -8, -9, 10],
                       'C': ['a', 'b', 'c', 'd', 'e']})
    ddf = dd.from_pandas(df, npartitions=2)
    assert_eq(ddf.A.abs(), df.A.abs())
    assert_eq(ddf[['A', 'B']].abs(), df[['A', 'B']].abs())
    pytest.raises(ValueError, lambda: ddf.C.abs())
    pytest.raises(TypeError, lambda: ddf.abs())


def test_round():
    df = pd.DataFrame({'col1': [1.123, 2.123, 3.123],
                       'col2': [1.234, 2.234, 3.234]})
    ddf = dd.from_pandas(df, npartitions=2)
    assert_eq(ddf.round(), df.round())
    assert_eq(ddf.round(2), df.round(2))


def test_cov():
    # DataFrame
    df = pd.util.testing.makeMissingDataframe(0.3, 42)
    ddf = dd.from_pandas(df, npartitions=6)

    res = ddf.cov()
    res2 = ddf.cov(split_every=2)
    res3 = ddf.cov(10)
    res4 = ddf.cov(10, split_every=2)
    sol = df.cov()
    sol2 = df.cov(10)
    assert_eq(res, sol)
    assert_eq(res2, sol)
    assert_eq(res3, sol2)
    assert_eq(res4, sol2)
    assert res._name == ddf.cov()._name
    assert res._name != res2._name
    assert res3._name != res4._name
    assert res._name != res3._name

    # Series
    a = df.A
    b = df.B
    da = dd.from_pandas(a, npartitions=6)
    db = dd.from_pandas(b, npartitions=7)

    res = da.cov(db)
    res2 = da.cov(db, split_every=2)
    res3 = da.cov(db, 10)
    res4 = da.cov(db, 10, split_every=2)
    sol = a.cov(b)
    sol2 = a.cov(b, 10)
    assert_eq(res, sol)
    assert_eq(res2, sol)
    assert_eq(res3, sol2)
    assert_eq(res4, sol2)
    assert res._name == da.cov(db)._name
    assert res._name != res2._name
    assert res3._name != res4._name
    assert res._name != res3._name


def test_corr():
    # DataFrame
    df = pd.util.testing.makeMissingDataframe(0.3, 42)
    ddf = dd.from_pandas(df, npartitions=6)

    res = ddf.corr()
    res2 = ddf.corr(split_every=2)
    res3 = ddf.corr(min_periods=10)
    res4 = ddf.corr(min_periods=10, split_every=2)
    sol = df.corr()
    sol2 = df.corr(min_periods=10)
    assert_eq(res, sol)
    assert_eq(res2, sol)
    assert_eq(res3, sol2)
    assert_eq(res4, sol2)
    assert res._name == ddf.corr()._name
    assert res._name != res2._name
    assert res3._name != res4._name
    assert res._name != res3._name

    pytest.raises(NotImplementedError, lambda: ddf.corr(method='spearman'))

    # Series
    a = df.A
    b = df.B
    da = dd.from_pandas(a, npartitions=6)
    db = dd.from_pandas(b, npartitions=7)

    res = da.corr(db)
    res2 = da.corr(db, split_every=2)
    res3 = da.corr(db, min_periods=10)
    res4 = da.corr(db, min_periods=10, split_every=2)
    sol = da.corr(db)
    sol2 = da.corr(db, min_periods=10)
    assert_eq(res, sol)
    assert_eq(res2, sol)
    assert_eq(res3, sol2)
    assert_eq(res4, sol2)
    assert res._name == da.corr(db)._name
    assert res._name != res2._name
    assert res3._name != res4._name
    assert res._name != res3._name

    pytest.raises(NotImplementedError, lambda: da.corr(db, method='spearman'))
    pytest.raises(TypeError, lambda: da.corr(ddf))


def test_cov_corr_meta():
    df = pd.DataFrame({'a': np.array([1, 2, 3]),
                       'b': np.array([1.0, 2.0, 3.0], dtype='f4'),
                       'c': np.array([1.0, 2.0, 3.0])},
                      index=pd.Index([1, 2, 3], name='myindex'))
    ddf = dd.from_pandas(df, npartitions=2)
    assert_eq(ddf.corr(), df.corr())
    assert_eq(ddf.cov(), df.cov())
    assert ddf.a.cov(ddf.b)._meta.dtype == 'f8'
    assert ddf.a.corr(ddf.b)._meta.dtype == 'f8'


@pytest.mark.slow
def test_cov_corr_stable():
    df = pd.DataFrame(np.random.uniform(-1, 1, (20000000, 2)), columns=['a', 'b'])
    ddf = dd.from_pandas(df, npartitions=50)
    assert_eq(ddf.cov(split_every=8), df.cov())
    assert_eq(ddf.corr(split_every=8), df.corr())


def test_cov_corr_mixed():
    size = 1000
    d = {'dates' : pd.date_range('2015-01-01', periods=size, freq='1T'),
         'unique_id' : np.arange(0, size),
         'ints' : np.random.randint(0, size, size=size),
         'floats' : np.random.randn(size),
         'bools' : np.random.choice([0, 1], size=size),
         'int_nans' : np.random.choice([0, 1, np.nan], size=size),
         'float_nans' : np.random.choice([0.0, 1.0, np.nan], size=size),
         'constant' : 1,
         'int_categorical' : np.random.choice([10, 20, 30, 40, 50], size=size) ,
         'categorical_binary' : np.random.choice(['a', 'b'], size=size),
         'categorical_nans' : np.random.choice(['a', 'b', 'c'], size=size)}
    df = pd.DataFrame(d)
    df['hardbools'] = df['bools'] == 1
    df['categorical_nans'] = df['categorical_nans'].replace('c', np.nan)
    df['categorical_binary'] = df['categorical_binary'].astype('category')
    df['unique_id'] = df['unique_id'].astype(str)

    ddf = dd.from_pandas(df, npartitions=20)
    assert_eq(ddf.corr(split_every=4), df.corr(), check_divisions=False)
    assert_eq(ddf.cov(split_every=4), df.cov(), check_divisions=False)


def test_autocorr():
    x = pd.Series(np.random.random(100))
    dx = dd.from_pandas(x, npartitions=10)
    assert_eq(dx.autocorr(2), x.autocorr(2))
    assert_eq(dx.autocorr(0), x.autocorr(0))
    assert_eq(dx.autocorr(-2), x.autocorr(-2))
    assert_eq(dx.autocorr(2, split_every=3), x.autocorr(2))
    pytest.raises(TypeError, lambda: dx.autocorr(1.5))


def test_apply_infer_columns():
    df = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)

    def return_df(x):
        # will create new DataFrame which columns is ['sum', 'mean']
        return pd.Series([x.sum(), x.mean()], index=['sum', 'mean'])

    # DataFrame to completely different DataFrame
    with pytest.warns(None):
        result = ddf.apply(return_df, axis=1)
    assert isinstance(result, dd.DataFrame)
    tm.assert_index_equal(result.columns, pd.Index(['sum', 'mean']))
    assert_eq(result, df.apply(return_df, axis=1))

    # DataFrame to Series
    with pytest.warns(None):
        result = ddf.apply(lambda x: 1, axis=1)
    assert isinstance(result, dd.Series)
    assert result.name is None
    assert_eq(result, df.apply(lambda x: 1, axis=1))

    def return_df2(x):
        return pd.Series([x * 2, x * 3], index=['x2', 'x3'])

    # Series to completely different DataFrame
    with pytest.warns(None):
        result = ddf.x.apply(return_df2)
    assert isinstance(result, dd.DataFrame)
    tm.assert_index_equal(result.columns, pd.Index(['x2', 'x3']))
    assert_eq(result, df.x.apply(return_df2))

    # Series to Series
    with pytest.warns(None):
        result = ddf.x.apply(lambda x: 1)
    assert isinstance(result, dd.Series)
    assert result.name == 'x'
    assert_eq(result, df.x.apply(lambda x: 1))


def test_index_time_properties():
    i = tm.makeTimeSeries()
    a = dd.from_pandas(i, npartitions=3)

    assert 'day' in dir(a.index)
    # returns a numpy array in pandas, but a Index in dask
    assert_eq(a.index.day, pd.Index(i.index.day))
    assert_eq(a.index.month, pd.Index(i.index.month))


def test_nlargest_nsmallest():
    from string import ascii_lowercase
    df = pd.DataFrame({'a': np.random.permutation(20),
                       'b': list(ascii_lowercase[:20]),
                       'c': np.random.permutation(20).astype('float64')})
    ddf = dd.from_pandas(df, npartitions=3)

    for m in ['nlargest', 'nsmallest']:
        f = lambda df, *args, **kwargs: getattr(df, m)(*args, **kwargs)

        res = f(ddf, 5, 'a')
        res2 = f(ddf, 5, 'a', split_every=2)
        sol = f(df, 5, 'a')
        assert_eq(res, sol)
        assert_eq(res2, sol)
        assert res._name != res2._name

        res = f(ddf, 5, ['a', 'c'])
        res2 = f(ddf, 5, ['a', 'c'], split_every=2)
        sol = f(df, 5, ['a', 'c'])
        assert_eq(res, sol)
        assert_eq(res2, sol)
        assert res._name != res2._name

        res = f(ddf.a, 5)
        res2 = f(ddf.a, 5, split_every=2)
        sol = f(df.a, 5)
        assert_eq(res, sol)
        assert_eq(res2, sol)
        assert res._name != res2._name


def test_reset_index():
    df = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)

    sol = df.reset_index()
    res = ddf.reset_index()
    assert all(d is None for d in res.divisions)
    assert_eq(res, sol, check_index=False)

    sol = df.reset_index(drop=True)
    res = ddf.reset_index(drop=True)
    assert all(d is None for d in res.divisions)
    assert_eq(res, sol, check_index=False)

    sol = df.x.reset_index()
    res = ddf.x.reset_index()
    assert all(d is None for d in res.divisions)
    assert_eq(res, sol, check_index=False)

    sol = df.x.reset_index(drop=True)
    res = ddf.x.reset_index(drop=True)
    assert all(d is None for d in res.divisions)
    assert_eq(res, sol, check_index=False)


def test_dataframe_compute_forward_kwargs():
    x = dd.from_pandas(pd.DataFrame({'a': range(10)}), npartitions=2).a.sum()
    x.compute(bogus_keyword=10)


def test_series_iteritems():
    df = pd.DataFrame({'x': [1, 2, 3, 4]})
    ddf = dd.from_pandas(df, npartitions=2)
    for (a, b) in zip(df['x'].iteritems(), ddf['x'].iteritems()):
        assert a == b


def test_dataframe_iterrows():
    df = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)

    for (a, b) in zip(df.iterrows(), ddf.iterrows()):
        tm.assert_series_equal(a[1], b[1])


def test_dataframe_itertuples():
    df = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)

    for (a, b) in zip(df.itertuples(), ddf.itertuples()):
        assert a == b


def test_astype():
    df = pd.DataFrame({'x': [1, 2, 3, None], 'y': [10, 20, 30, 40]},
                      index=[10, 20, 30, 40])
    a = dd.from_pandas(df, 2)

    assert_eq(a.astype(float), df.astype(float))
    assert_eq(a.x.astype(float), df.x.astype(float))


def test_astype_categoricals():
    df = pd.DataFrame({'x': ['a', 'b', 'c', 'b', 'c'],
                       'y': ['x', 'y', 'z', 'x', 'y'],
                       'z': [1, 2, 3, 4, 5]})
    df = df.astype({'y': 'category'})
    ddf = dd.from_pandas(df, 2)
    assert ddf.y.cat.known

    ddf2 = ddf.astype({'x': 'category'})
    assert not ddf2.x.cat.known
    assert ddf2.y.cat.known
    assert ddf2.x.dtype == 'category'
    assert ddf2.compute().x.dtype == 'category'

    dx = ddf.x.astype('category')
    assert not dx.cat.known
    assert dx.dtype == 'category'
    assert dx.compute().dtype == 'category'


@pytest.mark.skipif(PANDAS_VERSION < '0.21.0',
                    reason="No CategoricalDtype with categories")
def test_astype_categoricals_known():
    df = pd.DataFrame({'x': ['a', 'b', 'c', 'b', 'c'],
                       'y': ['x', 'y', 'z', 'y', 'z'],
                       'z': ['b', 'b', 'b', 'c', 'b'],
                       'other': [1, 2, 3, 4, 5]})
    ddf = dd.from_pandas(df, 2)

    abc = pd.api.types.CategoricalDtype(['a', 'b', 'c'])
    category = pd.api.types.CategoricalDtype()

    # DataFrame
    ddf2 = ddf.astype({'x': abc,
                       'y': category,
                       'z': 'category',
                       'other': 'f8'})

    for col, known in [('x', True), ('y', False), ('z', False)]:
        x = getattr(ddf2, col)
        assert pd.api.types.is_categorical_dtype(x.dtype)
        assert x.cat.known == known

    # Series
    for dtype, known in [('category', False), (category, False), (abc, True)]:
        dx2 = ddf.x.astype(dtype)
        assert pd.api.types.is_categorical_dtype(dx2.dtype)
        assert dx2.cat.known == known


def test_groupby_callable():
    a = pd.DataFrame({'x': [1, 2, 3, None], 'y': [10, 20, 30, 40]},
                     index=[1, 2, 3, 4])
    b = dd.from_pandas(a, 2)

    def iseven(x):
        return x % 2 == 0

    assert_eq(a.groupby(iseven).y.sum(),
              b.groupby(iseven).y.sum())
    assert_eq(a.y.groupby(iseven).sum(),
              b.y.groupby(iseven).sum())


def test_methods_tokenize_differently():
    df = pd.DataFrame({'x': [1, 2, 3, 4]})
    df = dd.from_pandas(df, npartitions=1)
    assert (df.x.map_partitions(lambda x: pd.Series(x.min()))._name !=
            df.x.map_partitions(lambda x: pd.Series(x.max()))._name)


def _assert_info(df, ddf, memory_usage=True):
    from io import StringIO
    assert isinstance(df, pd.DataFrame)
    assert isinstance(ddf, dd.DataFrame)

    buf_pd, buf_da = StringIO(), StringIO()

    df.info(buf=buf_pd, memory_usage=memory_usage)
    ddf.info(buf=buf_da, verbose=True, memory_usage=memory_usage)

    stdout_pd = buf_pd.getvalue()
    stdout_da = buf_da.getvalue()
    stdout_da = stdout_da.replace(str(type(ddf)), str(type(df)))
    assert stdout_pd == stdout_da


def test_info():
    from io import StringIO
    from dask.compatibility import unicode
    pandas_format._put_lines = put_lines

    test_frames = [
        pd.DataFrame({'x': [1, 2, 3, 4], 'y': [1, 0, 1, 0]}, index=pd.Int64Index(range(4))),  # No RangeIndex in dask
        pd.DataFrame()
    ]

    for df in test_frames:
        ddf = dd.from_pandas(df, npartitions=4)
        _assert_info(df, ddf)

    buf = StringIO()
    ddf = dd.from_pandas(pd.DataFrame({'x': [1, 2, 3, 4], 'y': [1, 0, 1, 0]}, index=range(4)), npartitions=4)

    # Verbose=False
    ddf.info(buf=buf, verbose=False)
    assert buf.getvalue() == unicode("<class 'dask.dataframe.core.DataFrame'>\n"
                                     "Columns: 2 entries, x to y\n"
                                     "dtypes: int64(2)")

    # buf=None
    assert ddf.info(buf=None) is None


def test_groupby_multilevel_info():
    # GH 1844
    from io import StringIO
    from dask.compatibility import unicode

    pandas_format._put_lines = put_lines

    df = pd.DataFrame({'A': [1, 1, 2, 2],
                       'B': [1, 2, 3, 4],
                       'C': [1, 2, 3, 4]})
    ddf = dd.from_pandas(df, npartitions=2)

    g = ddf.groupby(['A', 'B']).sum()
    # slight difference between memory repr (single additional space)
    _assert_info(g.compute(), g, memory_usage=False)

    buf = StringIO()
    g.info(buf, verbose=False)
    assert buf.getvalue() == unicode("""<class 'dask.dataframe.core.DataFrame'>
Columns: 1 entries, C to C
dtypes: int64(1)""")

    # multilevel
    g = ddf.groupby(['A', 'B']).agg(['count', 'sum'])
    _assert_info(g.compute(), g, memory_usage=False)

    buf = StringIO()
    g.info(buf, verbose=False)
    expected = unicode(textwrap.dedent("""\
    <class 'dask.dataframe.core.DataFrame'>
    Columns: 2 entries, ('C', 'count') to ('C', 'sum')
    dtypes: int64(2)"""))
    assert buf.getvalue() == expected


def test_categorize_info():
    # assert that we can call info after categorize
    # workaround for: https://github.com/pydata/pandas/issues/14368
    from io import StringIO
    from dask.compatibility import unicode
    pandas_format._put_lines = put_lines

    df = pd.DataFrame({'x': [1, 2, 3, 4],
                       'y': pd.Series(list('aabc')),
                       'z': pd.Series(list('aabc'))},
                      index=pd.Int64Index(range(4)))  # No RangeIndex in dask
    ddf = dd.from_pandas(df, npartitions=4).categorize(['y'])

    # Verbose=False
    buf = StringIO()
    ddf.info(buf=buf, verbose=True)
    expected = unicode("<class 'dask.dataframe.core.DataFrame'>\n"
                       "Int64Index: 4 entries, 0 to 3\n"
                       "Data columns (total 3 columns):\n"
                       "x    4 non-null int64\n"
                       "y    4 non-null category\n"
                       "z    4 non-null object\n"
                       "dtypes: category(1), object(1), int64(1)")
    assert buf.getvalue() == expected


def test_gh_1301():
    df = pd.DataFrame([['1', '2'], ['3', '4']])
    ddf = dd.from_pandas(df, npartitions=2)
    ddf2 = ddf.assign(y=ddf[1].astype(int))
    assert_eq(ddf2, df.assign(y=df[1].astype(int)))

    assert ddf2.dtypes['y'] == np.dtype(int)


def test_timeseries_sorted():
    df = tm.makeTimeDataFrame()
    ddf = dd.from_pandas(df.reset_index(), npartitions=2)
    df.index.name = 'index'
    assert_eq(ddf.set_index('index', sorted=True, drop=True), df)


def test_column_assignment():
    df = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [1, 0, 1, 0]})
    ddf = dd.from_pandas(df, npartitions=2)
    orig = ddf.copy()
    ddf['z'] = ddf.x + ddf.y
    df['z'] = df.x + df.y

    assert_eq(df, ddf)
    assert 'z' not in orig.columns


def test_columns_assignment():
    df = pd.DataFrame({'x': [1, 2, 3, 4]})
    ddf = dd.from_pandas(df, npartitions=2)

    df2 = df.assign(y=df.x + 1, z=df.x - 1)
    df[['a', 'b']] = df2[['y', 'z']]

    ddf2 = ddf.assign(y=ddf.x + 1, z=ddf.x - 1)
    ddf[['a', 'b']] = ddf2[['y', 'z']]

    assert_eq(df, ddf)


def test_attribute_assignment():
    df = pd.DataFrame({'x': [1, 2, 3, 4, 5],
                       'y': [1., 2., 3., 4., 5.]})
    ddf = dd.from_pandas(df, npartitions=2)

    ddf.y = ddf.x + ddf.y
    assert_eq(ddf, df.assign(y=df.x + df.y))


def test_setitem_triggering_realign():
    a = dd.from_pandas(pd.DataFrame({"A": range(12)}), npartitions=3)
    b = dd.from_pandas(pd.Series(range(12), name='B'), npartitions=4)
    a['C'] = b
    assert len(a) == 12


def test_inplace_operators():
    df = pd.DataFrame({'x': [1, 2, 3, 4, 5],
                       'y': [1., 2., 3., 4., 5.]})
    ddf = dd.from_pandas(df, npartitions=2)

    ddf.y **= 0.5

    assert_eq(ddf.y, df.y ** 0.5)
    assert_eq(ddf, df.assign(y=df.y ** 0.5))


@pytest.mark.parametrize("skipna", [True, False])
@pytest.mark.parametrize("idx", [
    np.arange(100),
    sorted(np.random.random(size=100)),
    pd.date_range('20150101', periods=100)
])
def test_idxmaxmin(idx, skipna):
    df = pd.DataFrame(np.random.randn(100, 5), columns=list('abcde'), index=idx)
    df.b.iloc[31] = np.nan
    df.d.iloc[78] = np.nan
    ddf = dd.from_pandas(df, npartitions=3)

    assert_eq(df.idxmax(axis=1, skipna=skipna),
              ddf.idxmax(axis=1, skipna=skipna))
    assert_eq(df.idxmin(axis=1, skipna=skipna),
              ddf.idxmin(axis=1, skipna=skipna))

    assert_eq(df.idxmax(skipna=skipna), ddf.idxmax(skipna=skipna))
    assert_eq(df.idxmax(skipna=skipna),
              ddf.idxmax(skipna=skipna, split_every=2))
    assert (ddf.idxmax(skipna=skipna)._name !=
            ddf.idxmax(skipna=skipna, split_every=2)._name)

    assert_eq(df.idxmin(skipna=skipna), ddf.idxmin(skipna=skipna))
    assert_eq(df.idxmin(skipna=skipna),
              ddf.idxmin(skipna=skipna, split_every=2))
    assert (ddf.idxmin(skipna=skipna)._name !=
            ddf.idxmin(skipna=skipna, split_every=2)._name)

    assert_eq(df.a.idxmax(skipna=skipna), ddf.a.idxmax(skipna=skipna))
    assert_eq(df.a.idxmax(skipna=skipna),
              ddf.a.idxmax(skipna=skipna, split_every=2))
    assert (ddf.a.idxmax(skipna=skipna)._name !=
            ddf.a.idxmax(skipna=skipna, split_every=2)._name)

    assert_eq(df.a.idxmin(skipna=skipna), ddf.a.idxmin(skipna=skipna))
    assert_eq(df.a.idxmin(skipna=skipna),
              ddf.a.idxmin(skipna=skipna, split_every=2))
    assert (ddf.a.idxmin(skipna=skipna)._name !=
            ddf.a.idxmin(skipna=skipna, split_every=2)._name)


def test_idxmaxmin_empty_partitions():
    df = pd.DataFrame({'a': [1, 2, 3],
                       'b': [1.5, 2, 3],
                       'c': [np.NaN] * 3,
                       'd': [1, 2, np.NaN]})
    empty = df.iloc[:0]

    ddf = dd.concat([dd.from_pandas(df, npartitions=1)] +
                    [dd.from_pandas(empty, npartitions=1)] * 10)

    for skipna in [True, False]:
        assert_eq(ddf.idxmin(skipna=skipna, split_every=3),
                  df.idxmin(skipna=skipna))

    assert_eq(ddf[['a', 'b', 'd']].idxmin(skipna=skipna, split_every=3),
              df[['a', 'b', 'd']].idxmin(skipna=skipna))

    assert_eq(ddf.b.idxmax(split_every=3), df.b.idxmax())

    # Completely empty raises
    ddf = dd.concat([dd.from_pandas(empty, npartitions=1)] * 10)
    with pytest.raises(ValueError):
        ddf.idxmax().compute()
    with pytest.raises(ValueError):
        ddf.b.idxmax().compute()


def test_getitem_meta():
    data = {'col1': ['a', 'a', 'b'],
            'col2': [0, 1, 0]}

    df = pd.DataFrame(data=data, columns=['col1', 'col2'])
    ddf = dd.from_pandas(df, npartitions=1)

    assert_eq(df.col2[df.col1 == 'a'], ddf.col2[ddf.col1 == 'a'])


def test_getitem_multilevel():
    pdf = pd.DataFrame({('A', '0') : [1,2,2], ('B', '1') : [1,2,3]})
    ddf = dd.from_pandas(pdf, npartitions=3)

    assert_eq(pdf['A', '0'], ddf['A', '0'])
    assert_eq(pdf[[('A', '0'), ('B', '1')]], ddf[[('A', '0'), ('B', '1')]])


def test_getitem_string_subclass():
    df = pd.DataFrame({'column_1': list(range(10))})
    ddf = dd.from_pandas(df, npartitions=3)

    class string_subclass(str):
        pass
    column_1 = string_subclass('column_1')

    assert_eq(df[column_1], ddf[column_1])


def test_diff():
    df = pd.DataFrame(np.random.randn(100, 5), columns=list('abcde'))
    ddf = dd.from_pandas(df, 5)

    assert_eq(ddf.diff(), df.diff())
    assert_eq(ddf.diff(0), df.diff(0))
    assert_eq(ddf.diff(2), df.diff(2))
    assert_eq(ddf.diff(-2), df.diff(-2))

    assert_eq(ddf.diff(2, axis=1), df.diff(2, axis=1))

    assert_eq(ddf.a.diff(), df.a.diff())
    assert_eq(ddf.a.diff(0), df.a.diff(0))
    assert_eq(ddf.a.diff(2), df.a.diff(2))
    assert_eq(ddf.a.diff(-2), df.a.diff(-2))

    assert ddf.diff(2)._name == ddf.diff(2)._name
    assert ddf.diff(2)._name != ddf.diff(3)._name
    pytest.raises(TypeError, lambda: ddf.diff(1.5))


def test_shift():
    df = tm.makeTimeDataFrame()
    ddf = dd.from_pandas(df, npartitions=4)

    # DataFrame
    assert_eq(ddf.shift(), df.shift())
    assert_eq(ddf.shift(0), df.shift(0))
    assert_eq(ddf.shift(2), df.shift(2))
    assert_eq(ddf.shift(-2), df.shift(-2))

    assert_eq(ddf.shift(2, axis=1), df.shift(2, axis=1))

    # Series
    assert_eq(ddf.A.shift(), df.A.shift())
    assert_eq(ddf.A.shift(0), df.A.shift(0))
    assert_eq(ddf.A.shift(2), df.A.shift(2))
    assert_eq(ddf.A.shift(-2), df.A.shift(-2))

    with pytest.raises(TypeError):
        ddf.shift(1.5)


def test_shift_with_freq():
    df = tm.makeTimeDataFrame(30)
    # DatetimeIndex
    for data_freq, divs1 in [('B', False), ('D', True), ('H', True)]:
        df = df.set_index(tm.makeDateIndex(30, freq=data_freq))
        ddf = dd.from_pandas(df, npartitions=4)
        for freq, divs2 in [('S', True), ('W', False),
                            (pd.Timedelta(10, unit='h'), True)]:
            for d, p in [(ddf, df), (ddf.A, df.A), (ddf.index, df.index)]:
                res = d.shift(2, freq=freq)
                assert_eq(res, p.shift(2, freq=freq))
                assert res.known_divisions == divs2
        # Index shifts also work with freq=None
        res = ddf.index.shift(2)
        assert_eq(res, df.index.shift(2))
        assert res.known_divisions == divs1

    # PeriodIndex
    for data_freq, divs in [('B', False), ('D', True), ('H', True)]:
        df = df.set_index(pd.period_range('2000-01-01', periods=30,
                                          freq=data_freq))
        ddf = dd.from_pandas(df, npartitions=4)
        for d, p in [(ddf, df), (ddf.A, df.A)]:
            res = d.shift(2, freq=data_freq)
            assert_eq(res, p.shift(2, freq=data_freq))
            assert res.known_divisions == divs
        # PeriodIndex.shift doesn't have `freq` parameter
        res = ddf.index.shift(2)
        assert_eq(res, df.index.shift(2))
        assert res.known_divisions == divs

    with pytest.raises(ValueError):
        ddf.index.shift(2, freq='D')  # freq keyword not supported

    # TimedeltaIndex
    for data_freq in ['T', 'D', 'H']:
        df = df.set_index(tm.makeTimedeltaIndex(30, freq=data_freq))
        ddf = dd.from_pandas(df, npartitions=4)
        for freq in ['S', pd.Timedelta(10, unit='h')]:
            for d, p in [(ddf, df), (ddf.A, df.A), (ddf.index, df.index)]:
                res = d.shift(2, freq=freq)
                assert_eq(res, p.shift(2, freq=freq))
                assert res.known_divisions
        # Index shifts also work with freq=None
        res = ddf.index.shift(2)
        assert_eq(res, df.index.shift(2))
        assert res.known_divisions

    # Other index types error
    df = tm.makeDataFrame()
    ddf = dd.from_pandas(df, npartitions=4)
    pytest.raises(NotImplementedError, lambda: ddf.shift(2, freq='S'))
    pytest.raises(NotImplementedError, lambda: ddf.A.shift(2, freq='S'))
    pytest.raises(NotImplementedError, lambda: ddf.index.shift(2))


@pytest.mark.parametrize('method', ['first', 'last'])
def test_first_and_last(method):
    f = lambda x, offset: getattr(x, method)(offset)
    freqs = ['12h', 'D']
    offsets = ['0d', '100h', '20d', '20B', '3W', '3M', '400d', '13M']
    for freq in freqs:
        index = pd.date_range('1/1/2000', '1/1/2001', freq=freq)[::4]
        df = pd.DataFrame(np.random.random((len(index), 4)), index=index,
                          columns=['A', 'B', 'C', 'D'])
        ddf = dd.from_pandas(df, npartitions=10)
        for offset in offsets:
            assert_eq(f(ddf, offset), f(df, offset))
            assert_eq(f(ddf.A, offset), f(df.A, offset))


@pytest.mark.parametrize('npartitions', [1, 4, 20])
@pytest.mark.parametrize('split_every', [2, 5])
@pytest.mark.parametrize('split_out', [None, 1, 5, 20])
def test_hash_split_unique(npartitions, split_every, split_out):
    from string import ascii_lowercase
    s = pd.Series(np.random.choice(list(ascii_lowercase), 1000, replace=True))
    ds = dd.from_pandas(s, npartitions=npartitions)

    dropped = ds.unique(split_every=split_every, split_out=split_out)

    dsk = dropped.__dask_optimize__(dropped.dask, dropped.__dask_keys__())
    from dask.core import get_deps
    dependencies, dependents = get_deps(dsk)

    assert len([k for k, v in dependencies.items() if not v]) == npartitions
    assert dropped.npartitions == (split_out or 1)
    assert sorted(dropped.compute(get=dask.get)) == sorted(s.unique())


@pytest.mark.parametrize('split_every', [None, 2])
def test_split_out_drop_duplicates(split_every):
    x = np.concatenate([np.arange(10)] * 100)[:, None]
    y = x.copy()
    z = np.concatenate([np.arange(20)] * 50)[:, None]
    rs = np.random.RandomState(1)
    rs.shuffle(x)
    rs.shuffle(y)
    rs.shuffle(z)
    df = pd.DataFrame(np.concatenate([x, y, z], axis=1), columns=['x', 'y', 'z'])
    ddf = dd.from_pandas(df, npartitions=20)

    for subset, keep in product([None, ['x', 'z']], ['first', 'last']):
        sol = df.drop_duplicates(subset=subset, keep=keep)
        res = ddf.drop_duplicates(subset=subset, keep=keep,
                                  split_every=split_every, split_out=10)
        assert res.npartitions == 10
        assert_eq(sol, res)


@pytest.mark.parametrize('split_every', [None, 2])
def test_split_out_value_counts(split_every):
    df = pd.DataFrame({'x': [1, 2, 3] * 100})
    ddf = dd.from_pandas(df, npartitions=5)

    assert ddf.x.value_counts(split_out=10, split_every=split_every).npartitions == 10
    assert_eq(ddf.x.value_counts(split_out=10, split_every=split_every), df.x.value_counts())


def test_values():
    from dask.array.utils import assert_eq
    df = pd.DataFrame({'x': ['a', 'b', 'c', 'd'],
                       'y': [2, 3, 4, 5]},
                      index=pd.Index([1., 2., 3., 4.], name='ind'))
    ddf = dd.from_pandas(df, 2)

    assert_eq(df.values, ddf.values)
    assert_eq(df.x.values, ddf.x.values)
    assert_eq(df.y.values, ddf.y.values)
    assert_eq(df.index.values, ddf.index.values)


def test_copy():
    df = pd.DataFrame({'x': [1, 2, 3]})

    a = dd.from_pandas(df, npartitions=2)
    b = a.copy()

    a['y'] = a.x * 2

    assert_eq(b, df)

    df['y'] = df.x * 2


def test_del():
    df = pd.DataFrame({'x': ['a', 'b', 'c', 'd'],
                       'y': [2, 3, 4, 5]},
                      index=pd.Index([1., 2., 3., 4.], name='ind'))
    a = dd.from_pandas(df, 2)
    b = a.copy()

    del a['x']
    assert_eq(b, df)

    del df['x']
    assert_eq(a, df)


@pytest.mark.parametrize('index', [True, False])
@pytest.mark.parametrize('deep', [True, False])
def test_memory_usage(index, deep):
    df = pd.DataFrame({'x': [1, 2, 3],
                       'y': [1.0, 2.0, 3.0],
                       'z': ['a', 'b', 'c']})
    ddf = dd.from_pandas(df, npartitions=2)

    assert_eq(df.memory_usage(index=index, deep=deep),
              ddf.memory_usage(index=index, deep=deep))
    assert (df.x.memory_usage(index=index, deep=deep) ==
            ddf.x.memory_usage(index=index, deep=deep).compute())


@pytest.mark.parametrize('reduction', ['sum', 'mean', 'std', 'var', 'count',
                                       'min', 'max', 'idxmin', 'idxmax',
                                       'prod', 'all', 'sem'])
def test_dataframe_reductions_arithmetic(reduction):
    df = pd.DataFrame({'x': [1, 2, 3, 4, 5],
                       'y': [1.1, 2.2, 3.3, 4.4, 5.5]})
    ddf = dd.from_pandas(df, npartitions=3)

    assert_eq(ddf - (getattr(ddf, reduction)() + 1),
              df - (getattr(df, reduction)() + 1))


def test_datetime_loc_open_slicing():
    dtRange = pd.date_range('01.01.2015','05.05.2015')
    df = pd.DataFrame(np.random.random((len(dtRange), 2)), index=dtRange)
    ddf = dd.from_pandas(df, npartitions=5)
    assert_eq(df.loc[:'02.02.2015'], ddf.loc[:'02.02.2015'])
    assert_eq(df.loc['02.02.2015':], ddf.loc['02.02.2015':])
    assert_eq(df[0].loc[:'02.02.2015'], ddf[0].loc[:'02.02.2015'])
    assert_eq(df[0].loc['02.02.2015':], ddf[0].loc['02.02.2015':])


def test_to_datetime():
    df = pd.DataFrame({'year': [2015, 2016],
                       'month': [2, 3],
                       'day': [4, 5]})
    ddf = dd.from_pandas(df, npartitions=2)

    assert_eq(pd.to_datetime(df), dd.to_datetime(ddf))

    s = pd.Series(['3/11/2000', '3/12/2000', '3/13/2000'] * 100)
    ds = dd.from_pandas(s, npartitions=10)

    assert_eq(pd.to_datetime(s, infer_datetime_format=True),
              dd.to_datetime(ds, infer_datetime_format=True))


def test_to_timedelta():
    s = pd.Series(range(10))
    ds = dd.from_pandas(s, npartitions=2)

    assert_eq(pd.to_timedelta(s), dd.to_timedelta(ds))
    assert_eq(pd.to_timedelta(s, unit='h'), dd.to_timedelta(ds, unit='h'))

    s = pd.Series([1, 2, 'this will error'])
    ds = dd.from_pandas(s, npartitions=2)
    assert_eq(pd.to_timedelta(s, errors='coerce'),
              dd.to_timedelta(ds, errors='coerce'))


@pytest.mark.skipif(PANDAS_VERSION < '0.22.0',
                    reason="No isna method")
@pytest.mark.parametrize('values', [[np.NaN, 0], [1, 1]])
def test_isna(values):
    s = pd.Series(values)
    ds = dd.from_pandas(s, npartitions=2)

    assert_eq(pd.isna(s), dd.isna(ds))


@pytest.mark.parametrize('drop', [0, 9])
def test_slice_on_filtered_boundary(drop):
    # https://github.com/dask/dask/issues/2211
    x = np.arange(10)
    x[[5, 6]] -= 2
    df = pd.DataFrame({"A": x, "B": np.arange(len(x))})
    pdf = df.set_index("A").query("B != {}".format(drop))
    ddf = dd.from_pandas(df, 1).set_index("A").query("B != {}".format(drop))

    result = dd.concat([ddf, ddf.rename(columns={"B": "C"})], axis=1)
    expected = pd.concat([pdf, pdf.rename(columns={"B": "C"})], axis=1)
    assert_eq(result, expected)


def test_boundary_slice_nonmonotonic():
    x = np.array([-1, -2, 2, 4, 3])
    df = pd.DataFrame({"B": range(len(x))}, index=x)
    result = methods.boundary_slice(df, 0, 4)
    expected = df.iloc[2:]
    tm.assert_frame_equal(result, expected)

    result = methods.boundary_slice(df, -1, 4)
    expected = df.drop(-2)
    tm.assert_frame_equal(result, expected)

    result = methods.boundary_slice(df, -2, 3)
    expected = df.drop(4)
    tm.assert_frame_equal(result, expected)

    result = methods.boundary_slice(df, -2, 3.5)
    expected = df.drop(4)
    tm.assert_frame_equal(result, expected)

    result = methods.boundary_slice(df, -2, 4)
    expected = df
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize('start, stop, right_boundary, left_boundary, drop', [
    (-1, None, False, False, [-1, -2]),
    (-1, None, False, True, [-2]),
    (None, 3, False, False, [3, 4]),
    (None, 3, True, False, [4]),
    # Missing keys
    (-.5, None, False, False, [-1, -2]),
    (-.5, None, False, True, [-1, -2]),
    (-1.5, None, False, True, [-2]),
    (None, 3.5, False, False, [4]),
    (None, 3.5, True, False, [4]),
    (None, 2.5, False, False, [3, 4]),
])
def test_with_boundary(start, stop, right_boundary, left_boundary, drop):
    x = np.array([-1, -2, 2, 4, 3])
    df = pd.DataFrame({"B": range(len(x))}, index=x)
    result = methods.boundary_slice(df, start, stop, right_boundary, left_boundary)
    expected = df.drop(drop)
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize('index, left, right', [
    (range(10), 0, 9),
    (range(10), -1, None),
    (range(10), None, 10),
    ([-1, 0, 2, 1], None, None),
    ([-1, 0, 2, 1], -1, None),
    ([-1, 0, 2, 1], None, 2),
    ([-1, 0, 2, 1], -2, 3),
    (pd.date_range("2017", periods=10), None, None),
    (pd.date_range("2017", periods=10), pd.Timestamp("2017"), None),
    (pd.date_range("2017", periods=10), None, pd.Timestamp("2017-01-10")),
    (pd.date_range("2017", periods=10), pd.Timestamp("2016"), None),
    (pd.date_range("2017", periods=10), None, pd.Timestamp("2018")),
])
def test_boundary_slice_same(index, left, right):
    df = pd.DataFrame({"A": range(len(index))}, index=index)
    result = methods.boundary_slice(df, left, right)
    tm.assert_frame_equal(result, df)


def test_better_errors_object_reductions():
    # GH2452
    s = pd.Series(['a', 'b', 'c', 'd'])
    ds = dd.from_pandas(s, npartitions=2)
    with pytest.raises(ValueError) as err:
        ds.mean()
    assert str(err.value) == "`mean` not supported with object series"


def test_sample_empty_partitions():
    @dask.delayed
    def make_df(n):
        return pd.DataFrame(np.zeros((n, 4)), columns=list('abcd'))
    ddf = dd.from_delayed([make_df(0), make_df(100), make_df(0)])
    ddf2 = ddf.sample(frac=0.2)
    # smoke test sample on empty partitions
    res = ddf2.compute()
    assert res.dtypes.equals(ddf2.dtypes)


def test_coerce():
    df = pd.DataFrame(np.arange(100).reshape((10,10)))
    ddf = dd.from_pandas(df, npartitions=2)
    funcs = (int, float, complex)
    for d,t in product(funcs,(ddf, ddf[0])):
        pytest.raises(TypeError, lambda: t(d))


def test_bool():
    df = pd.DataFrame(np.arange(100).reshape((10,10)))
    ddf = dd.from_pandas(df, npartitions=2)
    conditions = [ddf, ddf[0], ddf == ddf, ddf[0] == ddf[0]]
    for cond in conditions:
        with pytest.raises(ValueError):
            bool(cond)


def test_cumulative_multiple_columns():
    # GH 3037
    df = pd.DataFrame(np.random.randn(100, 5), columns=list('abcde'))
    ddf = dd.from_pandas(df, 5)

    for d in [ddf, df]:
        for c in df.columns:
            d[c + 'cs'] = d[c].cumsum()
            d[c + 'cmin'] = d[c].cummin()
            d[c + 'cmax'] = d[c].cummax()
            d[c + 'cp'] = d[c].cumprod()

    assert_eq(ddf, df)


@pytest.mark.parametrize('func', [np.asarray, M.to_records])
def test_map_partition_array(func):
    from dask.array.utils import assert_eq
    df = pd.DataFrame({'x': [1, 2, 3, 4, 5],
                       'y': [6.0, 7.0, 8.0, 9.0, 10.0]},
                      index=['a', 'b', 'c', 'd', 'e'])
    ddf = dd.from_pandas(df, npartitions=2)

    for pre in [lambda a: a,
                lambda a: a.x,
                lambda a: a.y,
                lambda a: a.index]:

        try:
            expected = func(pre(df))
        except Exception:
            continue
        x = pre(ddf).map_partitions(func)
        assert_eq(x, expected)

        assert isinstance(x, da.Array)
        assert x.chunks[0] == (np.nan, np.nan)


def test_map_partition_sparse():
    sparse = pytest.importorskip('sparse')
    df = pd.DataFrame({'x': [1, 2, 3, 4, 5],
                       'y': [6.0, 7.0, 8.0, 9.0, 10.0]},
                      index=['a', 'b', 'c', 'd', 'e'])
    ddf = dd.from_pandas(df, npartitions=2)

    def f(d):
        return sparse.COO(np.array(d))

    for pre in [lambda a: a,
                lambda a: a.x]:
        expected = f(pre(df))
        result = pre(ddf).map_partitions(f)
        assert isinstance(result, da.Array)
        computed = result.compute()
        assert (computed.data == expected.data).all()
        assert (computed.coords == expected.coords).all()


def test_mixed_dask_array_operations():
    df = pd.DataFrame({'x': [1, 2, 3]}, index=[4, 5, 6])
    ddf = dd.from_pandas(df, npartitions=2)

    assert_eq(df.x + df.x.values,
              ddf.x + ddf.x.values)
    assert_eq(df.x.values + df.x,
              ddf.x.values + ddf.x)

    assert_eq(df.x + df.index.values,
              ddf.x + ddf.index.values)
    assert_eq(df.index.values + df.x,
              ddf.index.values + ddf.x)


def test_mixed_dask_array_operations_errors():
    df = pd.DataFrame({'x': [1, 2, 3, 4, 5]}, index=[4, 5, 6, 7, 8])
    ddf = dd.from_pandas(df, npartitions=2)

    x = da.arange(5, chunks=((1, 4),))
    x._chunks = ((np.nan, np.nan),)

    with pytest.raises(ValueError):
        (ddf.x + x).compute()

    x = da.arange(5, chunks=((2, 2, 1),))
    with pytest.raises(ValueError) as info:
        ddf.x + x

    assert 'add' in str(info.value)


def test_mixed_dask_array_multi_dimensional():
    df = pd.DataFrame({'x': [1, 2, 3, 4, 5],
                       'y': [5., 6., 7., 8., 9.]},
                      columns=['x', 'y'])
    ddf = dd.from_pandas(df, npartitions=2)

    x = (df.values + 1).astype(float)
    dx = (ddf.values + 1).astype(float)

    assert_eq(ddf + dx + 1, df + x + 1)
    assert_eq(ddf + dx.rechunk((None, 1)) + 1, df + x + 1)
    assert_eq(ddf[['y', 'x']] + dx + 1, df[['y', 'x']] + x + 1)


def test_meta_raises():
    # Raise when we use a user defined fucntion
    s = pd.Series(['abcd', 'abcd'])
    ds = dd.from_pandas(s, npartitions=2)
    try:
        ds.map(lambda x: x[3])
    except ValueError as e:
        assert "meta=" in str(e)

    # But not otherwise
    df = pd.DataFrame({'a': ['x', 'y', 'y'],
                       'b': ['x', 'y', 'z'],
                       'c': [1, 2, 3]})
    ddf = dd.from_pandas(df, npartitions=1)

    with pytest.raises(Exception) as info:
        ddf.a + ddf.c

    assert "meta=" not in str(info.value)
