import collections

import numpy as np
import pandas as pd
import pandas.util.testing as tm

import pytest

import dask
import dask.dataframe as dd
from dask.dataframe.utils import assert_eq, assert_dask_graph, assert_max_deps, PANDAS_VERSION


AGG_FUNCS = ['sum', 'mean', 'min', 'max', 'count', 'size', 'std', 'var', 'nunique', 'first', 'last']


@pytest.fixture(params=AGG_FUNCS)
def agg_func(request):
    """
    Aggregations supported for groups
    """
    return request.param


def groupby_internal_repr():
    pdf = pd.DataFrame({'x': [1, 2, 3, 4, 6, 7, 8, 9, 10],
                        'y': list('abcbabbcda')})
    ddf = dd.from_pandas(pdf, 3)

    gp = pdf.groupby('y')
    dp = ddf.groupby('y')
    assert isinstance(dp, dd.groupby.DataFrameGroupBy)
    assert isinstance(dp._meta, pd.core.groupby.DataFrameGroupBy)
    assert isinstance(dp.obj, dd.DataFrame)
    assert_eq(dp.obj, gp.obj)

    gp = pdf.groupby('y')['x']
    dp = ddf.groupby('y')['x']
    assert isinstance(dp, dd.groupby.SeriesGroupBy)
    assert isinstance(dp._meta, pd.core.groupby.SeriesGroupBy)
    # slicing should not affect to internal
    assert isinstance(dp.obj, dd.Series)
    assert_eq(dp.obj, gp.obj)

    gp = pdf.groupby('y')[['x']]
    dp = ddf.groupby('y')[['x']]
    assert isinstance(dp, dd.groupby.DataFrameGroupBy)
    assert isinstance(dp._meta, pd.core.groupby.DataFrameGroupBy)
    # slicing should not affect to internal
    assert isinstance(dp.obj, dd.DataFrame)
    assert_eq(dp.obj, gp.obj)

    gp = pdf.groupby(pdf.y)['x']
    dp = ddf.groupby(ddf.y)['x']
    assert isinstance(dp, dd.groupby.SeriesGroupBy)
    assert isinstance(dp._meta, pd.core.groupby.SeriesGroupBy)
    # slicing should not affect to internal
    assert isinstance(dp.obj, dd.Series)
    assert_eq(dp.obj, gp.obj)

    gp = pdf.groupby(pdf.y)[['x']]
    dp = ddf.groupby(ddf.y)[['x']]
    assert isinstance(dp, dd.groupby.DataFrameGroupBy)
    assert isinstance(dp._meta, pd.core.groupby.DataFrameGroupBy)
    # slicing should not affect to internal
    assert isinstance(dp.obj, dd.DataFrame)
    assert_eq(dp.obj, gp.obj)


def groupby_error():
    pdf = pd.DataFrame({'x': [1, 2, 3, 4, 6, 7, 8, 9, 10],
                        'y': list('abcbabbcda')})
    ddf = dd.from_pandas(pdf, 3)

    with pytest.raises(KeyError):
        ddf.groupby('A')

    with pytest.raises(KeyError):
        ddf.groupby(['x', 'A'])

    dp = ddf.groupby('y')

    msg = 'Column not found: '
    with pytest.raises(KeyError) as err:
        dp['A']
    assert msg in str(err.value)

    with pytest.raises(KeyError) as err:
        dp[['x', 'A']]
    assert msg in str(err.value)


def groupby_internal_head():
    pdf = pd.DataFrame({'A': [1, 2] * 10,
                        'B': np.random.randn(20),
                        'C': np.random.randn(20)})
    ddf = dd.from_pandas(pdf, 3)

    assert_eq(ddf.groupby('A')._head().sum(),
              pdf.head().groupby('A').sum())

    assert_eq(ddf.groupby(ddf['A'])._head().sum(),
              pdf.head().groupby(pdf['A']).sum())

    assert_eq(ddf.groupby(ddf['A'] + 1)._head().sum(),
              pdf.head().groupby(pdf['A'] + 1).sum())


def test_full_groupby():
    df = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7, 8, 9],
                       'b': [4, 5, 6, 3, 2, 1, 0, 0, 0]},
                      index=[0, 1, 3, 5, 6, 8, 9, 9, 9])
    ddf = dd.from_pandas(df, npartitions=3)

    pytest.raises(KeyError, lambda: ddf.groupby('does_not_exist'))
    pytest.raises(AttributeError, lambda: ddf.groupby('a').does_not_exist)
    assert 'b' in dir(ddf.groupby('a'))

    def func(df):
        return df.assign(b=df.b - df.b.mean())

    assert_eq(df.groupby('a').apply(func),
              ddf.groupby('a').apply(func))


def test_full_groupby_apply_multiarg():
    df = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7, 8, 9],
                       'b': [4, 5, 6, 3, 2, 1, 0, 0, 0]},
                      index=[0, 1, 3, 5, 6, 8, 9, 9, 9])
    ddf = dd.from_pandas(df, npartitions=3)

    def func(df, c, d=3):
        return df.assign(b=df.b - df.b.mean() + c * d)

    c = df.a.sum()
    d = df.b.mean()

    c_scalar = ddf.a.sum()
    d_scalar = ddf.b.mean()
    c_delayed = dask.delayed(lambda: c)()
    d_delayed = dask.delayed(lambda: d)()

    meta = df.groupby('a').apply(func, c)

    for c_lazy, d_lazy in [(c_scalar, d_scalar),
                           (c_delayed, d_delayed)]:
        assert_eq(df.groupby('a').apply(func, c),
                  ddf.groupby('a').apply(func, c))

        assert_eq(df.groupby('a').apply(func, c, d=d),
                  ddf.groupby('a').apply(func, c, d=d))

        assert_eq(df.groupby('a').apply(func, c),
                  ddf.groupby('a').apply(func, c_lazy), check_dtype=False)

        assert_eq(df.groupby('a').apply(func, c),
                  ddf.groupby('a').apply(func, c_lazy, meta=meta))

        assert_eq(df.groupby('a').apply(func, c, d=d),
                  ddf.groupby('a').apply(func, c, d=d_lazy))

        assert_eq(df.groupby('a').apply(func, c, d=d),
                  ddf.groupby('a').apply(func, c, d=d_lazy, meta=meta))


@pytest.mark.parametrize('grouper', [
    lambda df: ['a'],
    lambda df: ['a', 'b'],
    lambda df: df['a'],
    lambda df: [df['a'], df['b']],
    pytest.mark.xfail(reason="not yet supported")(lambda df: [df['a'] > 2, df['b'] > 1])
])
@pytest.mark.parametrize('reverse', [True, False])
def test_full_groupby_multilevel(grouper, reverse):
    index = [0, 1, 3, 5, 6, 8, 9, 9, 9]
    if reverse:
        index = index[::-1]
    df = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7, 8, 9],
                       'd': [1, 2, 3, 4, 5, 6, 7, 8, 9],
                       'b': [4, 5, 6, 3, 2, 1, 0, 0, 0]},
                      index=index)
    ddf = dd.from_pandas(df, npartitions=3)

    def func(df):
        return df.assign(b=df.b - df.b.mean())

    # last one causes a DeprcationWarning from pandas.
    # See https://github.com/pandas-dev/pandas/issues/16481
    assert_eq(df.groupby(grouper(df)).apply(func),
              ddf.groupby(grouper(ddf)).apply(func))


def test_groupby_dir():
    df = pd.DataFrame({'a': range(10), 'b c d e': range(10)})
    ddf = dd.from_pandas(df, npartitions=2)
    g = ddf.groupby('a')
    assert 'a' in dir(g)
    assert 'b c d e' not in dir(g)


@pytest.mark.parametrize('get', [dask.get, dask.threaded.get])
def test_groupby_on_index(get):
    pdf = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7, 8, 9],
                        'b': [4, 5, 6, 3, 2, 1, 0, 0, 0]},
                       index=[0, 1, 3, 5, 6, 8, 9, 9, 9])
    ddf = dd.from_pandas(pdf, npartitions=3)

    ddf2 = ddf.set_index('a')
    pdf2 = pdf.set_index('a')
    assert_eq(ddf.groupby('a').b.mean(), ddf2.groupby(ddf2.index).b.mean())

    def func(df):
        return df.assign(b=df.b - df.b.mean())

    def func2(df):
        return df[['b']] - df[['b']].mean()

    with dask.set_options(get=get):
        with pytest.warns(None):
            assert_eq(ddf.groupby('a').apply(func),
                      pdf.groupby('a').apply(func))

            assert_eq(ddf.groupby('a').apply(func).set_index('a'),
                      pdf.groupby('a').apply(func).set_index('a'))

            assert_eq(pdf2.groupby(pdf2.index).apply(func2),
                      ddf2.groupby(ddf2.index).apply(func2))


@pytest.mark.parametrize('grouper',
                         [lambda df: df.groupby('a')['b'],
                          lambda df: df.groupby(['a', 'b']),
                          lambda df: df.groupby(['a', 'b'])['c'],
                          lambda df: df.groupby(df['a'])[['b', 'c']],
                          lambda df: df.groupby('a')[['b', 'c']],
                          lambda df: df.groupby('a')[['b']],
                          lambda df: df.groupby(['a', 'b', 'c'])])
def test_groupby_multilevel_getitem(grouper, agg_func):
    # nunique is not implemented for DataFrameGroupBy
    if agg_func == 'nunique':
        return

    df = pd.DataFrame({'a': [1, 2, 3, 1, 2, 3],
                       'b': [1, 2, 1, 4, 2, 1],
                       'c': [1, 3, 2, 1, 1, 2],
                       'd': [1, 2, 1, 1, 2, 2]})
    ddf = dd.from_pandas(df, 2)

    dask_group = grouper(ddf)
    pandas_group = grouper(df)

    dask_agg = getattr(dask_group, agg_func)
    pandas_agg = getattr(pandas_group, agg_func)

    assert isinstance(dask_group, dd.groupby._GroupBy)
    assert isinstance(pandas_group, pd.core.groupby.GroupBy)

    if agg_func == 'mean':
        assert_eq(dask_agg(), pandas_agg().astype(float))
    else:
        assert_eq(dask_agg(), pandas_agg())


def test_groupby_multilevel_agg():
    df = pd.DataFrame({'a': [1, 2, 3, 1, 2, 3],
                       'b': [1, 2, 1, 4, 2, 1],
                       'c': [1, 3, 2, 1, 1, 2],
                       'd': [1, 2, 1, 1, 2, 2]})
    ddf = dd.from_pandas(df, 2)

    sol = df.groupby(['a']).mean()
    res = ddf.groupby(['a']).mean()
    assert_eq(res, sol)

    sol = df.groupby(['a', 'c']).mean()
    res = ddf.groupby(['a', 'c']).mean()
    assert_eq(res, sol)

    sol = df.groupby([df['a'], df['c']]).mean()
    res = ddf.groupby([ddf['a'], ddf['c']]).mean()
    assert_eq(res, sol)


def test_groupby_get_group():
    dsk = {('x', 0): pd.DataFrame({'a': [1, 2, 6], 'b': [4, 2, 7]},
                                  index=[0, 1, 3]),
           ('x', 1): pd.DataFrame({'a': [4, 2, 6], 'b': [3, 3, 1]},
                                  index=[5, 6, 8]),
           ('x', 2): pd.DataFrame({'a': [4, 3, 7], 'b': [1, 1, 3]},
                                  index=[9, 9, 9])}
    meta = dsk[('x', 0)]
    d = dd.DataFrame(dsk, 'x', meta, [0, 4, 9, 9])
    full = d.compute()

    for ddkey, pdkey in [('b', 'b'), (d.b, full.b),
                         (d.b + 1, full.b + 1)]:
        ddgrouped = d.groupby(ddkey)
        pdgrouped = full.groupby(pdkey)
        # DataFrame
        assert_eq(ddgrouped.get_group(2), pdgrouped.get_group(2))
        assert_eq(ddgrouped.get_group(3), pdgrouped.get_group(3))
        # Series
        assert_eq(ddgrouped.a.get_group(3), pdgrouped.a.get_group(3))
        assert_eq(ddgrouped.a.get_group(2), pdgrouped.a.get_group(2))


def test_dataframe_groupby_nunique():
    strings = list('aaabbccccdddeee')
    data = np.random.randn(len(strings))
    ps = pd.DataFrame(dict(strings=strings, data=data))
    s = dd.from_pandas(ps, npartitions=3)
    expected = ps.groupby('strings')['data'].nunique()
    assert_eq(s.groupby('strings')['data'].nunique(), expected)


def test_dataframe_groupby_nunique_across_group_same_value():
    strings = list('aaabbccccdddeee')
    data = list(map(int, '123111223323412'))
    ps = pd.DataFrame(dict(strings=strings, data=data))
    s = dd.from_pandas(ps, npartitions=3)
    expected = ps.groupby('strings')['data'].nunique()
    assert_eq(s.groupby('strings')['data'].nunique(), expected)


def test_series_groupby_propagates_names():
    df = pd.DataFrame({'x': [1, 2, 3], 'y': [4, 5, 6]})
    ddf = dd.from_pandas(df, 2)
    func = lambda df: df['y'].sum()
    with pytest.warns(UserWarning):  # meta inference
        result = ddf.groupby('x').apply(func)
    expected = df.groupby('x').apply(func)
    assert_eq(result, expected)


def test_series_groupby():
    s = pd.Series([1, 2, 2, 1, 1])
    pd_group = s.groupby(s)

    ss = dd.from_pandas(s, npartitions=2)
    dask_group = ss.groupby(ss)

    pd_group2 = s.groupby(s + 1)
    dask_group2 = ss.groupby(ss + 1)

    for dg, pdg in [(dask_group, pd_group), (pd_group2, dask_group2)]:
        assert_eq(dg.count(), pdg.count())
        assert_eq(dg.sum(), pdg.sum())
        assert_eq(dg.min(), pdg.min())
        assert_eq(dg.max(), pdg.max())
        assert_eq(dg.size(), pdg.size())
        assert_eq(dg.first(), pdg.first())
        assert_eq(dg.last(), pdg.last())


def test_series_groupby_errors():
    s = pd.Series([1, 2, 2, 1, 1])

    ss = dd.from_pandas(s, npartitions=2)

    msg = "No group keys passed!"
    with pytest.raises(ValueError) as err:
        s.groupby([])    # pandas
    assert msg in str(err.value)
    with pytest.raises(ValueError) as err:
        ss.groupby([])   # dask should raise the same error
    assert msg in str(err.value)

    sss = dd.from_pandas(s, npartitions=3)
    pytest.raises(NotImplementedError, lambda: ss.groupby(sss))

    with pytest.raises(KeyError):
        s.groupby('x')    # pandas
    with pytest.raises(KeyError):
        ss.groupby('x')   # dask should raise the same error


def test_groupby_index_array():
    df = tm.makeTimeDataFrame()
    ddf = dd.from_pandas(df, npartitions=2)

    # first select column, then group
    assert_eq(df.A.groupby(df.index.month).nunique(),
              ddf.A.groupby(ddf.index.month).nunique(), check_names=False)

    # first group, then select column
    assert_eq(df.groupby(df.index.month).A.nunique(),
              ddf.groupby(ddf.index.month).A.nunique(), check_names=False)


def test_groupby_set_index():
    df = tm.makeTimeDataFrame()
    ddf = dd.from_pandas(df, npartitions=2)
    pytest.raises(TypeError,
                  lambda: ddf.groupby(df.index.month, as_index=False))


def test_split_apply_combine_on_series():
    pdf = pd.DataFrame({'a': [1, 2, 6, 4, 4, 6, 4, 3, 7],
                        'b': [4, 2, 7, 3, 3, 1, 1, 1, 2]},
                       index=[0, 1, 3, 5, 6, 8, 9, 9, 9])
    ddf = dd.from_pandas(pdf, npartitions=3)

    for ddkey, pdkey in [('b', 'b'), (ddf.b, pdf.b), (ddf.b + 1, pdf.b + 1)]:
        assert_eq(ddf.groupby(ddkey).a.min(), pdf.groupby(pdkey).a.min())
        assert_eq(ddf.groupby(ddkey).a.max(), pdf.groupby(pdkey).a.max())
        assert_eq(ddf.groupby(ddkey).a.count(), pdf.groupby(pdkey).a.count())
        assert_eq(ddf.groupby(ddkey).a.mean(), pdf.groupby(pdkey).a.mean())
        assert_eq(ddf.groupby(ddkey).a.nunique(), pdf.groupby(pdkey).a.nunique())
        assert_eq(ddf.groupby(ddkey).a.size(), pdf.groupby(pdkey).a.size())
        assert_eq(ddf.groupby(ddkey).a.first(), pdf.groupby(pdkey).a.first())
        assert_eq(ddf.groupby(ddkey).a.last(), pdf.groupby(pdkey).a.last())
        for ddof in [0, 1, 2]:
            assert_eq(ddf.groupby(ddkey).a.var(ddof),
                      pdf.groupby(pdkey).a.var(ddof))
            assert_eq(ddf.groupby(ddkey).a.std(ddof),
                      pdf.groupby(pdkey).a.std(ddof))

        assert_eq(ddf.groupby(ddkey).sum(), pdf.groupby(pdkey).sum())
        assert_eq(ddf.groupby(ddkey).min(), pdf.groupby(pdkey).min())
        assert_eq(ddf.groupby(ddkey).max(), pdf.groupby(pdkey).max())
        assert_eq(ddf.groupby(ddkey).count(), pdf.groupby(pdkey).count())
        assert_eq(ddf.groupby(ddkey).mean(), pdf.groupby(pdkey).mean())
        assert_eq(ddf.groupby(ddkey).size(), pdf.groupby(pdkey).size())
        assert_eq(ddf.groupby(ddkey).first(), pdf.groupby(pdkey).first())
        assert_eq(ddf.groupby(ddkey).last(), pdf.groupby(pdkey).last())

        for ddof in [0, 1, 2]:
            assert_eq(ddf.groupby(ddkey).var(ddof),
                      pdf.groupby(pdkey).var(ddof), check_dtype=False)
            assert_eq(ddf.groupby(ddkey).std(ddof),
                      pdf.groupby(pdkey).std(ddof), check_dtype=False)

    for ddkey, pdkey in [(ddf.b, pdf.b), (ddf.b + 1, pdf.b + 1)]:
        assert_eq(ddf.a.groupby(ddkey).sum(), pdf.a.groupby(pdkey).sum(), check_names=False)
        assert_eq(ddf.a.groupby(ddkey).max(), pdf.a.groupby(pdkey).max(), check_names=False)
        assert_eq(ddf.a.groupby(ddkey).count(), pdf.a.groupby(pdkey).count(), check_names=False)
        assert_eq(ddf.a.groupby(ddkey).mean(), pdf.a.groupby(pdkey).mean(), check_names=False)
        assert_eq(ddf.a.groupby(ddkey).nunique(), pdf.a.groupby(pdkey).nunique(), check_names=False)
        assert_eq(ddf.a.groupby(ddkey).first(), pdf.a.groupby(pdkey).first(), check_names=False)
        assert_eq(ddf.a.groupby(ddkey).last(), pdf.a.groupby(pdkey).last(), check_names=False)

        for ddof in [0, 1, 2]:
            assert_eq(ddf.a.groupby(ddkey).var(ddof),
                      pdf.a.groupby(pdkey).var(ddof))
            assert_eq(ddf.a.groupby(ddkey).std(ddof),
                      pdf.a.groupby(pdkey).std(ddof))

    for i in [0, 4, 7]:
        assert_eq(ddf.groupby(ddf.b > i).a.sum(), pdf.groupby(pdf.b > i).a.sum())
        assert_eq(ddf.groupby(ddf.b > i).a.min(), pdf.groupby(pdf.b > i).a.min())
        assert_eq(ddf.groupby(ddf.b > i).a.max(), pdf.groupby(pdf.b > i).a.max())
        assert_eq(ddf.groupby(ddf.b > i).a.count(), pdf.groupby(pdf.b > i).a.count())
        assert_eq(ddf.groupby(ddf.b > i).a.mean(), pdf.groupby(pdf.b > i).a.mean())
        assert_eq(ddf.groupby(ddf.b > i).a.nunique(), pdf.groupby(pdf.b > i).a.nunique())
        assert_eq(ddf.groupby(ddf.b > i).a.size(), pdf.groupby(pdf.b > i).a.size())
        assert_eq(ddf.groupby(ddf.b > i).a.first(), pdf.groupby(pdf.b > i).a.first())
        assert_eq(ddf.groupby(ddf.b > i).a.last(), pdf.groupby(pdf.b > i).a.last())

        assert_eq(ddf.groupby(ddf.a > i).b.sum(), pdf.groupby(pdf.a > i).b.sum())
        assert_eq(ddf.groupby(ddf.a > i).b.min(), pdf.groupby(pdf.a > i).b.min())
        assert_eq(ddf.groupby(ddf.a > i).b.max(), pdf.groupby(pdf.a > i).b.max())
        assert_eq(ddf.groupby(ddf.a > i).b.count(), pdf.groupby(pdf.a > i).b.count())
        assert_eq(ddf.groupby(ddf.a > i).b.mean(), pdf.groupby(pdf.a > i).b.mean())
        assert_eq(ddf.groupby(ddf.a > i).b.nunique(), pdf.groupby(pdf.a > i).b.nunique())
        assert_eq(ddf.groupby(ddf.b > i).b.size(), pdf.groupby(pdf.b > i).b.size())
        assert_eq(ddf.groupby(ddf.b > i).b.first(), pdf.groupby(pdf.b > i).b.first())
        assert_eq(ddf.groupby(ddf.b > i).b.last(), pdf.groupby(pdf.b > i).b.last())

        assert_eq(ddf.groupby(ddf.b > i).sum(), pdf.groupby(pdf.b > i).sum())
        assert_eq(ddf.groupby(ddf.b > i).min(), pdf.groupby(pdf.b > i).min())
        assert_eq(ddf.groupby(ddf.b > i).max(), pdf.groupby(pdf.b > i).max())
        assert_eq(ddf.groupby(ddf.b > i).count(), pdf.groupby(pdf.b > i).count())
        assert_eq(ddf.groupby(ddf.b > i).mean(), pdf.groupby(pdf.b > i).mean())
        assert_eq(ddf.groupby(ddf.b > i).size(), pdf.groupby(pdf.b > i).size())
        assert_eq(ddf.groupby(ddf.b > i).first(), pdf.groupby(pdf.b > i).first())
        assert_eq(ddf.groupby(ddf.b > i).last(), pdf.groupby(pdf.b > i).last())

        assert_eq(ddf.groupby(ddf.a > i).sum(), pdf.groupby(pdf.a > i).sum())
        assert_eq(ddf.groupby(ddf.a > i).min(), pdf.groupby(pdf.a > i).min())
        assert_eq(ddf.groupby(ddf.a > i).max(), pdf.groupby(pdf.a > i).max())
        assert_eq(ddf.groupby(ddf.a > i).count(), pdf.groupby(pdf.a > i).count())
        assert_eq(ddf.groupby(ddf.a > i).mean(), pdf.groupby(pdf.a > i).mean())
        assert_eq(ddf.groupby(ddf.a > i).size(), pdf.groupby(pdf.a > i).size())
        assert_eq(ddf.groupby(ddf.a > i).first(), pdf.groupby(pdf.a > i).first())
        assert_eq(ddf.groupby(ddf.a > i).last(), pdf.groupby(pdf.a > i).last())

        for ddof in [0, 1, 2]:
            assert_eq(ddf.groupby(ddf.b > i).std(ddof),
                      pdf.groupby(pdf.b > i).std(ddof))

    for ddkey, pdkey in [('a', 'a'), (ddf.a, pdf.a),
                         (ddf.a + 1, pdf.a + 1), (ddf.a > 3, pdf.a > 3)]:
        assert_eq(ddf.groupby(ddkey).b.sum(), pdf.groupby(pdkey).b.sum())
        assert_eq(ddf.groupby(ddkey).b.min(), pdf.groupby(pdkey).b.min())
        assert_eq(ddf.groupby(ddkey).b.max(), pdf.groupby(pdkey).b.max())
        assert_eq(ddf.groupby(ddkey).b.count(), pdf.groupby(pdkey).b.count())
        assert_eq(ddf.groupby(ddkey).b.mean(), pdf.groupby(pdkey).b.mean())
        assert_eq(ddf.groupby(ddkey).b.nunique(), pdf.groupby(pdkey).b.nunique())
        assert_eq(ddf.groupby(ddkey).b.size(), pdf.groupby(pdkey).b.size())
        assert_eq(ddf.groupby(ddkey).b.first(), pdf.groupby(pdkey).b.first())
        assert_eq(ddf.groupby(ddkey).last(), pdf.groupby(pdkey).last())

        assert_eq(ddf.groupby(ddkey).sum(), pdf.groupby(pdkey).sum())
        assert_eq(ddf.groupby(ddkey).min(), pdf.groupby(pdkey).min())
        assert_eq(ddf.groupby(ddkey).max(), pdf.groupby(pdkey).max())
        assert_eq(ddf.groupby(ddkey).count(), pdf.groupby(pdkey).count())
        assert_eq(ddf.groupby(ddkey).mean(), pdf.groupby(pdkey).mean().astype(float))
        assert_eq(ddf.groupby(ddkey).size(), pdf.groupby(pdkey).size())
        assert_eq(ddf.groupby(ddkey).first(), pdf.groupby(pdkey).first())
        assert_eq(ddf.groupby(ddkey).last(), pdf.groupby(pdkey).last())

        for ddof in [0, 1, 2]:
            assert_eq(ddf.groupby(ddkey).b.std(ddof),
                      pdf.groupby(pdkey).b.std(ddof))

    assert (sorted(ddf.groupby('b').a.sum().dask) ==
            sorted(ddf.groupby('b').a.sum().dask))
    assert (sorted(ddf.groupby(ddf.a > 3).b.mean().dask) ==
            sorted(ddf.groupby(ddf.a > 3).b.mean().dask))

    # test raises with incorrect key
    pytest.raises(KeyError, lambda: ddf.groupby('x'))
    pytest.raises(KeyError, lambda: ddf.groupby(['a', 'x']))
    pytest.raises(KeyError, lambda: ddf.groupby('a')['x'])
    pytest.raises(KeyError, lambda: ddf.groupby('a')['b', 'x'])
    pytest.raises(KeyError, lambda: ddf.groupby('a')[['b', 'x']])

    # test graph node labels
    assert_dask_graph(ddf.groupby('b').a.sum(), 'series-groupby-sum')
    assert_dask_graph(ddf.groupby('b').a.min(), 'series-groupby-min')
    assert_dask_graph(ddf.groupby('b').a.max(), 'series-groupby-max')
    assert_dask_graph(ddf.groupby('b').a.count(), 'series-groupby-count')
    assert_dask_graph(ddf.groupby('b').a.var(), 'series-groupby-var')
    assert_dask_graph(ddf.groupby('b').a.first(), 'series-groupby-first')
    assert_dask_graph(ddf.groupby('b').a.last(), 'series-groupby-last')
    # mean consists from sum and count operations
    assert_dask_graph(ddf.groupby('b').a.mean(), 'series-groupby-sum')
    assert_dask_graph(ddf.groupby('b').a.mean(), 'series-groupby-count')
    assert_dask_graph(ddf.groupby('b').a.nunique(), 'series-groupby-nunique')
    assert_dask_graph(ddf.groupby('b').a.size(), 'series-groupby-size')

    assert_dask_graph(ddf.groupby('b').sum(), 'dataframe-groupby-sum')
    assert_dask_graph(ddf.groupby('b').min(), 'dataframe-groupby-min')
    assert_dask_graph(ddf.groupby('b').max(), 'dataframe-groupby-max')
    assert_dask_graph(ddf.groupby('b').count(), 'dataframe-groupby-count')
    assert_dask_graph(ddf.groupby('b').first(), 'dataframe-groupby-first')
    assert_dask_graph(ddf.groupby('b').last(), 'dataframe-groupby-last')
    # mean consists from sum and count operations
    assert_dask_graph(ddf.groupby('b').mean(), 'dataframe-groupby-sum')
    assert_dask_graph(ddf.groupby('b').mean(), 'dataframe-groupby-count')
    assert_dask_graph(ddf.groupby('b').size(), 'dataframe-groupby-size')


@pytest.mark.parametrize('keyword', ['split_every', 'split_out'])
def test_groupby_reduction_split(keyword):
    pdf = pd.DataFrame({'a': [1, 2, 6, 4, 4, 6, 4, 3, 7] * 100,
                        'b': [4, 2, 7, 3, 3, 1, 1, 1, 2] * 100})
    ddf = dd.from_pandas(pdf, npartitions=15)

    def call(g, m, **kwargs):
        return getattr(g, m)(**kwargs)

    # DataFrame
    for m in AGG_FUNCS:
        # nunique is not implemented for DataFrameGroupBy
        if m == 'nunique':
            continue
        res = call(ddf.groupby('b'), m, **{keyword: 2})
        sol = call(pdf.groupby('b'), m)
        assert_eq(res, sol)
        assert call(ddf.groupby('b'), m)._name != res._name

    res = call(ddf.groupby('b'), 'var', ddof=2, **{keyword: 2})
    sol = call(pdf.groupby('b'), 'var', ddof=2)
    assert_eq(res, sol)
    assert call(ddf.groupby('b'), 'var', ddof=2)._name != res._name

    # Series, post select
    for m in AGG_FUNCS:
        res = call(ddf.groupby('b').a, m, **{keyword: 2})
        sol = call(pdf.groupby('b').a, m)
        assert_eq(res, sol)
        assert call(ddf.groupby('b').a, m)._name != res._name

    res = call(ddf.groupby('b').a, 'var', ddof=2, **{keyword: 2})
    sol = call(pdf.groupby('b').a, 'var', ddof=2)
    assert_eq(res, sol)
    assert call(ddf.groupby('b').a, 'var', ddof=2)._name != res._name

    # Series, pre select
    for m in AGG_FUNCS:
        res = call(ddf.a.groupby(ddf.b), m, **{keyword: 2})
        sol = call(pdf.a.groupby(pdf.b), m)
        # There's a bug in pandas 0.18.0 with `pdf.a.groupby(pdf.b).count()`
        # not forwarding the series name. Skip name checks here for now.
        assert_eq(res, sol, check_names=False)
        assert call(ddf.a.groupby(ddf.b), m)._name != res._name

    res = call(ddf.a.groupby(ddf.b), 'var', ddof=2, **{keyword: 2})
    sol = call(pdf.a.groupby(pdf.b), 'var', ddof=2)
    assert_eq(res, sol)
    assert call(ddf.a.groupby(ddf.b), 'var', ddof=2)._name != res._name


def test_apply_shuffle():
    pdf = pd.DataFrame({'A': [1, 2, 3, 4] * 5,
                        'B': np.random.randn(20),
                        'C': np.random.randn(20),
                        'D': np.random.randn(20)})
    ddf = dd.from_pandas(pdf, 3)

    with pytest.warns(UserWarning):  # meta inference
        assert_eq(ddf.groupby('A').apply(lambda x: x.sum()),
                  pdf.groupby('A').apply(lambda x: x.sum()))

        assert_eq(ddf.groupby(ddf['A']).apply(lambda x: x.sum()),
                  pdf.groupby(pdf['A']).apply(lambda x: x.sum()))

        assert_eq(ddf.groupby(ddf['A'] + 1).apply(lambda x: x.sum()),
                  pdf.groupby(pdf['A'] + 1).apply(lambda x: x.sum()))

        # SeriesGroupBy
        assert_eq(ddf.groupby('A')['B'].apply(lambda x: x.sum()),
                  pdf.groupby('A')['B'].apply(lambda x: x.sum()))

        assert_eq(ddf.groupby(ddf['A'])['B'].apply(lambda x: x.sum()),
                  pdf.groupby(pdf['A'])['B'].apply(lambda x: x.sum()))

        assert_eq(ddf.groupby(ddf['A'] + 1)['B'].apply(lambda x: x.sum()),
                  pdf.groupby(pdf['A'] + 1)['B'].apply(lambda x: x.sum()))

        # DataFrameGroupBy with column slice
        assert_eq(ddf.groupby('A')[['B', 'C']].apply(lambda x: x.sum()),
                  pdf.groupby('A')[['B', 'C']].apply(lambda x: x.sum()))

        assert_eq(ddf.groupby(ddf['A'])[['B', 'C']].apply(lambda x: x.sum()),
                  pdf.groupby(pdf['A'])[['B', 'C']].apply(lambda x: x.sum()))

        assert_eq(ddf.groupby(ddf['A'] + 1)[['B', 'C']].apply(lambda x: x.sum()),
                  pdf.groupby(pdf['A'] + 1)[['B', 'C']].apply(lambda x: x.sum()))


@pytest.mark.parametrize('grouper', [
    lambda df: 'AA',
    lambda df: ['AA', 'AB'],
    lambda df: df['AA'],
    lambda df: [df['AA'], df['AB']],
    lambda df: df['AA'] + 1,
    pytest.mark.xfail("NotImplemented")(lambda df: [df['AA'] + 1, df['AB'] + 1]),
])
def test_apply_shuffle_multilevel(grouper):
    pdf = pd.DataFrame({'AB': [1, 2, 3, 4] * 5,
                        'AA': [1, 2, 3, 4] * 5,
                        'B': np.random.randn(20),
                        'C': np.random.randn(20),
                        'D': np.random.randn(20)})
    ddf = dd.from_pandas(pdf, 3)

    with pytest.warns(UserWarning):
        # DataFrameGroupBy
        assert_eq(ddf.groupby(grouper(ddf)).apply(lambda x: x.sum()),
                  pdf.groupby(grouper(pdf)).apply(lambda x: x.sum()))

        # SeriesGroupBy
        assert_eq(ddf.groupby(grouper(ddf))['B'].apply(lambda x: x.sum()),
                  pdf.groupby(grouper(pdf))['B'].apply(lambda x: x.sum()))

        # DataFrameGroupBy with column slice
        assert_eq(ddf.groupby(grouper(ddf))[['B', 'C']].apply(lambda x: x.sum()),
                  pdf.groupby(grouper(pdf))[['B', 'C']].apply(lambda x: x.sum()))


def test_numeric_column_names():
    # df.groupby(0)[df.columns] fails if all columns are numbers (pandas bug)
    # This ensures that we cast all column iterables to list beforehand.
    df = pd.DataFrame({0: [0, 1, 0, 1],
                       1: [1, 2, 3, 4],
                       2: [0, 1, 0, 1],})
    ddf = dd.from_pandas(df, npartitions=2)
    assert_eq(ddf.groupby(0).sum(), df.groupby(0).sum())
    assert_eq(ddf.groupby([0, 2]).sum(), df.groupby([0, 2]).sum())
    assert_eq(ddf.groupby(0).apply(lambda x: x, meta={0: int, 1: int, 2: int}),
              df.groupby(0).apply(lambda x: x))


def test_groupby_apply_tasks():
    df = pd.util.testing.makeTimeDataFrame()
    df['A'] = df.A // 0.1
    df['B'] = df.B // 0.1
    ddf = dd.from_pandas(df, npartitions=10)

    with dask.set_options(shuffle='tasks'):
        for ind in [lambda x: 'A', lambda x: x.A]:
            a = df.groupby(ind(df)).apply(len)
            with pytest.warns(UserWarning):
                b = ddf.groupby(ind(ddf)).apply(len)
            assert_eq(a, b.compute())
            assert not any('partd' in k[0] for k in b.dask)

            a = df.groupby(ind(df)).B.apply(len)
            with pytest.warns(UserWarning):
                b = ddf.groupby(ind(ddf)).B.apply(len)
            assert_eq(a, b.compute())
            assert not any('partd' in k[0] for k in b.dask)


def test_groupby_multiprocessing():
    from dask.multiprocessing import get
    df = pd.DataFrame({'A': [1, 2, 3, 4, 5],
                       'B': ['1','1','a','a','a']})
    ddf = dd.from_pandas(df, npartitions=3)
    with dask.set_options(get=get):
        assert_eq(ddf.groupby('B').apply(lambda x: x, meta={"A": int,
                                                            "B": object}),
                  df.groupby('B').apply(lambda x: x))


def test_groupby_normalize_index():
    full = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7, 8, 9],
                         'b': [4, 5, 6, 3, 2, 1, 0, 0, 0]},
                        index=[0, 1, 3, 5, 6, 8, 9, 9, 9])
    d = dd.from_pandas(full, npartitions=3)

    assert d.groupby('a').index == 'a'
    assert d.groupby(d['a']).index == 'a'
    assert d.groupby(d['a'] > 2).index._name == (d['a'] > 2)._name
    assert d.groupby(['a', 'b']).index == ['a', 'b']

    assert d.groupby([d['a'], d['b']]).index == ['a', 'b']
    assert d.groupby([d['a'], 'b']).index == ['a', 'b']


@pytest.mark.parametrize('spec', [
    {'b': {'c': 'mean'}, 'c': {'a': 'max', 'b': 'min'}},
    {'b': 'mean', 'c': ['min', 'max']},
    {'b': np.sum, 'c': ['min', np.max, np.std, np.var]},
    ['sum', 'mean', 'min', 'max', 'count', 'size', 'std', 'var', 'first', 'last'],
    'var',
    {'b':'mean', 'c': 'first', 'd': 'last', 'a': ['first', 'last']},
    {'b': {'c': 'mean'}, 'c':{'a': 'first', 'b': 'last'}},
])
@pytest.mark.parametrize('split_every', [False, None])
@pytest.mark.parametrize('grouper', [
    lambda df: 'a',
    lambda df: ['a', 'd'],
    lambda df: [df['a'], df['d']],
    lambda df: df['a'],
    lambda df: df['a'] > 2,
])
def test_aggregate__examples(spec, split_every, grouper):
    pdf = pd.DataFrame({'a': [1, 2, 3, 1, 1, 2, 4, 3, 7] * 10,
                        'b': [4, 2, 7, 3, 3, 1, 1, 1, 2] * 10,
                        'c': [0, 1, 2, 3, 4, 5, 6, 7, 8] * 10,
                        'd': [3, 2, 1, 3, 2, 1, 2, 6, 4] * 10},
                       columns=['c', 'b', 'a', 'd'])
    ddf = dd.from_pandas(pdf, npartitions=10)

    # Warning from pandas deprecation .agg(dict[dict])
    # it's from pandas, so no reason to assert the deprecation warning,
    # but we should still test it for now
    with pytest.warns(None):
        assert_eq(pdf.groupby(grouper(pdf)).agg(spec),
                  ddf.groupby(grouper(ddf)).agg(spec, split_every=split_every))


@pytest.mark.parametrize('spec', [
    {'b': 'sum', 'c': 'min', 'd': 'max'},
    ['sum'],
    ['sum', 'mean', 'min', 'max', 'count', 'size', 'std', 'var', 'first', 'last'],
    'sum', 'size',
])
@pytest.mark.parametrize('split_every', [False, None])
@pytest.mark.parametrize('grouper', [
    lambda df: [df['a'], df['d']],
    lambda df: df['a'],
    lambda df: df['a'] > 2,
])
def test_series_aggregate__examples(spec, split_every, grouper):
    pdf = pd.DataFrame({'a': [1, 2, 3, 1, 1, 2, 4, 3, 7] * 10,
                        'b': [4, 2, 7, 3, 3, 1, 1, 1, 2] * 10,
                        'c': [0, 1, 2, 3, 4, 5, 6, 7, 8] * 10,
                        'd': [3, 2, 1, 3, 2, 1, 2, 6, 4] * 10},
                       columns=['c', 'b', 'a', 'd'])
    ps = pdf['c']

    ddf = dd.from_pandas(pdf, npartitions=10)
    ds = ddf['c']
    # Warning from pandas deprecation .agg(dict[dict])
    # it's from pandas, so no reason to assert the deprecation warning,
    # but we should still test it for now
    with pytest.warns(None):
        assert_eq(ps.groupby(grouper(pdf)).agg(spec),
                  ds.groupby(grouper(ddf)).agg(spec, split_every=split_every),
                  # pandas < 0.20.0 does not propagate the name for size
                  check_names=(spec != 'size'))


def test_aggregate__single_element_groups(agg_func):
    spec = agg_func

    # nunique is not supported in specs
    if spec == 'nunique':
        return

    pdf = pd.DataFrame({'a': [1, 1, 3, 3],
                        'b': [4, 4, 16, 16],
                        'c': [1, 1, 4, 4],
                        'd': [1, 1, 3, 3]},
                       columns=['c', 'b', 'a', 'd'])
    ddf = dd.from_pandas(pdf, npartitions=3)

    expected = pdf.groupby(['a', 'd']).agg(spec)

    # NOTE: for std the result is not recast ot the original dtype
    if spec in {'mean', 'var'}:
        expected = expected.astype(float)

    assert_eq(expected,
              ddf.groupby(['a', 'd']).agg(spec))


def test_aggregate_build_agg_args__reuse_of_intermediates():
    """Aggregate reuses intermediates. For example, with sum, count, and mean
    the sums and counts are only calculated once accross the graph and reused to
    compute the mean.
    """
    from dask.dataframe.groupby import _build_agg_args

    no_mean_spec = [
        ('foo', 'sum', 'input'),
        ('bar', 'count', 'input'),
    ]

    with_mean_spec = [
        ('foo', 'sum', 'input'),
        ('bar', 'count', 'input'),
        ('baz', 'mean', 'input'),
    ]

    no_mean_chunks, no_mean_aggs, no_mean_finalizers = _build_agg_args(no_mean_spec)
    with_mean_chunks, with_mean_aggs, with_mean_finalizers = _build_agg_args(with_mean_spec)

    assert len(no_mean_chunks) == len(with_mean_chunks)
    assert len(no_mean_aggs) == len(with_mean_aggs)

    assert len(no_mean_finalizers) == len(no_mean_spec)
    assert len(with_mean_finalizers) == len(with_mean_spec)


def test_aggregate__dask():
    dask_holder = collections.namedtuple('dask_holder', ['dask'])
    get_agg_dask = lambda obj: dask_holder({
        k: v for (k, v) in obj.dask.items() if k[0].startswith('aggregate')
    })

    specs = [
        {'b': {'c': 'mean'}, 'c': {'a': 'max', 'b': 'min'}},
        {'b': 'mean', 'c': ['min', 'max']},
        ['sum', 'mean', 'min', 'max', 'count', 'size', 'std', 'var', 'first', 'last'],
        'sum', 'mean', 'min', 'max', 'count', 'std', 'var', 'first', 'last'

        # NOTE: the 'size' spec is special since it bypasses aggregate
        # 'size'
    ]

    pdf = pd.DataFrame({'a': [1, 2, 3, 1, 1, 2, 4, 3, 7] * 100,
                        'b': [4, 2, 7, 3, 3, 1, 1, 1, 2] * 100,
                        'c': [0, 1, 2, 3, 4, 5, 6, 7, 8] * 100,
                        'd': [3, 2, 1, 3, 2, 1, 2, 6, 4] * 100},
                       columns=['c', 'b', 'a', 'd'])
    ddf = dd.from_pandas(pdf, npartitions=100)

    for spec in specs:
        result1 = ddf.groupby(['a', 'b']).agg(spec, split_every=2)
        result2 = ddf.groupby(['a', 'b']).agg(spec, split_every=2)

        agg_dask1 = get_agg_dask(result1)
        agg_dask2 = get_agg_dask(result2)

        # check that the number of partitions used is fixed by split_every
        assert_max_deps(agg_dask1, 2)
        assert_max_deps(agg_dask2, 2)

        # check for deterministic key names and values
        assert agg_dask1 == agg_dask2

        # the length of the dask does not depend on the passed spec
        for other_spec in specs:
            other = ddf.groupby(['a', 'b']).agg(other_spec, split_every=2)
            assert len(other.dask) == len(result1.dask)
            assert len(other.dask) == len(result2.dask)


@pytest.mark.parametrize('grouper', [
    lambda df: ['a'],
    lambda df: ['a', 'b'],
    lambda df: df['a'],
    lambda df: [df['a'], df['b']],
    lambda df: [df['a'] > 2, df['b'] > 1]
])
def test_dataframe_aggregations_multilevel(grouper, agg_func):
    def call(g, m, **kwargs):
        return getattr(g, m)(**kwargs)

    pdf = pd.DataFrame({'a': [1, 2, 6, 4, 4, 6, 4, 3, 7] * 10,
                        'b': [4, 2, 7, 3, 3, 1, 1, 1, 2] * 10,
                        'd': [0, 1, 2, 3, 4, 5, 6, 7, 8] * 10,
                        'c': [0, 1, 2, 3, 4, 5, 6, 7, 8] * 10},
                       columns=['c', 'b', 'a', 'd'])

    ddf = dd.from_pandas(pdf, npartitions=10)

    assert_eq(call(pdf.groupby(grouper(pdf))['c'], agg_func),
              call(ddf.groupby(grouper(ddf))['c'], agg_func, split_every=2))

    # not supported by pandas
    if agg_func != 'nunique':
        assert_eq(call(pdf.groupby(grouper(pdf))[['c', 'd']], agg_func),
                  call(ddf.groupby(grouper(ddf))[['c', 'd']], agg_func, split_every=2))

        assert_eq(call(pdf.groupby(grouper(pdf)), agg_func),
                  call(ddf.groupby(grouper(ddf)), agg_func, split_every=2))


@pytest.mark.parametrize('grouper', [
    lambda df: df['a'],
    lambda df: [df['a'], df['b']],
    lambda df: [df['a'] > 2, df['b'] > 1]
])
def test_series_aggregations_multilevel(grouper, agg_func):
    """
    similar to ``test_dataframe_aggregations_multilevel``, but series do not
    support all groupby args.
    """
    def call(g, m, **kwargs):
        return getattr(g, m)(**kwargs)

    pdf = pd.DataFrame({'a': [1, 2, 6, 4, 4, 6, 4, 3, 7] * 10,
                        'b': [4, 2, 7, 3, 3, 1, 1, 1, 2] * 10,
                        'c': [0, 1, 2, 3, 4, 5, 6, 7, 8] * 10},
                       columns=['c', 'b', 'a'])

    ddf = dd.from_pandas(pdf, npartitions=10)

    assert_eq(call(pdf['c'].groupby(grouper(pdf)), agg_func),
              call(ddf['c'].groupby(grouper(ddf)), agg_func, split_every=2),
              # for pandas ~ 0.18, the name is not not properly propagated for
              # the mean aggregation
              check_names=(agg_func not in {'mean', 'nunique'}))


@pytest.mark.parametrize('grouper', [
    lambda df: df['a'],
    lambda df: df['a'] > 2,
    lambda df: [df['a'], df['b']],
    lambda df: [df['a'] > 2],
    pytest.mark.xfail(reason="index dtype does not coincide: boolean != empty")(lambda df: [df['a'] > 2, df['b'] > 1])
])
@pytest.mark.parametrize('group_and_slice', [
    lambda df, grouper: df.groupby(grouper(df)),
    lambda df, grouper: df['c'].groupby(grouper(df)),
    lambda df, grouper: df.groupby(grouper(df))['c'],
])
def test_groupby_meta_content(group_and_slice, grouper):
    pdf = pd.DataFrame({'a': [1, 2, 6, 4, 4, 6, 4, 3, 7] * 10,
                        'b': [4, 2, 7, 3, 3, 1, 1, 1, 2] * 10,
                        'c': [0, 1, 2, 3, 4, 5, 6, 7, 8] * 10},
                       columns=['c', 'b', 'a'])

    ddf = dd.from_pandas(pdf, npartitions=10)

    expected = group_and_slice(pdf, grouper).first().head(0)
    meta = group_and_slice(ddf, grouper)._meta.first()
    meta_nonempty = group_and_slice(ddf, grouper)._meta_nonempty.first().head(0)

    assert_eq(expected, meta)
    assert_eq(expected, meta_nonempty)


def test_groupy_non_aligned_index():
    pdf = pd.DataFrame({'a': [1, 2, 6, 4, 4, 6, 4, 3, 7] * 10,
                        'b': [4, 2, 7, 3, 3, 1, 1, 1, 2] * 10,
                        'c': [0, 1, 2, 3, 4, 5, 6, 7, 8] * 10},
                       columns=['c', 'b', 'a'])

    ddf3 = dd.from_pandas(pdf, npartitions=3)
    ddf7 = dd.from_pandas(pdf, npartitions=7)

    # working examples
    ddf3.groupby(['a', 'b'])
    ddf3.groupby([ddf3['a'], ddf3['b']])

    # misaligned divisions
    with pytest.raises(NotImplementedError):
        ddf3.groupby(ddf7['a'])

    with pytest.raises(NotImplementedError):
        ddf3.groupby([ddf7['a'], ddf7['b']])

    with pytest.raises(NotImplementedError):
        ddf3.groupby([ddf7['a'], ddf3['b']])

    with pytest.raises(NotImplementedError):
        ddf3.groupby([ddf3['a'], ddf7['b']])

    with pytest.raises(NotImplementedError):
        ddf3.groupby([ddf7['a'], 'b'])


def test_groupy_series_wrong_grouper():
    df = pd.DataFrame({'a': [1, 2, 6, 4, 4, 6, 4, 3, 7] * 10,
                       'b': [4, 2, 7, 3, 3, 1, 1, 1, 2] * 10,
                       'c': [0, 1, 2, 3, 4, 5, 6, 7, 8] * 10},
                      columns=['c', 'b', 'a'])

    df = dd.from_pandas(df, npartitions=3)
    s = df['a']

    # working index values
    s.groupby(s)
    s.groupby([s, s])

    # non working index values
    with pytest.raises(KeyError):
        s.groupby('foo')

    with pytest.raises(KeyError):
        s.groupby([s, 'foo'])

    with pytest.raises(ValueError):
        s.groupby(df)

    with pytest.raises(ValueError):
        s.groupby([s, df])


@pytest.mark.parametrize('npartitions', [1, 4, 20])
@pytest.mark.parametrize('split_every', [2, 5])
@pytest.mark.parametrize('split_out', [None, 1, 5, 20])
def test_hash_groupby_aggregate(npartitions, split_every, split_out):
    df = pd.DataFrame({'x': np.arange(100) % 10,
                       'y': np.ones(100)})
    ddf = dd.from_pandas(df, npartitions)

    result = ddf.groupby('x').y.var(split_every=split_every,
                                    split_out=split_out)

    dsk = result.__dask_optimize__(result.dask, result.__dask_keys__())
    from dask.core import get_deps
    dependencies, dependents = get_deps(dsk)

    assert result.npartitions == (split_out or 1)
    assert len([k for k, v in dependencies.items() if not v]) == npartitions

    assert_eq(result, df.groupby('x').y.var())


def test_split_out_multi_column_groupby():
    df = pd.DataFrame({'x': np.arange(100) % 10,
                       'y': np.ones(100),
                       'z': [1, 2, 3, 4, 5] * 20})

    ddf = dd.from_pandas(df, npartitions=10)

    result = ddf.groupby(['x', 'y']).z.mean(split_out=4)
    expected = df.groupby(['x', 'y']).z.mean()

    assert_eq(result, expected, check_dtype=False)


def test_groupby_split_out_num():
    # GH 1841
    ddf = dd.from_pandas(pd.DataFrame({'A': [1, 1, 2, 2],
                                       'B': [1, 2, 3, 4]}),
                         npartitions=2)
    assert ddf.groupby('A').sum().npartitions == 1
    assert ddf.groupby('A').sum(split_out=2).npartitions == 2
    assert ddf.groupby('A').sum(split_out=3).npartitions == 3

    with pytest.raises(TypeError):
        # groupby doesn't adcept split_out
        ddf.groupby('A', split_out=2)


def test_groupby_not_supported():
    ddf = dd.from_pandas(pd.DataFrame({'A': [1, 1, 2, 2],
                                       'B': [1, 2, 3, 4]}),
                         npartitions=2)
    with pytest.raises(TypeError):
        ddf.groupby('A', axis=1)
    with pytest.raises(TypeError):
        ddf.groupby('A', level=1)
    with pytest.raises(TypeError):
        ddf.groupby('A', as_index=False)
    with pytest.raises(TypeError):
        ddf.groupby('A', sort=False)
    with pytest.raises(TypeError):
        ddf.groupby('A', group_keys=False)
    with pytest.raises(TypeError):
        ddf.groupby('A', squeeze=True)


def test_groupby_numeric_column():
    df = pd.DataFrame({'A' : ['foo', 'foo', 'bar'], 0: [1,2,3]})
    ddf = dd.from_pandas(df, npartitions=3)

    assert_eq(ddf.groupby(ddf.A)[0].sum(),
              df.groupby(df.A)[0].sum())


@pytest.mark.parametrize('sel', ['c', 'd', ['c', 'd']])
@pytest.mark.parametrize('key', ['a', ['a', 'b']])
@pytest.mark.parametrize('func', ['cumsum', 'cumprod', 'cumcount'])
def test_cumulative(func, key, sel):
    df = pd.DataFrame({'a': [1, 2, 6, 4, 4, 6, 4, 3, 7] * 6,
                       'b': [4, 2, 7, 3, 3, 1, 1, 1, 2] * 6,
                       'c': np.random.randn(54),
                       'd': np.random.randn(54)},
                      columns=['a', 'b', 'c', 'd'])
    df.iloc[[-18, -12, -6], -1] = np.nan
    ddf = dd.from_pandas(df, npartitions=10)

    g, dg = [d.groupby(key)[sel] for d in (df, ddf)]
    assert_eq(getattr(g, func)(), getattr(dg, func)())


@pytest.mark.parametrize('func', ['cumsum', 'cumprod'])
def test_cumulative_axis1(func):
    df = pd.DataFrame({'a': [1, 2, 6, 4, 4, 6, 4, 3, 7] * 2,
                       'b': np.random.randn(18),
                       'c': np.random.randn(18)})
    df.iloc[-6, -1] = np.nan
    ddf = dd.from_pandas(df, npartitions=4)
    assert_eq(getattr(df.groupby('a'), func)(axis=1),
              getattr(ddf.groupby('a'), func)(axis=1))


def test_groupby_unaligned_index():
    df = pd.DataFrame({'a': np.random.randint(0, 10, 50),
                       'b': np.random.randn(50),
                       'c': np.random.randn(50)})
    ddf = dd.from_pandas(df, npartitions=5)
    filtered = df[df.b < 0.5]
    dfiltered = ddf[ddf.b < 0.5]

    ddf_group = dfiltered.groupby(ddf.a)
    ds_group = dfiltered.b.groupby(ddf.a)

    bad = [ddf_group.mean(),
           ddf_group.var(),
           ddf_group.b.nunique(),
           ddf_group.get_group(0),
           ds_group.mean(),
           ds_group.var(),
           ds_group.nunique(),
           ds_group.get_group(0)]

    for obj in bad:
        with pytest.raises(ValueError):
            obj.compute()

    def add1(x):
        return x + 1

    df_group = filtered.groupby(df.a)
    good = [(ddf_group.apply(add1, meta=ddf), df_group.apply(add1)),
            (ddf_group.b.apply(add1, meta=ddf.b), df_group.b.apply(add1))]

    for (res, sol) in good:
        assert_eq(res, sol)


def test_groupby_slice_agg_reduces():
    d = pd.DataFrame({"a": [1, 2, 3, 4], "b": [2, 3, 4, 5]})
    a = dd.from_pandas(d, npartitions=2)
    result = a.groupby("a")["b"].agg(['min', 'max'])
    expected = d.groupby("a")['b'].agg(['min', 'max'])
    assert_eq(result, expected)


def test_groupby_agg_grouper_single():
    # https://github.com/dask/dask/issues/2255
    d = pd.DataFrame({'a': [1, 2, 3, 4]})
    a = dd.from_pandas(d, npartitions=2)

    result = a.groupby('a')['a'].agg(['min', 'max'])
    expected = d.groupby('a')['a'].agg(['min', 'max'])
    assert_eq(result, expected)


@pytest.mark.parametrize('slice_', [
    'a', ['a'], ['a', 'b'], ['b'],
])
def test_groupby_agg_grouper_multiple(slice_):
    # https://github.com/dask/dask/issues/2255
    d = pd.DataFrame({'a': [1, 2, 3, 4], 'b': [1, 2, 3, 4]})
    a = dd.from_pandas(d, npartitions=2)

    result = a.groupby('a')[slice_].agg(['min', 'max'])
    expected = d.groupby('a')[slice_].agg(['min', 'max'])
    assert_eq(result, expected)


@pytest.mark.skipif(PANDAS_VERSION < '0.21.0',
                    reason="Need pandas groupby bug fix "
                           "(pandas-dev/pandas#16859)")
@pytest.mark.parametrize('agg_func', [
    'cumprod', 'cumcount', 'cumsum', 'var', 'sum', 'mean', 'count', 'size',
    'std', 'min', 'max', 'first', 'last'
])
def test_groupby_column_and_index_agg_funcs(agg_func):

    def call(g, m, **kwargs):
        return getattr(g, m)(**kwargs)

    df = pd.DataFrame({'idx': [1, 1, 1, 2, 2, 2],
                       'a': [1, 2, 1, 2, 1, 2],
                       'b': np.arange(6),
                       'c': [1, 1, 1, 2, 2, 2]}
                      ).set_index('idx')

    ddf = dd.from_pandas(df, npartitions=df.index.nunique())
    ddf_no_divs = dd.from_pandas(df, npartitions=df.index.nunique(), sort=False)

    # Index and then column

    # Compute expected result
    expected = call(df.groupby(['idx', 'a']), agg_func)
    if agg_func in {'mean', 'var'}:
        expected = expected.astype(float)

    result = call(ddf.groupby(['idx', 'a']), agg_func)
    assert_eq(expected, result)

    result = call(ddf_no_divs.groupby(['idx', 'a']), agg_func)
    assert_eq(expected, result)

    # apply-combine-apply aggregation functions
    aca_agg = {'sum', 'mean', 'var', 'size', 'std', 'count', 'first', 'last'}

    # Test aggregate strings
    if agg_func in aca_agg:
        result = ddf_no_divs.groupby(['idx', 'a']).agg(agg_func)
        assert_eq(expected, result)

    # Column and then index

    # Compute expected result
    expected = call(df.groupby(['a', 'idx']), agg_func)
    if agg_func in {'mean', 'var'}:
        expected = expected.astype(float)

    result = call(ddf.groupby(['a', 'idx']), agg_func)
    assert_eq(expected, result)

    result = call(ddf_no_divs.groupby(['a', 'idx']), agg_func)
    assert_eq(expected, result)

    # Test aggregate strings
    if agg_func in aca_agg:
        result = ddf_no_divs.groupby(['a', 'idx']).agg(agg_func)
        assert_eq(expected, result)

    # Index only

    # Compute expected result
    expected = call(df.groupby('idx'), agg_func)
    if agg_func in {'mean', 'var'}:
        expected = expected.astype(float)

    result = call(ddf.groupby('idx'), agg_func)
    assert_eq(expected, result)

    result = call(ddf_no_divs.groupby('idx'), agg_func)
    assert_eq(expected, result)

    # Test aggregate strings
    if agg_func in aca_agg:
        result = ddf_no_divs.groupby('idx').agg(agg_func)
        assert_eq(expected, result)


@pytest.mark.skipif(PANDAS_VERSION < '0.21.0',
                    reason="Need 0.21.0 for mixed column/index grouping")
@pytest.mark.parametrize(
    'group_args', [['idx', 'a'], ['a', 'idx'], ['idx'], 'idx'])
@pytest.mark.parametrize(
    'apply_func', [np.min, np.mean, lambda s: np.max(s) - np.mean(s)])
def test_groupby_column_and_index_apply(group_args, apply_func):
    df = pd.DataFrame({'idx': [1, 1, 1, 2, 2, 2],
                       'a': [1, 2, 1, 2, 1, 2],
                       'b': np.arange(6)}
                      ).set_index('idx')

    ddf = dd.from_pandas(df, npartitions=df.index.nunique())
    ddf_no_divs = dd.from_pandas(df, npartitions=df.index.nunique(), sort=False)

    # Expected result
    expected = df.groupby(group_args).apply(apply_func)

    # Compute on dask DataFrame with divisions (no shuffling)
    result = ddf.groupby(group_args).apply(apply_func)
    assert_eq(expected, result, check_divisions=False)

    # Check that partitioning is preserved
    assert ddf.divisions == result.divisions

    # Check that no shuffling occurred.
    # The groupby operation should add only 1 task per partition
    assert len(result.dask) == (len(ddf.dask) + ddf.npartitions)

    # Compute on dask DataFrame without divisions (requires shuffling)
    result = ddf_no_divs.groupby(group_args).apply(apply_func)
    assert_eq(expected, result, check_divisions=False)

    # Check that divisions were preserved (all None in this case)
    assert ddf_no_divs.divisions == result.divisions

    # Crude check to see if shuffling was performed.
    # The groupby operation should add only more than 1 task per partition
    assert len(result.dask) > (len(ddf_no_divs.dask) + ddf_no_divs.npartitions)


custom_mean = dd.Aggregation(
    'mean',
    lambda s: (s.count(), s.sum()),
    lambda s0, s1: (s0.sum(), s1.sum()),
    lambda s0, s1: s1 / s0,
)

custom_sum = dd.Aggregation('sum', lambda s: s.sum(), lambda s0: s0.sum())


@pytest.mark.parametrize('pandas_spec, dask_spec, check_dtype', [
    ({'b': 'mean'}, {'b': custom_mean}, False),
    ({'b': 'sum'}, {'b': custom_sum}, True),
    (['mean', 'sum'], [custom_mean, custom_sum], False),
    ({'b': ['mean', 'sum']}, {'b': [custom_mean, custom_sum]}, False),
])
def test_dataframe_groupby_agg_custom_sum(pandas_spec, dask_spec, check_dtype):
    df = pd.DataFrame({'g': [0, 0, 1] * 3, 'b': [1, 2, 3] * 3})
    ddf = dd.from_pandas(df, npartitions=2)

    expected = df.groupby('g').aggregate(pandas_spec)
    result = ddf.groupby('g').aggregate(dask_spec)

    assert_eq(result, expected, check_dtype=check_dtype)


@pytest.mark.parametrize('pandas_spec, dask_spec', [
    ('mean', custom_mean),
    (['mean'], [custom_mean]),
    (['mean', 'sum'], [custom_mean, custom_sum]),
])
def test_series_groupby_agg_custom_mean(pandas_spec, dask_spec):
    d = pd.DataFrame({'g': [0, 0, 1] * 3, 'b': [1, 2, 3] * 3})
    a = dd.from_pandas(d, npartitions=2)

    expected = d['b'].groupby(d['g']).aggregate(pandas_spec)
    result = a['b'].groupby(a['g']).aggregate(dask_spec)

    assert_eq(result, expected, check_dtype=False)


def test_groupby_agg_custom__name_clash_with_internal_same_column():
    """for a single input column only unique names are allowed"""
    d = pd.DataFrame({'g': [0, 0, 1] * 3, 'b': [1, 2, 3] * 3})
    a = dd.from_pandas(d, npartitions=2)

    agg_func = dd.Aggregation('sum', lambda s: s.sum(), lambda s0: s0.sum())

    with pytest.raises(ValueError):
        a.groupby('g').aggregate({'b': [agg_func, 'sum']})


def test_groupby_agg_custom__name_clash_with_internal_different_column():
    """custom aggregation functions can share the name of a builtin function"""
    d = pd.DataFrame({'g': [0, 0, 1] * 3, 'b': [1, 2, 3] * 3, 'c': [4, 5, 6] * 3})
    a = dd.from_pandas(d, npartitions=2)

    # NOTE: this function is purposefully misnamed
    agg_func = dd.Aggregation(
        'sum',
        lambda s: (s.count(), s.sum()),
        lambda s0, s1: (s0.sum(), s1.sum()),
        lambda s0, s1: s1 / s0,
    )

    # NOTE: the name of agg-func is suppressed in the output,
    # since only a single agg func per column was specified
    result = a.groupby('g').aggregate({'b': agg_func, 'c': 'sum'})
    expected = d.groupby('g').aggregate({'b': 'mean', 'c': 'sum'})

    assert_eq(result, expected, check_dtype=False)


def test_groupby_agg_custom__mode():
    # mode function passing intermediates as pure python objects around. to protect
    # results from pandas in apply use return results as single-item lists
    def agg_mode(s):
        def impl(s):
            res, = s.iloc[0]

            for i, in s.iloc[1:]:
                res = res.add(i, fill_value=0)

            return [res]

        return s.apply(impl)

    agg_func = dd.Aggregation(
        'custom_mode',
        lambda s: s.apply(lambda s: [s.value_counts()]),
        agg_mode,
        lambda s: s.map(lambda i: i[0].argmax()),
    )

    d = pd.DataFrame({
        'g0': [0, 0, 0, 1, 1] * 3,
        'g1': [0, 0, 0, 1, 1] * 3,
        'cc': [4, 5, 4, 6, 6] * 3,
    })
    a = dd.from_pandas(d, npartitions=5)

    actual = a['cc'].groupby([a['g0'], a['g1']]).agg(agg_func)

    # cheat to get the correct index
    expected = pd.DataFrame({'g0': [0, 1], 'g1': [0, 1], 'cc': [4, 6]})
    expected = expected['cc'].groupby([expected['g0'], expected['g1']]).agg('sum')

    assert_eq(actual, expected)


def test_groupby_select_column_agg():
    pdf = pd.DataFrame({'A': [1, 2, 3, 1, 2, 3, 1, 2, 4],
                        'B': [-0.776, -0.4, -0.873, 0.054, 1.419, -0.948,
                              -0.967, -1.714, -0.666]})
    ddf = dd.from_pandas(pdf, npartitions=4)
    actual = ddf.groupby('A')['B'].agg('var')
    expected = pdf.groupby('A')['B'].agg('var')
    assert_eq(actual, expected)
