import os
import pandas as pd
import pytest
import pickle
import numpy as np
import string
from copy import copy
import pandas.util.testing as tm

import dask
import dask.dataframe as dd
from dask import delayed
from dask.base import compute_as_if_collection
from dask.threaded import get as threaded_get
from dask.multiprocessing import get as mp_get
from dask.dataframe.shuffle import (shuffle,
                                    partitioning_index,
                                    rearrange_by_column,
                                    rearrange_by_divisions,
                                    maybe_buffered_partd,
                                    remove_nans)
from dask.dataframe.utils import assert_eq, make_meta


dsk = {('x', 0): pd.DataFrame({'a': [1, 2, 3], 'b': [1, 4, 7]},
                              index=[0, 1, 3]),
       ('x', 1): pd.DataFrame({'a': [4, 5, 6], 'b': [2, 5, 8]},
                              index=[5, 6, 8]),
       ('x', 2): pd.DataFrame({'a': [7, 8, 9], 'b': [3, 6, 9]},
                              index=[9, 9, 9])}
meta = make_meta({'a': 'i8', 'b': 'i8'}, index=pd.Index([], 'i8'))
d = dd.DataFrame(dsk, 'x', meta, [0, 4, 9, 9])
full = d.compute()


shuffle_func = shuffle  # conflicts with keyword argument


@pytest.mark.parametrize('shuffle', ['disk', 'tasks'])
def test_shuffle(shuffle):
    s = shuffle_func(d, d.b, shuffle=shuffle)
    assert isinstance(s, dd.DataFrame)
    assert s.npartitions == d.npartitions

    x = dask.get(s.dask, (s._name, 0))
    y = dask.get(s.dask, (s._name, 1))

    assert not (set(x.b) & set(y.b))  # disjoint
    assert set(s.dask).issuperset(d.dask)

    assert shuffle_func(d, d.b)._name == shuffle_func(d, d.b)._name


def test_default_partitions():
    assert shuffle(d, d.b).npartitions == d.npartitions


def test_shuffle_npartitions_task():
    df = pd.DataFrame({'x': np.random.random(100)})
    ddf = dd.from_pandas(df, npartitions=10)
    s = shuffle(ddf, ddf.x, shuffle='tasks', npartitions=17, max_branch=4)
    sc = s.compute(get=dask.get)
    assert s.npartitions == 17
    assert set(s.dask).issuperset(set(ddf.dask))

    assert len(sc) == len(df)
    assert list(s.columns) == list(df.columns)
    assert (set(map(tuple, sc.values.tolist())) ==
            set(map(tuple, df.values.tolist())))


@pytest.mark.parametrize('method', ['disk', 'tasks'])
def test_index_with_non_series(method):
    from dask.dataframe.tests.test_multi import list_eq
    list_eq(shuffle(d, d.b, shuffle=method),
            shuffle(d, 'b', shuffle=method))


@pytest.mark.parametrize('method', ['disk', 'tasks'])
def test_index_with_dataframe(method):
    res1 = shuffle(d, d[['b']], shuffle=method).compute()
    res2 = shuffle(d, ['b'], shuffle=method).compute()
    res3 = shuffle(d, 'b', shuffle=method).compute()

    assert sorted(res1.values.tolist()) == sorted(res2.values.tolist())
    assert sorted(res1.values.tolist()) == sorted(res3.values.tolist())


@pytest.mark.parametrize('method', ['disk', 'tasks'])
def test_shuffle_from_one_partition_to_one_other(method):
    df = pd.DataFrame({'x': [1, 2, 3]})
    a = dd.from_pandas(df, 1)

    for i in [1, 2]:
        b = shuffle(a, 'x', npartitions=i, shuffle=method)
        assert len(a.compute(get=dask.get)) == len(b.compute(get=dask.get))


@pytest.mark.parametrize('method', ['disk', 'tasks'])
def test_shuffle_empty_partitions(method):
    df = pd.DataFrame({'x': [1, 2, 3] * 10})
    ddf = dd.from_pandas(df, npartitions=3)
    s = shuffle(ddf, ddf.x, npartitions=6, shuffle=method)
    parts = compute_as_if_collection(dd.DataFrame, s.dask, s.__dask_keys__())
    for p in parts:
        assert s.columns == p.columns


df2 = pd.DataFrame({'i32': np.array([1, 2, 3] * 3, dtype='int32'),
                    'f32': np.array([None, 2.5, 3.5] * 3, dtype='float32'),
                    'cat': pd.Series(['a', 'b', 'c'] * 3).astype('category'),
                    'obj': pd.Series(['d', 'e', 'f'] * 3),
                    'bool': np.array([True, False, True] * 3),
                    'dt': pd.Series(pd.date_range('20130101', periods=9)),
                    'dt_tz': pd.Series(pd.date_range('20130101', periods=9, tz='US/Eastern')),
                    'td': pd.Series(pd.timedelta_range('2000', periods=9))})


def test_partitioning_index():
    res = partitioning_index(df2.i32, 3)
    assert ((res < 3) & (res >= 0)).all()
    assert len(np.unique(res)) > 1

    assert (partitioning_index(df2.i32, 3) == partitioning_index(df2.i32, 3)).all()

    res = partitioning_index(df2[['i32']], 3)
    assert ((res < 3) & (res >= 0)).all()
    assert len(np.unique(res)) > 1

    res = partitioning_index(df2[['cat', 'bool', 'f32']], 2)
    assert ((0 <= res) & (res < 2)).all()

    res = partitioning_index(df2.index, 4)
    assert ((res < 4) & (res >= 0)).all()
    assert len(np.unique(res)) > 1


def test_partitioning_index_categorical_on_values():
    df = pd.DataFrame({'a': list(string.ascii_letters),
                       'b': [1, 2, 3, 4] * 13})
    df.a = df.a.astype('category')
    df2 = df.copy()
    df2.a = df2.a.cat.set_categories(list(reversed(df2.a.cat.categories)))

    res = partitioning_index(df.a, 5)
    res2 = partitioning_index(df2.a, 5)
    assert (res == res2).all()

    res = partitioning_index(df, 5)
    res2 = partitioning_index(df2, 5)
    assert (res == res2).all()


@pytest.mark.parametrize('npartitions', [1, 4, 7, pytest.mark.slow(23)])
def test_set_index_tasks(npartitions):
    df = pd.DataFrame({'x': np.random.random(100),
                       'y': np.random.random(100) // 0.2},
                      index=np.random.random(100))

    ddf = dd.from_pandas(df, npartitions=npartitions)

    assert_eq(df.set_index('x'),
              ddf.set_index('x', shuffle='tasks'))

    assert_eq(df.set_index('y'),
              ddf.set_index('y', shuffle='tasks'))

    assert_eq(df.set_index(df.x),
              ddf.set_index(ddf.x, shuffle='tasks'))

    assert_eq(df.set_index(df.x + df.y),
              ddf.set_index(ddf.x + ddf.y, shuffle='tasks'))

    assert_eq(df.set_index(df.x + 1),
              ddf.set_index(ddf.x + 1, shuffle='tasks'))

    assert_eq(df.set_index(df.index),
              ddf.set_index(ddf.index, shuffle='tasks'))


@pytest.mark.parametrize('shuffle', ['disk', 'tasks'])
def test_set_index_self_index(shuffle):
    df = pd.DataFrame({'x': np.random.random(100),
                       'y': np.random.random(100) // 0.2},
                      index=np.random.random(100))

    a = dd.from_pandas(df, npartitions=4)
    b = a.set_index(a.index, shuffle=shuffle)
    assert a is b

    assert_eq(b, df.set_index(df.index))


@pytest.mark.parametrize('shuffle', ['tasks'])
def test_set_index_names(shuffle):
    df = pd.DataFrame({'x': np.random.random(100),
                       'y': np.random.random(100) // 0.2},
                      index=np.random.random(100))

    ddf = dd.from_pandas(df, npartitions=4)

    assert (set(ddf.set_index('x', shuffle=shuffle).dask) ==
            set(ddf.set_index('x', shuffle=shuffle).dask))
    assert (set(ddf.set_index('x', shuffle=shuffle).dask) !=
            set(ddf.set_index('y', shuffle=shuffle).dask))
    assert (set(ddf.set_index('x', max_branch=4, shuffle=shuffle).dask) !=
            set(ddf.set_index('x', max_branch=3, shuffle=shuffle).dask))
    assert (set(ddf.set_index('x', drop=True, shuffle=shuffle).dask) !=
            set(ddf.set_index('x', drop=False, shuffle=shuffle).dask))


@pytest.mark.parametrize('shuffle', ['disk', 'tasks'])
def test_set_index_tasks_2(shuffle):
    df = dd.demo.make_timeseries(
        '2000', '2004', {'value': float, 'name': str, 'id': int},
        freq='2H', partition_freq='1M', seed=1)

    df2 = df.set_index('name', shuffle=shuffle)
    df2.value.sum().compute(get=dask.get)


@pytest.mark.parametrize('shuffle', ['disk', 'tasks'])
def test_set_index_tasks_3(shuffle):
    df = pd.DataFrame(np.random.random((10, 2)), columns=['x', 'y'])
    ddf = dd.from_pandas(df, npartitions=5)

    ddf2 = ddf.set_index('x', shuffle=shuffle, max_branch=2,
                         npartitions=ddf.npartitions)
    df2 = df.set_index('x')
    assert_eq(df2, ddf2)
    assert ddf2.npartitions == ddf.npartitions


@pytest.mark.parametrize('shuffle', ['tasks', 'disk'])
def test_shuffle_sort(shuffle):
    df = pd.DataFrame({'x': [1, 2, 3, 2, 1], 'y': [9, 8, 7, 1, 5]})
    ddf = dd.from_pandas(df, npartitions=3)

    df2 = df.set_index('x').sort_index()
    ddf2 = ddf.set_index('x', shuffle=shuffle)

    assert_eq(ddf2.loc[2:3], df2.loc[2:3])


@pytest.mark.parametrize('shuffle', ['tasks', 'disk'])
@pytest.mark.parametrize('get', [threaded_get, mp_get])
def test_rearrange(shuffle, get):
    df = pd.DataFrame({'x': np.random.random(10)})
    ddf = dd.from_pandas(df, npartitions=4)
    ddf2 = ddf.assign(y=ddf.x % 4)

    result = rearrange_by_column(ddf2, 'y', max_branch=32, shuffle=shuffle)
    assert result.npartitions == ddf.npartitions
    assert set(ddf.dask).issubset(result.dask)

    # Every value in exactly one partition
    a = result.compute(get=get)
    parts = get(result.dask, result.__dask_keys__())
    for i in a.y.drop_duplicates():
        assert sum(i in part.y for part in parts) == 1


def test_rearrange_by_column_with_narrow_divisions():
    from dask.dataframe.tests.test_multi import list_eq
    A = pd.DataFrame({'x': [1, 2, 3, 4, 5, 6], 'y': [1, 1, 2, 2, 3, 4]})
    a = dd.repartition(A, [0, 4, 5])

    df = rearrange_by_divisions(a, 'x', (0, 2, 5))
    list_eq(df, a)


def test_maybe_buffered_partd():
    import partd
    f = maybe_buffered_partd()
    p1 = f()
    assert isinstance(p1.partd, partd.Buffer)
    f2 = pickle.loads(pickle.dumps(f))
    assert not f2.buffer
    p2 = f2()
    assert isinstance(p2.partd, partd.File)


def test_set_index_with_explicit_divisions():
    df = pd.DataFrame({'x': [4, 1, 2, 5]}, index=[10, 20, 30, 40])

    ddf = dd.from_pandas(df, npartitions=2)

    def throw(*args, **kwargs):
        raise Exception()

    with dask.set_options(get=throw):
        ddf2 = ddf.set_index('x', divisions=[1, 3, 5])
    assert ddf2.divisions == (1, 3, 5)

    df2 = df.set_index('x')
    assert_eq(ddf2, df2)

    # Divisions must be sorted
    with pytest.raises(ValueError):
        ddf.set_index('x', divisions=[3, 1, 5])


def test_set_index_divisions_2():
    df = pd.DataFrame({'x': [1, 2, 3, 4, 5, 6], 'y': list('abdabd')})
    ddf = dd.from_pandas(df, 2)

    result = ddf.set_index('y', divisions=['a', 'c', 'd'])
    assert result.divisions == ('a', 'c', 'd')

    assert list(result.compute(get=dask.get).index[-2:]) == ['d', 'd']


def test_set_index_divisions_compute():
    d2 = d.set_index('b', divisions=[0, 2, 9], compute=False)
    d3 = d.set_index('b', divisions=[0, 2, 9], compute=True)

    assert_eq(d2, d3)
    assert_eq(d2, full.set_index('b'))
    assert_eq(d3, full.set_index('b'))
    assert len(d2.dask) > len(d3.dask)

    d4 = d.set_index(d.b, divisions=[0, 2, 9], compute=False)
    d5 = d.set_index(d.b, divisions=[0, 2, 9], compute=True)
    exp = full.copy()
    exp.index = exp.b
    assert_eq(d4, d5)
    assert_eq(d4, exp)
    assert_eq(d5, exp)
    assert len(d4.dask) > len(d5.dask)


def test_set_index_divisions_sorted():
    p1 = pd.DataFrame({'x': [10, 11, 12], 'y': ['a', 'a', 'a']})
    p2 = pd.DataFrame({'x': [13, 14, 15], 'y': ['b', 'b', 'c']})
    p3 = pd.DataFrame({'x': [16, 17, 18], 'y': ['d', 'e', 'e']})

    ddf = dd.DataFrame({('x', 0): p1, ('x', 1): p2, ('x', 2): p3},
                       'x', p1, [None, None, None, None])
    df = ddf.compute()

    def throw(*args, **kwargs):
        raise Exception("Shouldn't have computed")

    with dask.set_options(get=throw):
        res = ddf.set_index('x', divisions=[10, 13, 16, 18], sorted=True)
    assert_eq(res, df.set_index('x'))

    with dask.set_options(get=throw):
        res = ddf.set_index('y', divisions=['a', 'b', 'd', 'e'], sorted=True)
    assert_eq(res, df.set_index('y'))

    # with sorted=True, divisions must be same length as df.divisions
    with pytest.raises(ValueError):
        ddf.set_index('y', divisions=['a', 'b', 'c', 'd', 'e'], sorted=True)

    # Divisions must be sorted
    with pytest.raises(ValueError):
        ddf.set_index('y', divisions=['a', 'b', 'd', 'c'], sorted=True)


@pytest.mark.parametrize('shuffle', ['disk', 'tasks'])
def test_set_index_reduces_partitions_small(shuffle):
    df = pd.DataFrame({'x': np.random.random(100)})
    ddf = dd.from_pandas(df, npartitions=50)

    ddf2 = ddf.set_index('x', shuffle=shuffle, npartitions='auto')
    assert ddf2.npartitions < 10


def make_part(n):
    return pd.DataFrame({'x': np.random.random(n),
                         'y': np.random.random(n)})


@pytest.mark.parametrize('shuffle', ['disk', 'tasks'])
def test_set_index_reduces_partitions_large(shuffle):
    nbytes = 1e6
    nparts = 50
    n = int(nbytes / (nparts * 8))
    ddf = dd.DataFrame({('x', i): (make_part, n) for i in range(nparts)},
                       'x', make_part(1), [None] * (nparts + 1))
    ddf2 = ddf.set_index('x', shuffle=shuffle, npartitions='auto',
                         partition_size=nbytes)
    assert 1 < ddf2.npartitions < 20


@pytest.mark.parametrize('shuffle', ['disk', 'tasks'])
def test_set_index_doesnt_increase_partitions(shuffle):
    nparts = 2
    nbytes = 1e6
    n = int(nbytes / (nparts * 8))
    ddf = dd.DataFrame({('x', i): (make_part, n) for i in range(nparts)},
                       'x', make_part(1), [None] * (nparts + 1))
    ddf2 = ddf.set_index('x', shuffle=shuffle, npartitions='auto',
                         partition_size=nbytes)
    assert ddf2.npartitions <= ddf.npartitions


@pytest.mark.parametrize('shuffle', ['disk', 'tasks'])
def test_set_index_detects_sorted_data(shuffle):
    df = pd.DataFrame({'x': range(100), 'y': range(100)})
    ddf = dd.from_pandas(df, npartitions=10, name='x', sort=False)

    ddf2 = ddf.set_index('x', shuffle=shuffle)
    assert len(ddf2.dask) < ddf.npartitions * 4


def test_set_index_sorts():
    # https://github.com/dask/dask/issues/2288
    vals = np.array([1348550149000000000, 1348550149000000000, 1348558142000000000,
                     1348558142000000000, 1348585928000000000, 1348585928000000000,
                     1348600739000000000, 1348601706000000000, 1348600739000000000,
                     1348601706000000000, 1348614789000000000, 1348614789000000000,
                     1348621037000000000, 1348621038000000000, 1348621040000000000,
                     1348621037000000000, 1348621038000000000, 1348621040000000000,
                     1348637628000000000, 1348638159000000000, 1348638160000000000,
                     1348638159000000000, 1348638160000000000, 1348637628000000000,
                     1348646354000000000, 1348646354000000000, 1348659107000000000,
                     1348657111000000000, 1348659107000000000, 1348657111000000000,
                     1348672876000000000, 1348672876000000000, 1348682787000000000,
                     1348681985000000000, 1348682787000000000, 1348681985000000000,
                     1348728167000000000, 1348728167000000000, 1348730745000000000,
                     1348730745000000000, 1348750198000000000, 1348750198000000000,
                     1348750198000000000, 1348753539000000000, 1348753539000000000,
                     1348753539000000000, 1348754449000000000, 1348754449000000000,
                     1348761333000000000, 1348761554000000000, 1348761610000000000,
                     1348761333000000000, 1348761554000000000, 1348761610000000000,
                     1348782624000000000, 1348782624000000000, 1348782624000000000,
                     1348782624000000000])
    vals = pd.to_datetime(vals, unit='ns')
    breaks = [10, 36, 58]
    dfs = []

    for i in range(len(breaks)):
        lo = sum(breaks[:i])
        hi = sum(breaks[i:i + 1])

        dfs.append(pd.DataFrame({"timestamp": vals[lo:hi]}, index=range(lo, hi)))

    ddf = dd.concat(dfs).clear_divisions()
    assert ddf.set_index("timestamp").index.compute().is_monotonic is True


def test_set_index():
    dsk = {('x', 0): pd.DataFrame({'a': [1, 2, 3], 'b': [4, 2, 6]},
                                  index=[0, 1, 3]),
           ('x', 1): pd.DataFrame({'a': [4, 5, 6], 'b': [3, 5, 8]},
                                  index=[5, 6, 8]),
           ('x', 2): pd.DataFrame({'a': [7, 8, 9], 'b': [9, 1, 8]},
                                  index=[9, 9, 9])}
    d = dd.DataFrame(dsk, 'x', meta, [0, 4, 9, 9])
    full = d.compute()

    d2 = d.set_index('b', npartitions=3)
    assert d2.npartitions == 3
    assert d2.index.name == 'b'
    assert_eq(d2, full.set_index('b'))

    d3 = d.set_index(d.b, npartitions=3)
    assert d3.npartitions == 3
    assert d3.index.name == 'b'
    assert_eq(d3, full.set_index(full.b))

    d4 = d.set_index('b')
    assert d4.index.name == 'b'
    assert_eq(d4, full.set_index('b'))


def test_set_index_interpolate():
    df = pd.DataFrame({'x': [4, 1, 1, 3, 3], 'y': [1., 1, 1, 1, 2]})
    d = dd.from_pandas(df, 2)

    d1 = d.set_index('x', npartitions=3)
    assert d1.npartitions == 3
    assert set(d1.divisions) == set([1, 2, 3, 4])

    d2 = d.set_index('y', npartitions=3)
    assert d2.divisions[0] == 1.
    assert 1. < d2.divisions[1] < d2.divisions[2] < 2.
    assert d2.divisions[3] == 2.


def test_set_index_interpolate_int():
    L = sorted(list(range(0, 200, 10)) * 2)
    df = pd.DataFrame({'x': 2 * L})
    d = dd.from_pandas(df, 2)
    d1 = d.set_index('x', npartitions=10)
    assert all(np.issubdtype(type(x), np.integer) for x in d1.divisions)


def test_set_index_timezone():
    s_naive = pd.Series(pd.date_range('20130101', periods=3))
    s_aware = pd.Series(pd.date_range('20130101', periods=3, tz='US/Eastern'))
    df = pd.DataFrame({'tz': s_aware, 'notz': s_naive})
    d = dd.from_pandas(df, 2)

    d1 = d.set_index('notz', npartitions=2)
    s1 = pd.DatetimeIndex(s_naive.values, dtype=s_naive.dtype)
    assert d1.divisions[0] == s_naive[0] == s1[0]
    assert d1.divisions[-1] == s_naive[2] == s1[2]

    # We currently lose "freq".  Converting data with pandas-defined dtypes
    # to numpy or pure Python can be lossy like this.
    d2 = d.set_index('tz', npartitions=2)
    s2 = pd.DatetimeIndex(s_aware, dtype=s_aware.dtype)
    assert d2.divisions[0] == s2[0]
    assert d2.divisions[-1] == s2[2]
    assert d2.divisions[0].tz == s2[0].tz
    assert d2.divisions[0].tz is not None
    s2badtype = pd.DatetimeIndex(s_aware.values, dtype=s_naive.dtype)
    with pytest.raises(TypeError):
        d2.divisions[0] == s2badtype[0]


@pytest.mark.parametrize('drop', [True, False])
def test_set_index_drop(drop):
    pdf = pd.DataFrame({'A': list('ABAABBABAA'),
                        'B': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                        'C': [1, 2, 3, 2, 1, 3, 2, 4, 2, 3]})
    ddf = dd.from_pandas(pdf, 3)

    assert_eq(ddf.set_index('A', drop=drop),
              pdf.set_index('A', drop=drop))
    assert_eq(ddf.set_index('B', drop=drop),
              pdf.set_index('B', drop=drop))
    assert_eq(ddf.set_index('C', drop=drop),
              pdf.set_index('C', drop=drop))
    assert_eq(ddf.set_index(ddf.A, drop=drop),
              pdf.set_index(pdf.A, drop=drop))
    assert_eq(ddf.set_index(ddf.B, drop=drop),
              pdf.set_index(pdf.B, drop=drop))
    assert_eq(ddf.set_index(ddf.C, drop=drop),
              pdf.set_index(pdf.C, drop=drop))

    # numeric columns
    pdf = pd.DataFrame({0: list('ABAABBABAA'),
                        1: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                        2: [1, 2, 3, 2, 1, 3, 2, 4, 2, 3]})
    ddf = dd.from_pandas(pdf, 3)
    assert_eq(ddf.set_index(0, drop=drop),
              pdf.set_index(0, drop=drop))
    assert_eq(ddf.set_index(2, drop=drop),
              pdf.set_index(2, drop=drop))


def test_set_index_raises_error_on_bad_input():
    df = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6, 7],
                       'b': [7, 6, 5, 4, 3, 2, 1]})
    ddf = dd.from_pandas(df, 2)

    msg = r"Dask dataframe does not yet support multi-indexes"
    with pytest.raises(NotImplementedError) as err:
        ddf.set_index(['a', 'b'])
    assert msg in str(err.value)


def test_set_index_sorted_true():
    df = pd.DataFrame({'x': [1, 2, 3, 4],
                       'y': [10, 20, 30, 40],
                       'z': [4, 3, 2, 1]})
    a = dd.from_pandas(df, 2, sort=False)
    assert not a.known_divisions

    b = a.set_index('x', sorted=True)
    assert b.known_divisions
    assert set(a.dask).issubset(set(b.dask))

    for drop in [True, False]:
        assert_eq(a.set_index('x', drop=drop),
                  df.set_index('x', drop=drop))
        assert_eq(a.set_index(a.x, sorted=True, drop=drop),
                  df.set_index(df.x, drop=drop))
        assert_eq(a.set_index(a.x + 1, sorted=True, drop=drop),
                  df.set_index(df.x + 1, drop=drop))

    with pytest.raises(ValueError):
        a.set_index(a.z, sorted=True)


def test_set_index_sorted_single_partition():
    df = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [1, 0, 1, 0]})
    ddf = dd.from_pandas(df, npartitions=1)
    assert_eq(ddf.set_index('x', sorted=True),
              df.set_index('x'))


def test_set_index_sorted_min_max_same():
    a = pd.DataFrame({'x': [1, 2, 3], 'y': [0, 0, 0]})
    b = pd.DataFrame({'x': [1, 2, 3], 'y': [1, 1, 1]})

    aa = delayed(a)
    bb = delayed(b)

    df = dd.from_delayed([aa, bb], meta=a)
    assert not df.known_divisions

    df2 = df.set_index('y', sorted=True)
    assert df2.divisions == (0, 1, 1)


def test_set_index_empty_partition():
    test_vals = [1, 2, 3]

    converters = [
        int,
        float,
        str,
        lambda x: pd.to_datetime(x, unit='ns'),
    ]

    for conv in converters:
        df = pd.DataFrame([{'x': conv(i), 'y': i} for i in test_vals], columns=['x', 'y'])
        ddf = dd.concat([
            dd.from_pandas(df, npartitions=1),
            dd.from_pandas(df[df.y > df.y.max()], npartitions=1),
        ])

        assert any(ddf.get_partition(p).compute().empty for p in range(ddf.npartitions))
        assert assert_eq(ddf.set_index('x'), df.set_index('x'))


def test_set_index_on_empty():
    test_vals = [1, 2, 3, 4]
    converters = [
        int,
        float,
        str,
        lambda x: pd.to_datetime(x, unit='ns'),
    ]

    for converter in converters:
        df = pd.DataFrame([{'x': converter(x), 'y': x} for x in test_vals])
        ddf = dd.from_pandas(df, npartitions=4)

        assert ddf.npartitions > 1

        ddf = ddf[ddf.y > df.y.max()].set_index('x')
        expected_df = df[df.y > df.y.max()].set_index('x')

        assert assert_eq(ddf, expected_df)
        assert ddf.npartitions == 1


def test_compute_divisions():
    from dask.dataframe.shuffle import compute_divisions
    df = pd.DataFrame({'x': [1, 2, 3, 4],
                       'y': [10, 20, 30, 40],
                       'z': [4, 3, 2, 1]},
                      index=[1, 3, 10, 20])
    a = dd.from_pandas(df, 2, sort=False)
    assert not a.known_divisions

    divisions = compute_divisions(a)
    b = copy(a)
    b.divisions = divisions

    assert_eq(a, b, check_divisions=False)
    assert b.known_divisions


def test_temporary_directory(tmpdir):
    df = pd.DataFrame({'x': np.random.random(100),
                       'y': np.random.random(100),
                       'z': np.random.random(100)})
    ddf = dd.from_pandas(df, npartitions=10, name='x', sort=False)

    with dask.set_options(temporary_directory=str(tmpdir),
                          get=dask.multiprocessing.get):
        ddf2 = ddf.set_index('x', shuffle='disk')
        ddf2.compute()
        assert any(fn.endswith('.partd') for fn in os.listdir(str(tmpdir)))


def test_empty_partitions():
    # See https://github.com/dask/dask/issues/2408
    df = pd.DataFrame({'a': list(range(10))})
    df['b'] = df['a'] % 3
    df['c'] = df['b'].astype(str)

    ddf = dd.from_pandas(df, npartitions=3)
    ddf = ddf.set_index('b')
    ddf = ddf.repartition(npartitions=3)
    ddf.get_partition(0).compute()
    assert_eq(ddf, df.set_index('b'))

    ddf = ddf.set_index('c')
    assert_eq(ddf, df.set_index('b').set_index('c'))


def test_remove_nans():
    tests = [
        ((1, 1, 2), (1, 1, 2)),
        ((None, 1, 2), (1, 1, 2)),
        ((1, None, 2), (1, 2, 2)),
        ((1, 2, None), (1, 2, 2)),
        ((1, 2, None, None), (1, 2, 2, 2)),
        ((None, None, 1, 2), (1, 1, 1, 2)),
        ((1, None, None, 2), (1, 2, 2, 2)),
        ((None, 1, None, 2, None, 3, None), (1, 1, 2, 2, 3, 3, 3)),
    ]

    converters = [
        (int, np.nan),
        (float, np.nan),
        (str, np.nan),
        (lambda x: pd.to_datetime(x, unit='ns'), np.datetime64('NaT')),
    ]

    for conv, none_val in converters:
        for inputs, expected in tests:
            params = [none_val if x is None else conv(x) for x in inputs]
            expected = [conv(x) for x in expected]
            assert remove_nans(params) == expected


@pytest.mark.slow
def test_gh_2730():
    large = pd.DataFrame({'KEY': np.arange(0, 50000)})
    small = pd.DataFrame({'KEY': np.arange(25, 500)})

    dd_left = dd.from_pandas(small, npartitions=3)
    dd_right = dd.from_pandas(large, npartitions=257)

    with dask.set_options(shuffle='tasks', get=dask.get):
        dd_merged = dd_left.merge(dd_right, how='inner', on='KEY')
        result = dd_merged.compute()

    expected = large.merge(small, how='inner', on='KEY')

    tm.assert_frame_equal(
        result.sort_values('KEY').reset_index(drop=True),
        expected
    )
