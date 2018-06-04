import numpy as np
import pandas as pd
import pandas.util.testing as tm
import dask.dataframe as dd
from dask.dataframe.utils import (shard_df_on_index, meta_nonempty, make_meta,
                                  raise_on_meta_error, check_meta,
                                  UNKNOWN_CATEGORIES, PANDAS_VERSION)

import pytest


def test_shard_df_on_index():
    df = pd.DataFrame({'x': [1, 2, 3, 4, 5, 6], 'y': list('abdabd')},
                      index=[10, 20, 30, 40, 50, 60])

    result = list(shard_df_on_index(df, [20, 50]))
    assert list(result[0].index) == [10]
    assert list(result[1].index) == [20, 30, 40]
    assert list(result[2].index) == [50, 60]


def test_make_meta():
    df = pd.DataFrame({'a': [1, 2, 3], 'b': list('abc'), 'c': [1., 2., 3.]},
                      index=[10, 20, 30])

    # Pandas dataframe
    meta = make_meta(df)
    assert len(meta) == 0
    assert (meta.dtypes == df.dtypes).all()
    assert isinstance(meta.index, type(df.index))

    # Pandas series
    meta = make_meta(df.a)
    assert len(meta) == 0
    assert meta.dtype == df.a.dtype
    assert isinstance(meta.index, type(df.index))

    # Pandas index
    meta = make_meta(df.index)
    assert isinstance(meta, type(df.index))
    assert len(meta) == 0

    # Dask object
    ddf = dd.from_pandas(df, npartitions=2)
    assert make_meta(ddf) is ddf._meta

    # Dict
    meta = make_meta({'a': 'i8', 'b': 'O', 'c': 'f8'})
    assert isinstance(meta, pd.DataFrame)
    assert len(meta) == 0
    assert (meta.dtypes == df.dtypes).all()
    assert isinstance(meta.index, pd.RangeIndex)

    # Iterable
    meta = make_meta([('a', 'i8'), ('c', 'f8'), ('b', 'O')])
    assert (meta.columns == ['a', 'c', 'b']).all()
    assert len(meta) == 0
    assert (meta.dtypes == df.dtypes[meta.dtypes.index]).all()
    assert isinstance(meta.index, pd.RangeIndex)

    # Tuple
    meta = make_meta(('a', 'i8'))
    assert isinstance(meta, pd.Series)
    assert len(meta) == 0
    assert meta.dtype == 'i8'
    assert meta.name == 'a'

    # With index
    meta = make_meta({'a': 'i8', 'b': 'i4'}, pd.Int64Index([1, 2], name='foo'))
    assert isinstance(meta.index, pd.Int64Index)
    assert len(meta.index) == 0
    meta = make_meta(('a', 'i8'), pd.Int64Index([1, 2], name='foo'))
    assert isinstance(meta.index, pd.Int64Index)
    assert len(meta.index) == 0

    # Categoricals
    meta = make_meta({'a': 'category'})
    assert len(meta.a.cat.categories) == 1
    assert meta.a.cat.categories[0] == UNKNOWN_CATEGORIES
    meta = make_meta(('a', 'category'))
    assert len(meta.cat.categories) == 1
    assert meta.cat.categories[0] == UNKNOWN_CATEGORIES

    # Numpy scalar
    meta = make_meta(np.float64(1.0))
    assert isinstance(meta, np.float64)

    # Python scalar
    meta = make_meta(1.0)
    assert isinstance(meta, np.float64)

    # Timestamp
    x = pd.Timestamp(2000, 1, 1)
    meta = make_meta(x)
    assert meta is x

    # Dtype expressions
    meta = make_meta('i8')
    assert isinstance(meta, np.int64)
    meta = make_meta(float)
    assert isinstance(meta, np.dtype(float).type)
    meta = make_meta(np.dtype('bool'))
    assert isinstance(meta, np.bool_)
    assert pytest.raises(TypeError, lambda: make_meta(None))


def test_meta_nonempty():
    df1 = pd.DataFrame({'A': pd.Categorical(['Alice', 'Bob', 'Carol']),
                        'B': list('abc'),
                        'C': 'bar',
                        'D': np.float32(1),
                        'E': np.int32(1),
                        'F': pd.Timestamp('2016-01-01'),
                        'G': pd.date_range('2016-01-01', periods=3,
                                           tz='America/New_York'),
                        'H': pd.Timedelta('1 hours', 'ms'),
                        'I': np.void(b' '),
                        'J': pd.Categorical([UNKNOWN_CATEGORIES] * 3)},
                       columns=list('DCBAHGFEIJ'))
    df2 = df1.iloc[0:0]
    df3 = meta_nonempty(df2)
    assert (df3.dtypes == df2.dtypes).all()
    assert df3['A'][0] == 'Alice'
    assert df3['B'][0] == 'foo'
    assert df3['C'][0] == 'foo'
    assert df3['D'][0] == np.float32(1)
    assert df3['D'][0].dtype == 'f4'
    assert df3['E'][0] == np.int32(1)
    assert df3['E'][0].dtype == 'i4'
    assert df3['F'][0] == pd.Timestamp('1970-01-01 00:00:00')
    assert df3['G'][0] == pd.Timestamp('1970-01-01 00:00:00',
                                       tz='America/New_York')
    assert df3['H'][0] == pd.Timedelta('1', 'ms')
    assert df3['I'][0] == 'foo'
    assert df3['J'][0] == UNKNOWN_CATEGORIES

    s = meta_nonempty(df2['A'])
    assert s.dtype == df2['A'].dtype
    assert (df3['A'] == s).all()


def test_meta_duplicated():
    df = pd.DataFrame(columns=['A', 'A', 'B'])
    res = meta_nonempty(df)

    exp = pd.DataFrame([['foo', 'foo', 'foo'],
                        ['foo', 'foo', 'foo']],
                       index=['a', 'b'],
                       columns=['A', 'A', 'B'])
    tm.assert_frame_equal(res, exp)


def test_meta_nonempty_empty_categories():
    for dtype in ['O', 'f8', 'M8']:
        # Index
        idx = pd.CategoricalIndex([], pd.Index([], dtype=dtype),
                                  ordered=True, name='foo')
        res = meta_nonempty(idx)
        assert type(res) is pd.CategoricalIndex
        assert type(res.categories) is type(idx.categories)
        assert res.ordered == idx.ordered
        assert res.name == idx.name
        # Series
        s = idx.to_series()
        res = meta_nonempty(s)
        assert res.dtype == 'category'
        assert s.dtype == 'category'
        assert type(res.cat.categories) is type(s.cat.categories)
        assert res.cat.ordered == s.cat.ordered
        assert res.name == s.name


def test_meta_nonempty_index():
    idx = pd.RangeIndex(1, name='foo')
    res = meta_nonempty(idx)
    assert type(res) is pd.RangeIndex
    assert res.name == idx.name

    idx = pd.Int64Index([1], name='foo')
    res = meta_nonempty(idx)
    assert type(res) is pd.Int64Index
    assert res.name == idx.name

    idx = pd.Index(['a'], name='foo')
    res = meta_nonempty(idx)
    assert type(res) is pd.Index
    assert res.name == idx.name

    idx = pd.DatetimeIndex(['1970-01-01'], freq='d',
                           tz='America/New_York', name='foo')
    res = meta_nonempty(idx)
    assert type(res) is pd.DatetimeIndex
    assert res.tz == idx.tz
    assert res.freq == idx.freq
    assert res.name == idx.name

    idx = pd.PeriodIndex(['1970-01-01'], freq='d', name='foo')
    res = meta_nonempty(idx)
    assert type(res) is pd.PeriodIndex
    assert res.freq == idx.freq
    assert res.name == idx.name

    idx = pd.TimedeltaIndex([np.timedelta64(1, 'D')], freq='d', name='foo')
    res = meta_nonempty(idx)
    assert type(res) is pd.TimedeltaIndex
    assert res.freq == idx.freq
    assert res.name == idx.name

    idx = pd.CategoricalIndex(['a'], ['a', 'b'], ordered=True, name='foo')
    res = meta_nonempty(idx)
    assert type(res) is pd.CategoricalIndex
    assert (res.categories == idx.categories).all()
    assert res.ordered == idx.ordered
    assert res.name == idx.name

    idx = pd.CategoricalIndex([], [UNKNOWN_CATEGORIES], ordered=True, name='foo')
    res = meta_nonempty(idx)
    assert type(res) is pd.CategoricalIndex
    assert res.ordered == idx.ordered
    assert res.name == idx.name

    levels = [pd.Int64Index([1], name='a'),
              pd.Float64Index([1.0], name='b')]
    idx = pd.MultiIndex(levels=levels, labels=[[0], [0]], names=['a', 'b'])
    res = meta_nonempty(idx)
    assert type(res) is pd.MultiIndex
    for idx1, idx2 in zip(idx.levels, res.levels):
        assert type(idx1) is type(idx2)
        assert idx1.name == idx2.name
    assert res.names == idx.names

    levels = [pd.Int64Index([1], name='a'),
              pd.CategoricalIndex(data=['b'], categories=['b'], name='b'),
              pd.TimedeltaIndex([np.timedelta64(1, 'D')], name='timedelta')]
    idx = pd.MultiIndex(levels=levels, labels=[[0], [0], [0]], names=['a', 'b', 'timedelta'])
    res = meta_nonempty(idx)
    assert type(res) is pd.MultiIndex
    for idx1, idx2 in zip(idx.levels, res.levels):
        assert type(idx1) is type(idx2)
        assert idx1.name == idx2.name
    assert res.names == idx.names


@pytest.mark.skipif(PANDAS_VERSION < '0.20.0',
                    reason="Pandas < 0.20.0 doesn't support UInt64Index")
def test_meta_nonempty_uint64index():
    idx = pd.UInt64Index([1], name='foo')
    res = meta_nonempty(idx)
    assert type(res) is pd.UInt64Index
    assert res.name == idx.name


def test_meta_nonempty_scalar():
    meta = meta_nonempty(np.float64(1.0))
    assert isinstance(meta, np.float64)

    x = pd.Timestamp(2000, 1, 1)
    meta = meta_nonempty(x)
    assert meta is x


def test_raise_on_meta_error():
    try:
        with raise_on_meta_error():
            raise RuntimeError("Bad stuff")
    except Exception as e:
        assert e.args[0].startswith("Metadata inference failed.\n")
        assert 'RuntimeError' in e.args[0]
    else:
        assert False, "should have errored"

    try:
        with raise_on_meta_error("myfunc"):
            raise RuntimeError("Bad stuff")
    except Exception as e:
        assert e.args[0].startswith("Metadata inference failed in `myfunc`.\n")
        assert 'RuntimeError' in e.args[0]
    else:
        assert False, "should have errored"


def test_check_meta():
    df = pd.DataFrame({'a': ['x', 'y', 'z'],
                       'b': [True, False, True],
                       'c': [1, 2.5, 3.5],
                       'd': [1, 2, 3],
                       'e': pd.Categorical(['x', 'y', 'z'])})
    meta = df.iloc[:0]

    # DataFrame metadata passthrough if correct
    assert check_meta(df, meta) is df
    # Series metadata passthrough if correct
    e = df.e
    assert check_meta(e, meta.e) is e
    # numeric_equal means floats and ints are equivalent
    d = df.d
    assert check_meta(d, meta.d.astype('f8'), numeric_equal=True) is d

    # Series metadata error
    with pytest.raises(ValueError) as err:
        check_meta(d, meta.d.astype('f8'), numeric_equal=False)
    assert str(err.value) == ('Metadata mismatch found.\n'
                              '\n'
                              'Partition type: `Series`\n'
                              '+----------+---------+\n'
                              '|          | dtype   |\n'
                              '+----------+---------+\n'
                              '| Found    | int64   |\n'
                              '| Expected | float64 |\n'
                              '+----------+---------+')

    # DataFrame metadata error
    meta2 = meta.astype({'a': 'category', 'd': 'f8'})[['a', 'b', 'c', 'd']]
    df2 = df[['a', 'b', 'd', 'e']]
    with pytest.raises(ValueError) as err:
        check_meta(df2, meta2, funcname='from_delayed')

    exp = (
        'Metadata mismatch found in `from_delayed`.\n'
        '\n'
        'Partition type: `DataFrame`\n'
        '+--------+----------+----------+\n'
        '| Column | Found    | Expected |\n'
        '+--------+----------+----------+\n'
        '| a      | object   | category |\n'
        '| c      | -        | float64  |\n'
        '| e      | category | -        |\n'
        '+--------+----------+----------+')
    assert str(err.value) == exp
