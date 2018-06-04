from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
from distutils.version import LooseVersion

import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest

import dask
import dask.multiprocessing
import dask.dataframe as dd
from dask.dataframe.utils import assert_eq
from dask.dataframe.io.parquet import _parse_pandas_metadata
from dask.utils import natural_sort_key

try:
    import fastparquet
except ImportError:
    fastparquet = False

try:
    import pyarrow.parquet as pq
except ImportError:
    pq = False


try:
    import pyarrow as pa
    check_pa_divs = pa.__version__ >= LooseVersion('0.9.0')
except ImportError:
    check_pa_divs = False


def should_check_divs(engine):
    if engine == 'fastparquet':
        return True
    elif engine == 'pyarrow' and check_pa_divs:
        return True
    return False


nrows = 40
npartitions = 15
df = pd.DataFrame({'x': [i * 7 % 5 for i in range(nrows)],  # Not sorted
                   'y': [i * 2.5 for i in range(nrows)]  # Sorted
                   },
                  index=pd.Index([10 * i for i in range(nrows)], name='myindex'))

ddf = dd.from_pandas(df, npartitions=npartitions)


@pytest.fixture(params=[pytest.mark.skipif(not fastparquet, 'fastparquet',
                                           reason='fastparquet not found'),
                        pytest.mark.skipif(not pq, 'pyarrow',
                                           reason='pyarrow not found')])
def engine(request):
    return request.param


def check_fastparquet():
    if not fastparquet:
        pytest.skip('fastparquet not found')


def check_pyarrow():
    if not pq:
        pytest.skip('pyarrow not found')


def write_read_engines(**kwargs):
    """Product of both engines for write/read:

    To add custom marks, pass keyword of the form: `mark_writer_reader=reason`,
    or `mark_engine=reason` to apply to all parameters with that engine."""
    backends = {'pyarrow', 'fastparquet'}
    marks = {(w, r): [] for w in backends for r in backends}

    # Skip if uninstalled
    for name, exists in [('fastparquet', fastparquet), ('pyarrow', pq)]:
        val = pytest.mark.skip(reason='%s not found' % name)
        if not exists:
            for k in marks:
                if name in k:
                    marks[k].append(val)

    # Custom marks
    for kw, val in kwargs.items():
        kind, rest = kw.split('_', 1)
        key = tuple(rest.split('_'))
        if (kind not in ('xfail', 'skip') or len(key) > 2 or
                set(key).difference(backends)):
            raise ValueError("unknown keyword %r" % kw)
        val = getattr(pytest.mark, kind)(reason=val)
        if len(key) == 2:
            marks[key].append(val)
        else:
            for k in marks:
                if key in k:
                    marks[k].append(val)

    return pytest.mark.parametrize(('write_engine', 'read_engine'),
                                   [pytest.param(*k, marks=tuple(v))
                                    for (k, v) in sorted(marks.items())])


pyarrow_fastparquet_msg = "fastparquet fails reading pyarrow written directories"
write_read_engines_xfail = write_read_engines(xfail_pyarrow_fastparquet=pyarrow_fastparquet_msg)


@write_read_engines_xfail
def test_local(tmpdir, write_engine, read_engine):
    tmp = str(tmpdir)
    data = pd.DataFrame({'i32': np.arange(1000, dtype=np.int32),
                         'i64': np.arange(1000, dtype=np.int64),
                         'f': np.arange(1000, dtype=np.float64),
                         'bhello': np.random.choice(['hello', 'yo', 'people'], size=1000).astype("O")})
    df = dd.from_pandas(data, chunksize=500)

    df.to_parquet(tmp, write_index=False, engine=write_engine)

    files = os.listdir(tmp)
    assert '_common_metadata' in files
    assert 'part.0.parquet' in files

    df2 = dd.read_parquet(tmp, index=False, engine=read_engine)

    assert len(df2.divisions) > 1

    out = df2.compute(get=dask.get).reset_index()

    for column in df.columns:
        assert (data[column] == out[column]).all()


@write_read_engines_xfail
def test_index(tmpdir, write_engine, read_engine):
    fn = str(tmpdir)
    ddf.to_parquet(fn, engine=write_engine)

    # Infer divisions for engines/versions that support it
    ddf2 = dd.read_parquet(fn, engine=read_engine, infer_divisions=should_check_divs(read_engine))
    assert_eq(ddf, ddf2, check_divisions=should_check_divs(read_engine))

    # infer_divisions False
    ddf2_no_divs = dd.read_parquet(fn, engine=read_engine, infer_divisions=False)
    assert_eq(ddf.clear_divisions(), ddf2_no_divs, check_divisions=True)

    # infer_divisions unspecified
    ddf2_default = dd.read_parquet(fn, engine=read_engine)

    if read_engine == 'fastparquet':
        # The fastparquet engine infers divisions by default because it only supports reading datasets that have a
        # global _metadata file
        assert_eq(ddf, ddf2_default, check_divisions=True)
    else:
        # pyarrow does not infer divisions by default because doing so requires reading metadata from each file in
        # the dataset, which could be expensive
        assert_eq(ddf.clear_divisions(), ddf2_default, check_divisions=True)


@pytest.mark.parametrize('index', [False, True])
@write_read_engines_xfail
def test_empty(tmpdir, write_engine, read_engine, index):
    fn = str(tmpdir)
    df = pd.DataFrame({'a': ['a', 'b', 'b'], 'b': [4, 5, 6]})[:0]
    if index:
        df.set_index('a', inplace=True, drop=True)
    ddf = dd.from_pandas(df, npartitions=2)

    ddf.to_parquet(fn, write_index=index, engine=write_engine)
    read_df = dd.read_parquet(fn, engine=read_engine)
    assert_eq(ddf, read_df)


@write_read_engines()
def test_read_glob(tmpdir, write_engine, read_engine):
    fn = str(tmpdir)
    ddf.to_parquet(fn, engine=write_engine)
    if os.path.exists(os.path.join(fn, '_metadata')):
        os.unlink(os.path.join(fn, '_metadata'))

    files = os.listdir(fn)
    assert '_metadata' not in files

    # Infer divisions for engines/versions that support it

    ddf2 = dd.read_parquet(os.path.join(fn, '*'), engine=read_engine,
                           infer_divisions=should_check_divs(write_engine) and should_check_divs(read_engine))
    assert_eq(ddf, ddf2, check_divisions=should_check_divs(write_engine) and should_check_divs(read_engine))

    # No divisions
    ddf2_no_divs = dd.read_parquet(os.path.join(fn, '*'), engine=read_engine, infer_divisions=False)
    assert_eq(ddf.clear_divisions(), ddf2_no_divs, check_divisions=True)


@write_read_engines()
def test_read_list(tmpdir, write_engine, read_engine):
    tmpdir = str(tmpdir)
    ddf.to_parquet(tmpdir, engine=write_engine)
    files = sorted([os.path.join(tmpdir, f)
                   for f in os.listdir(tmpdir)
                   if not f.endswith('_metadata')],
                   key=natural_sort_key)

    # Infer divisions for engines/versions that support it
    ddf2 = dd.read_parquet(files, engine=read_engine,
                           infer_divisions=should_check_divs(write_engine) and should_check_divs(read_engine))
    assert_eq(ddf, ddf2, check_divisions=should_check_divs(write_engine) and should_check_divs(read_engine))

    # No divisions
    ddf2_no_divs = dd.read_parquet(files, engine=read_engine, infer_divisions=False)
    assert_eq(ddf.clear_divisions(), ddf2_no_divs, check_divisions=True)


@write_read_engines_xfail
def test_columns_index(tmpdir, write_engine, read_engine):
    fn = str(tmpdir)
    ddf.to_parquet(fn, engine=write_engine)

    # With Index
    # ----------
    # ### Emtpy columns ###
    # With divisions if supported
    assert_eq(dd.read_parquet(fn, columns=[], engine=read_engine, infer_divisions=should_check_divs(read_engine)),
              ddf[[]], check_divisions=should_check_divs(read_engine))

    # No divisions
    assert_eq(dd.read_parquet(fn, columns=[], engine=read_engine, infer_divisions=False),
              ddf[[]].clear_divisions(), check_divisions=True)

    # ### Single column, auto select index ###
    # With divisions if supported
    assert_eq(dd.read_parquet(fn, columns=['x'], engine=read_engine, infer_divisions=should_check_divs(read_engine)),
              ddf[['x']], check_divisions=should_check_divs(read_engine))

    # No divisions
    assert_eq(dd.read_parquet(fn, columns=['x'], engine=read_engine, infer_divisions=False),
              ddf[['x']].clear_divisions(), check_divisions=True)

    # ### Single column, specify index ###
    # With divisions if supported
    assert_eq(dd.read_parquet(fn, index='myindex', columns=['x'], engine=read_engine,
                              infer_divisions=should_check_divs(read_engine)),
              ddf[['x']], check_divisions=should_check_divs(read_engine))

    # No divisions
    assert_eq(dd.read_parquet(fn, index='myindex', columns=['x'], engine=read_engine,
                              infer_divisions=False),
              ddf[['x']].clear_divisions(), check_divisions=True)

    # ### Two columns, specify index ###
    # With divisions if supported
    assert_eq(dd.read_parquet(fn, index='myindex', columns=['x', 'y'], engine=read_engine,
                              infer_divisions=should_check_divs(read_engine)),
              ddf, check_divisions=should_check_divs(read_engine))

    # No divisions
    assert_eq(dd.read_parquet(fn, index='myindex', columns=['x', 'y'], engine=read_engine,
                              infer_divisions=False),
              ddf.clear_divisions(), check_divisions=True)


@write_read_engines_xfail
def test_columns_no_index(tmpdir, write_engine, read_engine):
    fn = str(tmpdir)
    ddf.to_parquet(fn, engine=write_engine)
    ddf2 = ddf.reset_index()

    # No Index
    # --------
    # All columns, none as index
    assert_eq(dd.read_parquet(fn, index=False, engine=read_engine, infer_divisions=False),
              ddf2, check_index=False, check_divisions=True)

    # Two columns, none as index
    assert_eq(dd.read_parquet(fn, index=False, columns=['x', 'y'], engine=read_engine,
                              infer_divisions=False),
              ddf2[['x', 'y']], check_index=False, check_divisions=True)

    # One column and one index, all as columns
    assert_eq(dd.read_parquet(fn, index=False, columns=['myindex', 'x'], engine=read_engine,
                              infer_divisions=False),
              ddf2[['myindex', 'x']], check_index=False, check_divisions=True)


@write_read_engines_xfail
def test_infer_divisions_not_sorted(tmpdir, write_engine, read_engine):
    fn = str(tmpdir)
    ddf.to_parquet(fn, engine=write_engine)

    if read_engine == 'pyarrow' and not check_pa_divs:
        match = 'requires pyarrow >=0.9.0'
        ex = NotImplementedError
    else:
        match = 'not known to be sorted across partitions'
        ex = ValueError

    with pytest.raises(ex, match=match):
        dd.read_parquet(fn, index='x', engine=read_engine, infer_divisions=True)


@write_read_engines_xfail
def test_infer_divisions_no_index(tmpdir, write_engine, read_engine):
    fn = str(tmpdir)
    ddf.to_parquet(fn, engine=write_engine, write_index=False)

    if read_engine == 'pyarrow' and not check_pa_divs:
        match = 'requires pyarrow >=0.9.0'
        ex = NotImplementedError
    else:
        match = 'no index column was discovered'
        ex = ValueError

    with pytest.raises(ex, match=match):
        dd.read_parquet(fn, engine=read_engine, infer_divisions=True)


def test_columns_index_with_multi_index(tmpdir, engine):
    fn = os.path.join(str(tmpdir), 'test.parquet')
    index = pd.MultiIndex.from_arrays([np.arange(10), np.arange(10) + 1],
                                      names=['x0', 'x1'])
    df = pd.DataFrame(np.random.randn(10, 2), columns=['a', 'b'], index=index)
    df2 = df.reset_index(drop=False)

    if engine == 'fastparquet':
        fastparquet.write(fn, df, write_index=True)

        # fastparquet doesn't support multi-index
        with pytest.raises(ValueError):
            ddf = dd.read_parquet(fn, engine=engine)
    else:
        import pyarrow as pa
        pq.write_table(pa.Table.from_pandas(df), fn)

        # Pyarrow supports multi-index reads
        ddf = dd.read_parquet(fn, engine=engine)
        assert_eq(ddf, df)

        d = dd.read_parquet(fn, columns='a', engine=engine)
        assert_eq(d, df['a'])

        d = dd.read_parquet(fn, index=['a', 'b'], columns=['x0', 'x1'], engine=engine)
        assert_eq(d, df2.set_index(['a', 'b'])[['x0', 'x1']])

    # Just index
    d = dd.read_parquet(fn, index=False, engine=engine)
    assert_eq(d, df2)

    d = dd.read_parquet(fn, index=['a'], engine=engine)
    assert_eq(d, df2.set_index('a')[['b']])

    d = dd.read_parquet(fn, index=['x0'], engine=engine)
    assert_eq(d, df2.set_index('x0')[['a', 'b']])

    # Just columns
    d = dd.read_parquet(fn, columns=['x0', 'a'], engine=engine)
    assert_eq(d, df2.set_index('x1')[['x0', 'a']])

    # Both index and columns
    d = dd.read_parquet(fn, index=False, columns=['x0', 'b'], engine=engine)
    assert_eq(d, df2[['x0', 'b']])

    for index in ['x1', 'b']:
        d = dd.read_parquet(fn, index=index, columns=['x0', 'a'], engine=engine)
        assert_eq(d, df2.set_index(index)[['x0', 'a']])

    # Columns and index intersect
    for index in ['a', 'x0']:
        with pytest.raises(ValueError):
            d = dd.read_parquet(fn, index=index, columns=['x0', 'a'], engine=engine)

    # Series output
    for ind, col, sol_df in [(None, 'x0', df2.set_index('x1')),
                             (False, 'b', df2),
                             (False, 'x0', df2),
                             ('a', 'x0', df2.set_index('a')),
                             ('a', 'b', df2.set_index('a'))]:
        d = dd.read_parquet(fn, index=ind, columns=col, engine=engine)
        assert_eq(d, sol_df[col])


@write_read_engines_xfail
def test_no_index(tmpdir, write_engine, read_engine):
    fn = str(tmpdir)
    df = pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]})
    ddf = dd.from_pandas(df, npartitions=2)
    ddf.to_parquet(fn, write_index=False, engine=write_engine)
    ddf2 = dd.read_parquet(fn, engine=read_engine)
    assert_eq(df, ddf2, check_index=False)


def test_read_series(tmpdir, engine):
    fn = str(tmpdir)
    ddf.to_parquet(fn, engine=engine)
    ddf2 = dd.read_parquet(fn, columns=['x'], engine=engine, infer_divisions=should_check_divs(engine))
    assert_eq(ddf[['x']], ddf2, check_divisions=should_check_divs(engine))

    ddf2 = dd.read_parquet(fn, columns='x', index='myindex', engine=engine, infer_divisions=should_check_divs(engine))
    assert_eq(ddf.x, ddf2, check_divisions=should_check_divs(engine))


def test_names(tmpdir, engine):
    fn = str(tmpdir)
    ddf.to_parquet(fn, engine=engine)

    def read(fn, **kwargs):
        return dd.read_parquet(fn, engine=engine, **kwargs)

    assert (set(read(fn).dask) == set(read(fn).dask))

    assert (set(read(fn).dask) !=
            set(read(fn, columns=['x']).dask))

    assert (set(read(fn, columns=('x',)).dask) ==
            set(read(fn, columns=['x']).dask))


@pytest.mark.parametrize('c', [['x'], 'x', ['x', 'y'], []])
def test_optimize(tmpdir, c):
    check_fastparquet()
    fn = str(tmpdir)
    ddf.to_parquet(fn)
    ddf2 = dd.read_parquet(fn)
    assert_eq(df[c], ddf2[c])
    x = ddf2[c]

    dsk = x.__dask_optimize__(x.dask, x.__dask_keys__())
    assert len(dsk) == x.npartitions
    assert all(v[4] == c for v in dsk.values())


@pytest.mark.skipif(not hasattr(pd.DataFrame, 'to_parquet'),
                    reason="no to_parquet method")
@write_read_engines()
def test_roundtrip_from_pandas(tmpdir, write_engine, read_engine):
    fn = str(tmpdir.join('test.parquet'))
    df = pd.DataFrame({'x': [1, 2, 3]})
    df.to_parquet(fn, engine=write_engine)
    ddf = dd.read_parquet(fn, engine=read_engine)
    assert_eq(df, ddf)


@write_read_engines_xfail
def test_categorical(tmpdir, write_engine, read_engine):
    tmp = str(tmpdir)
    df = pd.DataFrame({'x': ['a', 'b', 'c'] * 100}, dtype='category')
    ddf = dd.from_pandas(df, npartitions=3)
    dd.to_parquet(ddf, tmp, engine=write_engine)

    ddf2 = dd.read_parquet(tmp, categories='x', engine=read_engine)
    assert ddf2.compute().x.cat.categories.tolist() == ['a', 'b', 'c']

    ddf2 = dd.read_parquet(tmp, categories=['x'], engine=read_engine)
    assert ddf2.compute().x.cat.categories.tolist() == ['a', 'b', 'c']

    # autocat
    if read_engine != 'pyarrow':
        ddf2 = dd.read_parquet(tmp, engine=read_engine)
        assert ddf2.compute().x.cat.categories.tolist() == ['a', 'b', 'c']

        ddf2.loc[:1000].compute()
        df.index.name = 'index'  # defaults to 'index' in this case
        assert assert_eq(df, ddf2)

    # dereference cats
    ddf2 = dd.read_parquet(tmp, categories=[], engine=read_engine)

    ddf2.loc[:1000].compute()
    assert (df.x == ddf2.x).all()


def test_append(tmpdir, engine):
    """Test that appended parquet equal to the original one."""
    check_fastparquet()
    tmp = str(tmpdir)
    df = pd.DataFrame({'i32': np.arange(1000, dtype=np.int32),
                       'i64': np.arange(1000, dtype=np.int64),
                       'f': np.arange(1000, dtype=np.float64),
                       'bhello': np.random.choice(['hello', 'yo', 'people'],
                                                  size=1000).astype("O")})
    df.index.name = 'index'

    half = len(df) // 2
    ddf1 = dd.from_pandas(df.iloc[:half], chunksize=100)
    ddf2 = dd.from_pandas(df.iloc[half:], chunksize=100)
    ddf1.to_parquet(tmp)
    ddf2.to_parquet(tmp, append=True)

    ddf3 = dd.read_parquet(tmp, engine=engine)
    assert_eq(df, ddf3)


def test_append_create(tmpdir):
    """Test that appended parquet equal to the original one."""
    check_fastparquet()
    tmp = str(tmpdir)
    df = pd.DataFrame({'i32': np.arange(1000, dtype=np.int32),
                       'i64': np.arange(1000, dtype=np.int64),
                       'f': np.arange(1000, dtype=np.float64),
                       'bhello': np.random.choice(['hello', 'yo', 'people'],
                                                  size=1000).astype("O")})
    df.index.name = 'index'

    half = len(df) // 2
    ddf1 = dd.from_pandas(df.iloc[:half], chunksize=100)
    ddf2 = dd.from_pandas(df.iloc[half:], chunksize=100)
    ddf1.to_parquet(tmp, append=True)
    ddf2.to_parquet(tmp, append=True)

    ddf3 = dd.read_parquet(tmp, engine='fastparquet')
    assert_eq(df, ddf3)


def test_append_with_partition(tmpdir):
    check_fastparquet()
    tmp = str(tmpdir)
    df0 = pd.DataFrame({'lat': np.arange(0, 10), 'lon': np.arange(10, 20),
                        'value': np.arange(100, 110)})
    df0.index.name = 'index'
    df1 = pd.DataFrame({'lat': np.arange(10, 20), 'lon': np.arange(10, 20),
                        'value': np.arange(120, 130)})
    df1.index.name = 'index'
    dd_df0 = dd.from_pandas(df0, npartitions=1)
    dd_df1 = dd.from_pandas(df1, npartitions=1)
    dd.to_parquet(dd_df0, tmp, partition_on=['lon'])
    dd.to_parquet(dd_df1, tmp, partition_on=['lon'], append=True,
                  ignore_divisions=True)

    out = dd.read_parquet(tmp).compute()
    out['lon'] = out.lon.astype('int64')  # just to pass assert
    # sort required since partitioning breaks index order
    assert_eq(out.sort_values('value'), pd.concat([df0, df1])[out.columns],
              check_index=False)


def test_partition_on_cats(tmpdir):
    check_fastparquet()
    tmp = str(tmpdir)
    d = pd.DataFrame({'a': np.random.rand(50),
                      'b': np.random.choice(['x', 'y', 'z'], size=50),
                      'c': np.random.choice(['x', 'y', 'z'], size=50)})
    d = dd.from_pandas(d, 2)
    d.to_parquet(tmp, partition_on=['b'], engine='fastparquet')
    df = dd.read_parquet(tmp, engine='fastparquet')
    assert set(df.b.cat.categories) == {'x', 'y', 'z'}

    d.to_parquet(tmp, partition_on=['b', 'c'], engine='fastparquet')
    df = dd.read_parquet(tmp, engine='fastparquet')
    assert set(df.b.cat.categories) == {'x', 'y', 'z'}
    assert set(df.c.cat.categories) == {'x', 'y', 'z'}
    df = dd.read_parquet(tmp, columns=['a', 'c'], engine='fastparquet')
    assert set(df.c.cat.categories) == {'x', 'y', 'z'}
    assert 'b' not in df.columns
    df = dd.read_parquet(tmp, index='c', engine='fastparquet')
    assert set(df.index.categories) == {'x', 'y', 'z'}
    assert 'c' not in df.columns
    # series
    df = dd.read_parquet(tmp, columns='b', engine='fastparquet')
    assert set(df.cat.categories) == {'x', 'y', 'z'}


def test_append_wo_index(tmpdir):
    """Test append with write_index=False."""
    check_fastparquet()
    tmp = str(tmpdir.join('tmp1.parquet'))
    df = pd.DataFrame({'i32': np.arange(1000, dtype=np.int32),
                       'i64': np.arange(1000, dtype=np.int64),
                       'f': np.arange(1000, dtype=np.float64),
                       'bhello': np.random.choice(['hello', 'yo', 'people'],
                                                  size=1000).astype("O")})
    half = len(df) // 2
    ddf1 = dd.from_pandas(df.iloc[:half], chunksize=100)
    ddf2 = dd.from_pandas(df.iloc[half:], chunksize=100)
    ddf1.to_parquet(tmp)

    with pytest.raises(ValueError) as excinfo:
        ddf2.to_parquet(tmp, write_index=False, append=True)
    assert 'Appended columns' in str(excinfo.value)

    tmp = str(tmpdir.join('tmp2.parquet'))
    ddf1.to_parquet(tmp, write_index=False)
    ddf2.to_parquet(tmp, write_index=False, append=True)

    ddf3 = dd.read_parquet(tmp, index='f')
    assert_eq(df.set_index('f'), ddf3)


def test_append_overlapping_divisions(tmpdir):
    """Test raising of error when divisions overlapping."""
    check_fastparquet()
    tmp = str(tmpdir)
    df = pd.DataFrame({'i32': np.arange(1000, dtype=np.int32),
                       'i64': np.arange(1000, dtype=np.int64),
                       'f': np.arange(1000, dtype=np.float64),
                       'bhello': np.random.choice(['hello', 'yo', 'people'],
                                                  size=1000).astype("O")})
    half = len(df) // 2
    ddf1 = dd.from_pandas(df.iloc[:half], chunksize=100)
    ddf2 = dd.from_pandas(df.iloc[half - 10:], chunksize=100)
    ddf1.to_parquet(tmp)

    with pytest.raises(ValueError) as excinfo:
        ddf2.to_parquet(tmp, append=True)
    assert 'Appended divisions' in str(excinfo.value)

    ddf2.to_parquet(tmp, append=True, ignore_divisions=True)


def test_append_different_columns(tmpdir):
    """Test raising of error when non equal columns."""
    check_fastparquet()
    tmp = str(tmpdir)
    df1 = pd.DataFrame({'i32': np.arange(100, dtype=np.int32)})
    df2 = pd.DataFrame({'i64': np.arange(100, dtype=np.int64)})
    df3 = pd.DataFrame({'i32': np.arange(100, dtype=np.int64)})

    ddf1 = dd.from_pandas(df1, chunksize=2)
    ddf2 = dd.from_pandas(df2, chunksize=2)
    ddf3 = dd.from_pandas(df3, chunksize=2)

    ddf1.to_parquet(tmp)

    with pytest.raises(ValueError) as excinfo:
        ddf2.to_parquet(tmp, append=True)
    assert 'Appended columns' in str(excinfo.value)

    with pytest.raises(ValueError) as excinfo:
        ddf3.to_parquet(tmp, append=True)
    assert 'Appended dtypes' in str(excinfo.value)


@write_read_engines_xfail
def test_ordering(tmpdir, write_engine, read_engine):
    tmp = str(tmpdir)
    df = pd.DataFrame({'a': [1, 2, 3],
                       'b': [10, 20, 30],
                       'c': [100, 200, 300]},
                      index=pd.Index([-1, -2, -3], name='myindex'),
                      columns=['c', 'a', 'b'])
    ddf = dd.from_pandas(df, npartitions=2)
    dd.to_parquet(ddf, tmp, engine=write_engine)

    if read_engine == 'fastparquet':
        pf = fastparquet.ParquetFile(tmp)
        assert pf.columns == ['myindex', 'c', 'a', 'b']

    ddf2 = dd.read_parquet(tmp, index='myindex', engine=read_engine)
    assert_eq(ddf, ddf2, check_divisions=False)


def test_read_parquet_custom_columns(tmpdir, engine):
    tmp = str(tmpdir)
    data = pd.DataFrame({'i32': np.arange(1000, dtype=np.int32),
                         'f': np.arange(1000, dtype=np.float64)})
    df = dd.from_pandas(data, chunksize=50)
    df.to_parquet(tmp)

    df2 = dd.read_parquet(tmp,
                          columns=['i32', 'f'],
                          engine=engine,
                          infer_divisions=should_check_divs(engine))
    assert_eq(df[['i32', 'f']], df2,
              check_index=False, check_divisions=should_check_divs(engine))

    df3 = dd.read_parquet(tmp,
                          columns=['f', 'i32'],
                          engine=engine,
                          infer_divisions=should_check_divs(engine))
    assert_eq(df[['f', 'i32']], df3,
              check_index=False, check_divisions=should_check_divs(engine))


@pytest.mark.parametrize('df,write_kwargs,read_kwargs', [
    (pd.DataFrame({'x': [3, 2, 1]}), {}, {}),
    (pd.DataFrame({'x': ['c', 'a', 'b']}), {'object_encoding': 'utf8'}, {}),
    (pd.DataFrame({'x': ['cc', 'a', 'bbb']}), {'object_encoding': 'utf8'}, {}),
    (pd.DataFrame({'x': [b'a', b'b', b'c']}), {'object_encoding': 'bytes'}, {}),
    (pd.DataFrame({'x': pd.Categorical(['a', 'b', 'a'])}),
     {'object_encoding': 'utf8'}, {'categories': ['x']}),
    (pd.DataFrame({'x': pd.Categorical([1, 2, 1])}), {}, {'categories': ['x']}),
    (pd.DataFrame({'x': list(map(pd.Timestamp, [3000, 2000, 1000]))}), {}, {}),
    (pd.DataFrame({'x': [3000, 2000, 1000]}).astype('M8[ns]'), {}, {}),
    pytest.mark.xfail((pd.DataFrame({'x': [3, 2, 1]}).astype('M8[ns]'), {}, {}),
                      reason="Parquet doesn't support nanosecond precision"),
    (pd.DataFrame({'x': [3, 2, 1]}).astype('M8[us]'), {}, {}),
    (pd.DataFrame({'x': [3, 2, 1]}).astype('M8[ms]'), {}, {}),
    (pd.DataFrame({'x': [3, 2, 1]}).astype('uint16'), {}, {}),
    (pd.DataFrame({'x': [3, 2, 1]}).astype('float32'), {}, {}),
    (pd.DataFrame({'x': [3, 1, 2]}, index=[3, 2, 1]), {}, {}),
    (pd.DataFrame({'x': [3, 1, 5]}, index=pd.Index([1, 2, 3], name='foo')), {}, {}),
    (pd.DataFrame({'x': [1, 2, 3],
                   'y': [3, 2, 1]}), {}, {}),
    (pd.DataFrame({'x': [1, 2, 3],
                   'y': [3, 2, 1]}, columns=['y', 'x']), {}, {}),
    (pd.DataFrame({'0': [3, 2, 1]}), {}, {}),
    (pd.DataFrame({'x': [3, 2, None]}), {}, {}),
    (pd.DataFrame({'-': [3., 2., None]}), {}, {}),
    (pd.DataFrame({'.': [3., 2., None]}), {}, {}),
    (pd.DataFrame({' ': [3., 2., None]}), {}, {}),
])
def test_roundtrip(tmpdir, df, write_kwargs, read_kwargs):
    check_fastparquet()
    tmp = str(tmpdir)
    if df.index.name is None:
        df.index.name = 'index'
    ddf = dd.from_pandas(df, npartitions=2)

    dd.to_parquet(ddf, tmp, **write_kwargs)
    ddf2 = dd.read_parquet(tmp, index=df.index.name, **read_kwargs)
    assert_eq(ddf, ddf2)


def test_categories(tmpdir):
    check_fastparquet()
    fn = str(tmpdir)
    df = pd.DataFrame({'x': [1, 2, 3, 4, 5],
                       'y': list('caaab')})
    ddf = dd.from_pandas(df, npartitions=2)
    ddf['y'] = ddf.y.astype('category')
    ddf.to_parquet(fn)
    ddf2 = dd.read_parquet(fn, categories=['y'])
    with pytest.raises(NotImplementedError):
        ddf2.y.cat.categories
    assert set(ddf2.y.compute().cat.categories) == {'a', 'b', 'c'}
    cats_set = ddf2.map_partitions(lambda x: x.y.cat.categories).compute()
    assert cats_set.tolist() == ['a', 'c', 'a', 'b']
    assert_eq(ddf.y, ddf2.y, check_names=False)
    with pytest.raises(TypeError):
        # attempt to load as category that which is not so encoded
        ddf2 = dd.read_parquet(fn, categories=['x']).compute()

    with pytest.raises(ValueError):
        # attempt to load as category unknown column
        ddf2 = dd.read_parquet(fn, categories=['foo'])


def test_empty_partition(tmpdir, engine):
    fn = str(tmpdir)
    df = pd.DataFrame({"a": range(10), "b": range(10)})
    ddf = dd.from_pandas(df, npartitions=5)

    ddf2 = ddf[ddf.a <= 5]
    ddf2.to_parquet(fn, engine=engine)

    ddf3 = dd.read_parquet(fn, engine=engine)
    sol = ddf2.compute()
    assert_eq(sol, ddf3, check_names=False, check_index=False)


def test_timestamp_index(tmpdir, engine):
    fn = str(tmpdir)
    df = tm.makeTimeDataFrame()
    df.index.name = 'foo'
    ddf = dd.from_pandas(df, npartitions=5)
    ddf.to_parquet(fn, engine=engine)
    ddf2 = dd.read_parquet(fn, engine=engine, infer_divisions=should_check_divs(engine))
    assert_eq(ddf, ddf2, check_divisions=should_check_divs(engine))


def test_to_parquet_default_writes_nulls(tmpdir):
    check_fastparquet()
    check_pyarrow()
    fn = str(tmpdir.join('test.parquet'))

    df = pd.DataFrame({'c1': [1., np.nan, 2, np.nan, 3]})
    ddf = dd.from_pandas(df, npartitions=1)

    ddf.to_parquet(fn)
    table = pq.read_table(fn)
    assert table[1].null_count == 2


@write_read_engines(
    xfail_pyarrow_fastparquet=pyarrow_fastparquet_msg,
    xfail_pyarrow_pyarrow=("Race condition writing using pyarrow with partition_on. "
                           "Fixed on master, but not on pyarrow 0.8.0")
)
def test_partition_on(tmpdir, write_engine, read_engine):
    tmpdir = str(tmpdir)
    df = pd.DataFrame({'a': np.random.choice(['A', 'B', 'C'], size=100),
                       'b': np.random.random(size=100),
                       'c': np.random.randint(1, 5, size=100)})
    d = dd.from_pandas(df, npartitions=2)
    d.to_parquet(tmpdir, partition_on=['a'], engine=write_engine)
    out = dd.read_parquet(tmpdir, engine=read_engine).compute()
    for val in df.a.unique():
        assert set(df.b[df.a == val]) == set(out.b[out.a == val])


def test_filters(tmpdir):
    check_fastparquet()
    fn = str(tmpdir)

    df = pd.DataFrame({'at': ['ab', 'aa', 'ba', 'da', 'bb']})
    ddf = dd.from_pandas(df, npartitions=1)

    # Ok with 1 partition and filters
    ddf.repartition(npartitions=1, force=True).to_parquet(fn, write_index=False)
    ddf2 = dd.read_parquet(fn, index=False,
                           filters=[('at', '==', 'aa')]).compute()
    assert_eq(ddf2, ddf)

    # with >1 partition and no filters
    ddf.repartition(npartitions=2, force=True).to_parquet(fn)
    dd.read_parquet(fn).compute()
    assert_eq(ddf2, ddf)

    # with >1 partition and filters using base fastparquet
    ddf.repartition(npartitions=2, force=True).to_parquet(fn)
    df2 = fastparquet.ParquetFile(fn).to_pandas(filters=[('at', '==', 'aa')])
    assert len(df2) > 0

    # with >1 partition and filters
    ddf.repartition(npartitions=2, force=True).to_parquet(fn)
    dd.read_parquet(fn, filters=[('at', '==', 'aa')]).compute()
    assert len(ddf2) > 0


@pytest.mark.parametrize('get', [dask.threaded.get, dask.multiprocessing.get])
def test_to_parquet_lazy(tmpdir, get, engine):
    tmpdir = str(tmpdir)
    df = pd.DataFrame({'a': [1, 2, 3, 4],
                       'b': [1., 2., 3., 4.]})
    df.index.name = 'index'
    ddf = dd.from_pandas(df, npartitions=2)
    value = ddf.to_parquet(tmpdir, compute=False, engine=engine)

    assert hasattr(value, 'dask')
    value.compute(get=get)
    assert os.path.exists(tmpdir)

    ddf2 = dd.read_parquet(tmpdir, engine=engine, infer_divisions=should_check_divs(engine))

    assert_eq(ddf, ddf2, check_divisions=should_check_divs(engine))


def test_timestamp96(tmpdir):
    check_fastparquet()
    fn = str(tmpdir)
    df = pd.DataFrame({'a': ['now']}, dtype='M8[ns]')
    ddf = dd.from_pandas(df, 1)
    ddf.to_parquet(fn, write_index=False, times='int96')
    pf = fastparquet.ParquetFile(fn)
    assert pf._schema[1].type == fastparquet.parquet_thrift.Type.INT96
    out = dd.read_parquet(fn).compute()
    assert_eq(out, df)


def test_drill_scheme(tmpdir):
    check_fastparquet()
    fn = str(tmpdir)
    N = 5
    df1 = pd.DataFrame({c: np.random.random(N)
                        for i, c in enumerate(['a', 'b', 'c'])})
    df2 = pd.DataFrame({c: np.random.random(N)
                        for i, c in enumerate(['a', 'b', 'c'])})
    files = []
    for d in ['test_data1', 'test_data2']:
        dn = os.path.join(fn, d)
        if not os.path.exists(dn):
            os.mkdir(dn)
        files.append(os.path.join(dn, 'data1.parq'))

    fastparquet.write(files[0], df1)
    fastparquet.write(files[1], df2)

    df = dd.read_parquet(files)
    assert 'dir0' in df.columns
    out = df.compute()
    assert 'dir0' in out
    assert (np.unique(out.dir0) == ['test_data1', 'test_data2']).all()


def test_parquet_select_cats(tmpdir):
    check_fastparquet()
    fn = str(tmpdir)
    df = pd.DataFrame({
        'categories': pd.Series(
            np.random.choice(['a', 'b', 'c', 'd', 'e', 'f'], size=100),
            dtype='category'),
        'ints': pd.Series(list(range(0, 100)), dtype='int'),
        'floats': pd.Series(list(range(0, 100)), dtype='float')})

    ddf = dd.from_pandas(df, 1)
    ddf.to_parquet(fn)
    rddf = dd.read_parquet(fn, columns=['ints'])
    assert list(rddf.columns) == ['ints']
    rddf = dd.read_parquet(fn)
    assert list(rddf.columns) == list(df)


@write_read_engines(
    xfail_pyarrow_fastparquet=pyarrow_fastparquet_msg,
    xfail_fastparquet_pyarrow="fastparquet gh#251"
)
def test_columns_name(tmpdir, write_engine, read_engine):
    if write_engine == 'fastparquet':
        pytest.skip('Fastparquet does not write column_indexes')

    if write_engine == 'pyarrow':
        import pyarrow as pa
        if pa.__version__ < LooseVersion('0.8.0'):
            pytest.skip("pyarrow<0.8.0 did not write column_indexes")

    df = pd.DataFrame({"A": [1, 2]}, index=pd.Index(['a', 'b'], name='idx'))
    df.columns.name = "cols"
    ddf = dd.from_pandas(df, 2)

    tmp = str(tmpdir)

    ddf.to_parquet(tmp, engine=write_engine)
    result = dd.read_parquet(tmp, engine=read_engine)
    assert_eq(result, df)


@pytest.mark.parametrize('compression,', ['default', None, 'gzip', 'snappy'])
def test_writing_parquet_with_compression(tmpdir, compression, engine):
    fn = str(tmpdir)

    if engine == 'fastparquet' and compression in ['snappy', 'default']:
        pytest.importorskip('snappy')

    df = pd.DataFrame({'x': ['a', 'b', 'c'] * 10,
                       'y': [1, 2, 3] * 10})
    ddf = dd.from_pandas(df, npartitions=3)

    ddf.to_parquet(fn, compression=compression, engine=engine)
    if engine == 'fastparquet' and compression == 'default':
        # ensure default compression for fastparquet is Snappy
        import fastparquet
        pf = fastparquet.ParquetFile(fn)
        assert pf.row_groups[0].columns[0].meta_data.codec == 1

    out = dd.read_parquet(fn, engine=engine, infer_divisions=should_check_divs(engine))
    assert_eq(out, ddf, check_index=(engine != 'fastparquet'), check_divisions=should_check_divs(engine))


@pytest.fixture(params=[
    # fastparquet 0.1.3
    {'columns': [{'metadata': None,
                  'name': 'idx',
                  'numpy_type': 'int64',
                  'pandas_type': 'int64'},
                 {'metadata': None,
                  'name': 'A',
                  'numpy_type': 'int64',
                  'pandas_type': 'int64'}],
     'index_columns': ['idx'],
     'pandas_version': '0.21.0'},

    # pyarrow 0.7.1
    {'columns': [{'metadata': None,
                  'name': 'A',
                  'numpy_type': 'int64',
                  'pandas_type': 'int64'},
                 {'metadata': None,
                  'name': 'idx',
                  'numpy_type': 'int64',
                  'pandas_type': 'int64'}],
     'index_columns': ['idx'],
     'pandas_version': '0.21.0'},

    # pyarrow 0.8.0
    {'column_indexes': [{'field_name': None,
                         'metadata': {'encoding': 'UTF-8'},
                         'name': None,
                         'numpy_type': 'object',
                         'pandas_type': 'unicode'}],
     'columns': [{'field_name': 'A',
                  'metadata': None,
                  'name': 'A',
                  'numpy_type': 'int64',
                  'pandas_type': 'int64'},
                 {'field_name': '__index_level_0__',
                  'metadata': None,
                  'name': 'idx',
                  'numpy_type': 'int64',
                  'pandas_type': 'int64'}],
     'index_columns': ['__index_level_0__'],
     'pandas_version': '0.21.0'},

    # TODO: fastparquet update
])
def pandas_metadata(request):
    return request.param


def test_parse_pandas_metadata(pandas_metadata):
    index_names, column_names, mapping, column_index_names = (
        _parse_pandas_metadata(pandas_metadata)
    )
    assert index_names == ['idx']
    assert column_names == ['A']
    assert column_index_names == [None]

    # for new pyarrow
    if pandas_metadata['index_columns'] == ['__index_level_0__']:
        assert mapping == {'__index_level_0__': 'idx', 'A': 'A'}
    else:
        assert mapping == {'idx': 'idx', 'A': 'A'}

    assert isinstance(mapping, dict)


def test_parse_pandas_metadata_null_index():
    # pyarrow 0.7.1 None for index
    e_index_names = [None]
    e_column_names = ['x']
    e_mapping = {'__index_level_0__': None, 'x': 'x'}
    e_column_index_names = [None]

    md = {'columns': [{'metadata': None,
                       'name': 'x',
                       'numpy_type': 'int64',
                       'pandas_type': 'int64'},
                      {'metadata': None,
                       'name': '__index_level_0__',
                       'numpy_type': 'int64',
                       'pandas_type': 'int64'}],
          'index_columns': ['__index_level_0__'],
          'pandas_version': '0.21.0'}
    index_names, column_names, mapping, column_index_names = (
        _parse_pandas_metadata(md)
    )
    assert index_names == e_index_names
    assert column_names == e_column_names
    assert mapping == e_mapping
    assert column_index_names == e_column_index_names

    # pyarrow 0.8.0 None for index
    md = {'column_indexes': [{'field_name': None,
                              'metadata': {'encoding': 'UTF-8'},
                              'name': None,
                              'numpy_type': 'object',
                              'pandas_type': 'unicode'}],
          'columns': [{'field_name': 'x',
                       'metadata': None,
                       'name': 'x',
                       'numpy_type': 'int64',
                       'pandas_type': 'int64'},
                      {'field_name': '__index_level_0__',
                       'metadata': None,
                       'name': None,
                       'numpy_type': 'int64',
                       'pandas_type': 'int64'}],
          'index_columns': ['__index_level_0__'],
          'pandas_version': '0.21.0'}
    index_names, column_names, mapping, column_index_names = (
        _parse_pandas_metadata(md)
    )
    assert index_names == e_index_names
    assert column_names == e_column_names
    assert mapping == e_mapping
    assert column_index_names == e_column_index_names


def test_pyarrow_raises_filters_categoricals(tmpdir):
    check_pyarrow()
    tmp = str(tmpdir)
    data = pd.DataFrame({"A": [1, 2]})
    df = dd.from_pandas(data, npartitions=2)

    df.to_parquet(tmp, write_index=False, engine="pyarrow")

    with pytest.raises(NotImplementedError):
        dd.read_parquet(tmp, engine="pyarrow", filters=["A>1"])


def test_read_no_metadata(tmpdir, engine):
    # use pyarrow.parquet to create a parquet file without
    # pandas metadata
    pa = pytest.importorskip("pyarrow")
    import pyarrow.parquet as pq
    tmp = str(tmpdir) + "table.parq"

    table = pa.Table.from_arrays([pa.array([1, 2, 3]),
                                  pa.array([3, 4, 5])],
                                 names=['A', 'B'])
    pq.write_table(table, tmp)
    result = dd.read_parquet(tmp, engine=engine)
    expected = pd.DataFrame({"A": [1, 2, 3], "B": [3, 4, 5]})
    assert_eq(result, expected)


def test_parse_pandas_metadata_duplicate_index_columns():
    md = {
        'column_indexes': [{
            'field_name': None,
            'metadata': {
                'encoding': 'UTF-8'
            },
            'name': None,
            'numpy_type': 'object',
            'pandas_type': 'unicode'
        }],
        'columns': [{
            'field_name': 'A',
            'metadata': None,
            'name': 'A',
            'numpy_type': 'int64',
            'pandas_type': 'int64'
        }, {
            'field_name': '__index_level_0__',
            'metadata': None,
            'name': 'A',
            'numpy_type': 'object',
            'pandas_type': 'unicode'
        }],
        'index_columns': ['__index_level_0__'],
        'pandas_version': '0.21.0'
    }
    index_names, column_names, storage_name_mapping, column_index_names = (
        _parse_pandas_metadata(md)
    )
    assert index_names == ['A']
    assert column_names == ['A']
    assert storage_name_mapping == {'__index_level_0__': 'A', 'A': 'A'}
    assert column_index_names == [None]


def test_parse_pandas_metadata_column_with_index_name():
    md = {
        'column_indexes': [{
            'field_name': None,
            'metadata': {
                'encoding': 'UTF-8'
            },
            'name': None,
            'numpy_type': 'object',
            'pandas_type': 'unicode'
        }],
        'columns': [{
            'field_name': 'A',
            'metadata': None,
            'name': 'A',
            'numpy_type': 'int64',
            'pandas_type': 'int64'
        }, {
            'field_name': '__index_level_0__',
            'metadata': None,
            'name': 'A',
            'numpy_type': 'object',
            'pandas_type': 'unicode'
        }],
        'index_columns': ['__index_level_0__'],
        'pandas_version': '0.21.0'
    }
    index_names, column_names, storage_name_mapping, column_index_names = (
        _parse_pandas_metadata(md)
    )
    assert index_names == ['A']
    assert column_names == ['A']
    assert storage_name_mapping == {'__index_level_0__': 'A', 'A': 'A'}
    assert column_index_names == [None]


def test_writing_parquet_with_kwargs(tmpdir, engine):
    fn = str(tmpdir)
    path1 = os.path.join(fn, 'normal')
    path2 = os.path.join(fn, 'partitioned')

    df = pd.DataFrame({'a': np.random.choice(['A', 'B', 'C'], size=100),
                       'b': np.random.random(size=100),
                       'c': np.random.randint(1, 5, size=100)})
    ddf = dd.from_pandas(df, npartitions=3)

    engine_kwargs = {
        'pyarrow': {
            'compression': 'snappy',
            'coerce_timestamps': None,
            'use_dictionary': True
        },
        'fastparquet': {
            'compression': 'snappy',
            'times': 'int64',
            'fixed_text': None
        }
    }

    ddf.to_parquet(path1,  engine=engine, **engine_kwargs[engine])
    out = dd.read_parquet(path1, engine=engine, infer_divisions=should_check_divs(engine))
    assert_eq(out, ddf, check_index=(engine != 'fastparquet'), check_divisions=should_check_divs(engine))

    # Avoid race condition in pyarrow 0.8.0 on writing partitioned datasets
    with dask.set_options(get=dask.get):
        ddf.to_parquet(path2, engine=engine, partition_on=['a'],
                       **engine_kwargs[engine])
    out = dd.read_parquet(path2, engine=engine).compute()
    for val in df.a.unique():
        assert set(df.b[df.a == val]) == set(out.b[out.a == val])


def test_writing_parquet_with_unknown_kwargs(tmpdir, engine):
    fn = str(tmpdir)

    with pytest.raises(TypeError):
        ddf.to_parquet(fn,  engine=engine, unknown_key='unknown_value')


def test_select_partitioned_column(tmpdir, engine):
    if engine == 'pyarrow':
        import pyarrow as pa
        if pa.__version__ < LooseVersion('0.9.0'):
            pytest.skip("pyarrow<0.9.0 did not support this")

    fn = str(tmpdir)
    size = 20
    d = {'signal1': np.random.normal(0, 0.3, size=size).cumsum() + 50,
         'fake_categorical1': np.random.choice(['A', 'B', 'C'], size=size),
         'fake_categorical2': np.random.choice(['D', 'E', 'F'], size=size)}
    df = dd.from_pandas(pd.DataFrame(d), 2)
    df.to_parquet(fn, compression='snappy', write_index=False, engine=engine,
                  partition_on=['fake_categorical1', 'fake_categorical2'])

    df_partitioned = dd.read_parquet(fn, engine=engine)
    df_partitioned[df_partitioned.fake_categorical1 == 'A'].compute()
