import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest

import dask.dataframe as dd

from dask.dataframe.utils import assert_eq, make_meta


@pytest.mark.parametrize('data', [
    pd.Series([1, 1, 1, 2, 2, 1, 3, 4], dtype='category'),
    pd.Series(pd.Categorical([1, 1, 1, 2, 2, 1, 3, 4], categories=[4, 3, 2, 1])),
    pd.DataFrame({'a': [1, 2, 3, 4, 4, 3, 2, 1],
                  'b': pd.Categorical(list('abcdabcd'))})]
)
def test_get_dummies(data):
    exp = pd.get_dummies(data)

    ddata = dd.from_pandas(data, 2)
    res = dd.get_dummies(ddata)
    assert_eq(res, exp)
    tm.assert_index_equal(res.columns, exp.columns)


def test_get_dummies_object():
    df = pd.DataFrame({'a': pd.Categorical([1, 2, 3, 4, 4, 3, 2, 1]),
                       'b': list('abcdabcd'),
                       'c': pd.Categorical(list('abcdabcd'))})
    ddf = dd.from_pandas(df, 2)

    # Explicitly exclude object columns
    exp = pd.get_dummies(df, columns=['a', 'c'])
    res = dd.get_dummies(ddf, columns=['a', 'c'])
    assert_eq(res, exp)
    tm.assert_index_equal(res.columns, exp.columns)

    with pytest.raises(NotImplementedError):
        dd.get_dummies(ddf)

    with pytest.raises(NotImplementedError):
        dd.get_dummies(ddf.b)

    with pytest.raises(NotImplementedError):
        dd.get_dummies(ddf, columns=['b'])


def test_get_dummies_kwargs():
    s = pd.Series([1, 1, 1, 2, 2, 1, 3, 4], dtype='category')
    exp = pd.get_dummies(s, prefix='X', prefix_sep='-')

    ds = dd.from_pandas(s, 2)
    res = dd.get_dummies(ds, prefix='X', prefix_sep='-')
    assert_eq(res, exp)
    tm.assert_index_equal(res.columns, pd.Index(['X-1', 'X-2', 'X-3', 'X-4']))

    exp = pd.get_dummies(s, drop_first=True)

    ds = dd.from_pandas(s, 2)
    res = dd.get_dummies(ds, drop_first=True)
    assert_eq(res, exp)
    tm.assert_index_equal(res.columns, exp.columns)

    # nan
    s = pd.Series([1, 1, 1, 2, np.nan, 3, np.nan, 5], dtype='category')
    exp = pd.get_dummies(s)

    ds = dd.from_pandas(s, 2)
    res = dd.get_dummies(ds)
    assert_eq(res, exp)
    tm.assert_index_equal(res.columns, exp.columns)

    # dummy_na
    exp = pd.get_dummies(s, dummy_na=True)

    ds = dd.from_pandas(s, 2)
    res = dd.get_dummies(ds, dummy_na=True)
    assert_eq(res, exp)
    tm.assert_index_equal(res.columns, pd.Index([1, 2, 3, 5, np.nan]))

    msg = 'sparse=True is not supported'
    with pytest.raises(NotImplementedError) as err:
        dd.get_dummies(ds, sparse=True)
    assert msg in str(err.value)


def test_get_dummies_errors():
    with pytest.raises(NotImplementedError):
        # not Categorical
        s = pd.Series([1, 1, 1, 2, 2, 1, 3, 4])
        ds = dd.from_pandas(s, 2)
        dd.get_dummies(ds)

    # unknown categories
    df = pd.DataFrame({'x': list('abcbc'), 'y': list('bcbcb')})
    ddf = dd.from_pandas(df, npartitions=2)
    ddf._meta = make_meta({'x': 'category', 'y': 'category'})

    with pytest.raises(NotImplementedError):
        dd.get_dummies(ddf)

    with pytest.raises(NotImplementedError):
        dd.get_dummies(ddf, columns=['x', 'y'])

    with pytest.raises(NotImplementedError):
        dd.get_dummies(ddf.x)


@pytest.mark.parametrize('aggfunc', ['mean', 'sum', 'count'])
def test_pivot_table(aggfunc):
    df = pd.DataFrame({'A': np.random.choice(list('XYZ'), size=100),
                       'B': np.random.randn(100),
                       'C': pd.Categorical(np.random.choice(list('abc'), size=100))})
    ddf = dd.from_pandas(df, 5)

    res = dd.pivot_table(ddf, index='A', columns='C', values='B',
                         aggfunc=aggfunc)
    exp = pd.pivot_table(df, index='A', columns='C', values='B',
                         aggfunc=aggfunc)
    if aggfunc == 'count':
        # dask result cannot be int64 dtype depending on divisions because of NaN
        exp = exp.astype(np.float64)

    assert_eq(res, exp)

    # method
    res = ddf.pivot_table(index='A', columns='C', values='B',
                          aggfunc=aggfunc)
    exp = df.pivot_table(index='A', columns='C', values='B',
                         aggfunc=aggfunc)
    if aggfunc == 'count':
        # dask result cannot be int64 dtype depending on divisions because of NaN
        exp = exp.astype(np.float64)
    assert_eq(res, exp)


def test_pivot_table_dtype():

    df = pd.DataFrame({'A': list('AABB'),
                       'B': pd.Categorical(list('ABAB')),
                       'C': [1, 2, 3, 4]})
    ddf = dd.from_pandas(df, 2)
    res = dd.pivot_table(ddf, index='A', columns='B',
                         values='C', aggfunc='count')

    exp_index = pd.CategoricalIndex(['A', 'B'], name='B')
    exp = pd.Series([np.float64] * 2, index=exp_index)
    tm.assert_series_equal(res.dtypes, exp)

    exp = pd.pivot_table(df, index='A', columns='B',
                         values='C', aggfunc='count').astype(np.float64)

    assert_eq(res, exp)


def test_pivot_table_errors():
    df = pd.DataFrame({'A': np.random.choice(list('abc'), size=10),
                       'B': np.random.randn(10),
                       'C': pd.Categorical(np.random.choice(list('abc'), size=10))})
    ddf = dd.from_pandas(df, 2)

    msg = "'index' must be the name of an existing column"
    with pytest.raises(ValueError) as err:
        dd.pivot_table(ddf, index=['A'], columns='C', values='B')
    assert msg in str(err.value)
    msg = "'columns' must be the name of an existing column"
    with pytest.raises(ValueError) as err:
        dd.pivot_table(ddf, index='A', columns=['C'], values='B')
    assert msg in str(err.value)
    msg = "'values' must be the name of an existing column"
    with pytest.raises(ValueError) as err:
        dd.pivot_table(ddf, index='A', columns='C', values=['B'])
    assert msg in str(err.value)

    msg = "aggfunc must be either 'mean', 'sum' or 'count'"
    with pytest.raises(ValueError) as err:
        dd.pivot_table(ddf, index='A', columns='C', values='B', aggfunc=['sum'])
    assert msg in str(err.value)

    with pytest.raises(ValueError) as err:
        dd.pivot_table(ddf, index='A', columns='C', values='B', aggfunc='xx')
    assert msg in str(err.value)

    # unknown categories
    ddf._meta = make_meta({'A': object, 'B': float, 'C': 'category'})
    msg = "'columns' must have known categories"
    with pytest.raises(ValueError) as err:
        dd.pivot_table(ddf, index='A', columns='C', values=['B'])
    assert msg in str(err.value)

    df = pd.DataFrame({'A': np.random.choice(list('abc'), size=10),
                       'B': np.random.randn(10),
                       'C': np.random.choice(list('abc'), size=10)})
    ddf = dd.from_pandas(df, 2)
    msg = "'columns' must be category dtype"
    with pytest.raises(ValueError) as err:
        dd.pivot_table(ddf, index='A', columns='C', values='B')
    assert msg in str(err.value)
