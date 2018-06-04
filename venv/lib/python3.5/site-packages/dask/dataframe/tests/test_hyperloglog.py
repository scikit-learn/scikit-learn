
import dask
import dask.dataframe as dd

import pandas as pd
import numpy as np
import pytest


rs = np.random.RandomState(96)


@pytest.mark.parametrize("df", [
    pd.DataFrame({
        'x': [1, 2, 3] * 3,
        'y': [1.2, 3.4, 5.6] * 3,
        'z': -np.arange(9, dtype=np.int8)}),
    pd.DataFrame({
        'x': rs.randint(0, 1000000, (10000,)),
        'y': rs.randn(10000),
        'z': rs.uniform(0, 9999999, (10000,))}),
    pd.DataFrame({
        'x': np.repeat(rs.randint(0, 1000000, (1000,)), 3),
        'y': np.repeat(rs.randn(1000), 3),
        'z': np.repeat(rs.uniform(0, 9999999, (1000,)), 3)}),
    pd.DataFrame({
        'x': rs.randint(0, 1000000, (10000,))}),
    pd.DataFrame({
        'x': rs.randint(0, 1000000, (7,)),
        'y': ['a', 'bet', 'is', 'a', 'tax', 'on', 'bs']}),
    pd.DataFrame({
        'w': np.zeros((20000,)),
        'x': np.zeros((20000,)),
        'y': np.zeros((20000,)) + 4803592,
        'z': np.zeros((20000,))}),
    pd.DataFrame({'x': [1, 2, 3] * 1000}),
    pd.DataFrame({'x': np.random.random(1000)}),
    pd.DataFrame({
        'a': [1, 2, 3] * 3,
        'b': [1.2, 3.4, 5.6] * 3,
        'c': [1 + 2j, 3 + 4j, 5 + 6j] * 3,
        'd': -np.arange(9, dtype=np.int8)}),
    pd.Series([1, 2, 3] * 1000),
    pd.Series(np.random.random(1000)),
    pd.Series(np.random.random(1000), index=np.ones(1000)),
    pd.Series(np.random.random(1000), index=np.random.random(1000)),
])
@pytest.mark.parametrize('npartitions', [2, 20])
def test_basic(df, npartitions):
    ddf = dd.from_pandas(df, npartitions=npartitions)

    approx = ddf.nunique_approx().compute(get=dask.local.get_sync)
    exact = len(df.drop_duplicates())
    assert abs(approx - exact) <= 2 or abs(approx - exact) / exact < 0.05


@pytest.mark.parametrize('split_every', [None, 2, 10])
@pytest.mark.parametrize('npartitions', [2, 20])
def test_split_every(split_every, npartitions):
    df = pd.Series([1, 2, 3] * 1000)
    ddf = dd.from_pandas(df, npartitions=npartitions)

    approx = ddf.nunique_approx(split_every=split_every).compute(get=dask.local.get_sync)
    exact = len(df.drop_duplicates())
    assert abs(approx - exact) <= 2 or abs(approx - exact) / exact < 0.05


def test_larger_data():
    df = dd.demo.make_timeseries('2000-01-01', '2000-04-01',
                                 {'value': float, 'id': int},
                                 freq='10s', partition_freq='1D', seed=1)
    assert df.nunique_approx().compute() > 1000
