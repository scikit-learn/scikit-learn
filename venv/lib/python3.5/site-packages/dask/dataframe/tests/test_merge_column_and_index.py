import dask.dataframe as dd
import numpy as np
import pandas as pd
import pytest

from dask.dataframe.utils import assert_eq, PANDAS_VERSION


# Fixtures
# ========
@pytest.fixture
def df_left():
    # Create frame with 10 partitions
    # Frame has 11 distinct idx values
    partition_sizes = np.array([3, 4, 2, 5, 3, 2, 5, 9, 4, 7, 4])
    idx = [i for i, s in enumerate(partition_sizes) for _ in range(s)]
    k = [i for s in partition_sizes for i in range(s)]
    vi = range(len(k))

    return pd.DataFrame(dict(
        idx=idx,
        k=k,
        v1=vi
    )).set_index(['idx'])


@pytest.fixture
def df_right():
    # Create frame with 10 partitions
    # Frame has 11 distinct idx values
    partition_sizes = np.array([4, 2, 5, 3, 2, 5, 9, 4, 7, 4, 8])
    idx = [i for i, s in enumerate(partition_sizes) for _ in range(s)]
    k = [i for s in partition_sizes for i in range(s)]
    vi = range(len(k))

    return pd.DataFrame(dict(
        idx=idx,
        k=k,
        v1=vi
    )).set_index(['idx'])


@pytest.fixture
def ddf_left(df_left):
    # Create frame with 10 partitions
    # Skip division on 2 so there is one mismatch with ddf_right
    return dd.repartition(df_left, [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11])


@pytest.fixture
def ddf_left_unknown(ddf_left):
    return ddf_left.clear_divisions()


@pytest.fixture
def ddf_left_single(df_left):
    return dd.from_pandas(df_left, npartitions=1, sort=False)


@pytest.fixture
def ddf_right(df_right):
    # Create frame with 10 partitions
    # Skip division on 3 so there is one mismatch with ddf_left
    return dd.repartition(df_right, [0, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11])


@pytest.fixture
def ddf_right_unknown(ddf_right):
    return ddf_right.clear_divisions()


@pytest.fixture
def ddf_right_single(df_right):
    return dd.from_pandas(df_right, npartitions=1, sort=False)


@pytest.fixture(params=['inner', 'left', 'right', 'outer'])
def how(request):
    return request.param


@pytest.fixture(params=[
    'idx',
    ['idx'],
    ['idx', 'k'],
    ['k', 'idx']])
def on(request):
    return request.param


# Tests
# =====
@pytest.mark.skipif(PANDAS_VERSION < '0.23.0',
                    reason="Need pandas col+index merge support (pandas-dev/pandas#14355)")
def test_merge_known_to_known(df_left, df_right, ddf_left, ddf_right, on, how):
    # Compute expected
    expected = df_left.merge(df_right, on=on, how=how)

    # Perform merge
    result = ddf_left.merge(ddf_right, on=on, how=how, shuffle='tasks')

    # Assertions
    assert_eq(result, expected)
    assert_eq(result.divisions, tuple(range(12)))
    assert len(result.__dask_graph__()) < 80


@pytest.mark.skipif(PANDAS_VERSION < '0.23.0',
                    reason="Need pandas col+index merge support (pandas-dev/pandas#14355)")
@pytest.mark.parametrize('how', ['inner', 'left'])
def test_merge_known_to_single(df_left, df_right, ddf_left, ddf_right_single, on, how):
    # Compute expected
    expected = df_left.merge(df_right, on=on, how=how)

    # Perform merge
    result = ddf_left.merge(ddf_right_single, on=on, how=how, shuffle='tasks')

    # Assertions
    assert_eq(result, expected)
    assert_eq(result.divisions, ddf_left.divisions)
    assert len(result.__dask_graph__()) < 30


@pytest.mark.skipif(PANDAS_VERSION < '0.23.0',
                    reason="Need pandas col+index merge support (pandas-dev/pandas#14355)")
@pytest.mark.parametrize('how', ['inner', 'right'])
def test_merge_single_to_known(df_left, df_right, ddf_left_single, ddf_right, on, how):
    # Compute expected
    expected = df_left.merge(df_right, on=on, how=how)

    # Perform merge
    result = ddf_left_single.merge(ddf_right, on=on, how=how, shuffle='tasks')

    # Assertions
    assert_eq(result, expected)
    assert_eq(result.divisions, ddf_right.divisions)
    assert len(result.__dask_graph__()) < 30


@pytest.mark.skipif(PANDAS_VERSION < '0.23.0',
                    reason="Need pandas col+index merge support (pandas-dev/pandas#14355)")
def test_merge_known_to_unknown(df_left, df_right, ddf_left, ddf_right_unknown, on, how):
    # Compute expected
    expected = df_left.merge(df_right, on=on, how=how)

    # Perform merge
    result = ddf_left.merge(ddf_right_unknown, on=on, how=how, shuffle='tasks')

    # Assertions
    assert_eq(result, expected)
    assert_eq(result.divisions, tuple(None for _ in range(11)))
    assert len(result.__dask_graph__()) >= 400


@pytest.mark.skipif(PANDAS_VERSION < '0.23.0',
                    reason="Need pandas col+index merge support (pandas-dev/pandas#14355)")
def test_merge_unknown_to_known(df_left, df_right, ddf_left_unknown, ddf_right, on, how):
    # Compute expected
    expected = df_left.merge(df_right, on=on, how=how)

    # Perform merge
    result = ddf_left_unknown.merge(ddf_right, on=on, how=how, shuffle='tasks')

    # Assertions
    assert_eq(result, expected)
    assert_eq(result.divisions, tuple(None for _ in range(11)))
    assert len(result.__dask_graph__()) > 400


@pytest.mark.skipif(PANDAS_VERSION < '0.23.0',
                    reason="Need pandas col+index merge support (pandas-dev/pandas#14355)")
def test_merge_unknown_to_unknown(df_left, df_right, ddf_left_unknown, ddf_right_unknown, on, how):
    # Compute expected
    expected = df_left.merge(df_right, on=on, how=how)

    # Merge unknown to unknown
    result = ddf_left_unknown.merge(ddf_right_unknown, on=on, how=how, shuffle='tasks')

    # Assertions
    assert_eq(result, expected)
    assert_eq(result.divisions, tuple(None for _ in range(11)))
    assert len(result.__dask_graph__()) > 400
