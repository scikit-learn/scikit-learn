import pytest
import pandas as pd
import numpy as np

from pandas.util.testing import assert_frame_equal, assert_series_equal
import pandas.util.testing as tm


@pytest.fixture(params=[['inner'], ['inner', 'outer']])
def frame(request):
    levels = request.param
    df = pd.DataFrame({'outer': ['a', 'a', 'a', 'b', 'b', 'b'],
                       'inner': [1, 2, 3, 1, 2, 3],
                       'A': np.arange(6),
                       'B': ['one', 'one', 'two', 'two', 'one', 'one']})
    if levels:
        df = df.set_index(levels)

    return df


@pytest.fixture()
def series():
    df = pd.DataFrame({'outer': ['a', 'a', 'a', 'b', 'b', 'b'],
                       'inner': [1, 2, 3, 1, 2, 3],
                       'A': np.arange(6),
                       'B': ['one', 'one', 'two', 'two', 'one', 'one']})
    s = df.set_index(['outer', 'inner', 'B'])['A']

    return s


@pytest.mark.parametrize('key_strs,groupers', [
    ('inner',  # Index name
     pd.Grouper(level='inner')
     ),
    (['inner'],  # List of index name
     [pd.Grouper(level='inner')]
     ),
    (['B', 'inner'],  # Column and index
     ['B', pd.Grouper(level='inner')]
     ),
    (['inner', 'B'],  # Index and column
     [pd.Grouper(level='inner'), 'B'])])
def test_grouper_index_level_as_string(frame, key_strs, groupers):
    result = frame.groupby(key_strs).mean()
    expected = frame.groupby(groupers).mean()
    assert_frame_equal(result, expected)


@pytest.mark.parametrize('levels', [
    'inner', 'outer', 'B',
    ['inner'], ['outer'], ['B'],
    ['inner', 'outer'], ['outer', 'inner'],
    ['inner', 'outer', 'B'], ['B', 'outer', 'inner']
])
def test_grouper_index_level_as_string_series(series, levels):

    # Compute expected result
    if isinstance(levels, list):
        groupers = [pd.Grouper(level=lv) for lv in levels]
    else:
        groupers = pd.Grouper(level=levels)

    expected = series.groupby(groupers).mean()

    # Compute and check result
    result = series.groupby(levels).mean()
    assert_series_equal(result, expected)


@pytest.mark.parametrize('key_strs,key_groupers,level_groupers', [
    ('inner',  # Index name
     pd.Grouper(key='inner'),
     pd.Grouper(level='inner'),
     ),
    (['inner'],  # List of index name
     [pd.Grouper(key='inner')],
     [pd.Grouper(level='inner')]
     ),
    (['B', 'inner'],  # Column and index
     ['B', pd.Grouper(key='inner')],
     ['B', pd.Grouper(level='inner')]
     ),
    (['inner', 'B'],  # Index and column
     [pd.Grouper(key='inner'), 'B'],
     [pd.Grouper(level='inner'), 'B'])])
def test_grouper_column_index_level_precedence(frame,
                                               key_strs,
                                               key_groupers,
                                               level_groupers):

    # GH 5677, when a string passed as the `by` parameter
    # matches a column and an index level the column takes
    # precedence and a FutureWarning is raised

    # Add 'inner' column to frame
    # (frame already has an 'inner' index)
    frame['inner'] = [1, 1, 1, 1, 1, 1]

    # Performing a groupby with strings should produce warning
    with tm.assert_produces_warning(FutureWarning):
        result = frame.groupby(key_strs).mean()

    # Grouping with key Grouper should produce the same result and no warning
    with tm.assert_produces_warning(False):
        expected = frame.groupby(key_groupers).mean()

    assert_frame_equal(result, expected)

    # Grouping with level Grouper should produce a different result but
    # still no warning
    with tm.assert_produces_warning(False):
        not_expected = frame.groupby(level_groupers).mean()

    assert not result.index.equals(not_expected.index)
