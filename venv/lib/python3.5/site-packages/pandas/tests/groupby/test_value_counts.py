"""
these are systematically testing all of the args to value_counts
with different size combinations. This is to ensure stability of the sorting
and proper parameter handling
"""

import pytest

from itertools import product
import numpy as np

from pandas.util import testing as tm
from pandas import MultiIndex, DataFrame, Series, date_range


# our starting frame
def seed_df(seed_nans, n, m):
    np.random.seed(1234)
    days = date_range('2015-08-24', periods=10)

    frame = DataFrame({
        '1st': np.random.choice(
            list('abcd'), n),
        '2nd': np.random.choice(days, n),
        '3rd': np.random.randint(1, m + 1, n)
    })

    if seed_nans:
        frame.loc[1::11, '1st'] = np.nan
        frame.loc[3::17, '2nd'] = np.nan
        frame.loc[7::19, '3rd'] = np.nan
        frame.loc[8::19, '3rd'] = np.nan
        frame.loc[9::19, '3rd'] = np.nan

    return frame


# create input df, keys, and the bins
binned = []
ids = []
for seed_nans in [True, False]:
    for n, m in product((100, 1000), (5, 20)):

        df = seed_df(seed_nans, n, m)
        bins = None, np.arange(0, max(5, df['3rd'].max()) + 1, 2)
        keys = '1st', '2nd', ['1st', '2nd']
        for k, b in product(keys, bins):
            binned.append((df, k, b, n, m))
            ids.append("{}-{}-{}".format(k, n, m))


@pytest.mark.slow
@pytest.mark.parametrize("df, keys, bins, n, m", binned, ids=ids)
def test_series_groupby_value_counts(df, keys, bins, n, m):

    def rebuild_index(df):
        arr = list(map(df.index.get_level_values, range(df.index.nlevels)))
        df.index = MultiIndex.from_arrays(arr, names=df.index.names)
        return df

    for isort, normalize, sort, ascending, dropna \
            in product((False, True), repeat=5):

        kwargs = dict(normalize=normalize, sort=sort,
                      ascending=ascending, dropna=dropna, bins=bins)

        gr = df.groupby(keys, sort=isort)
        left = gr['3rd'].value_counts(**kwargs)

        gr = df.groupby(keys, sort=isort)
        right = gr['3rd'].apply(Series.value_counts, **kwargs)
        right.index.names = right.index.names[:-1] + ['3rd']

        # have to sort on index because of unstable sort on values
        left, right = map(rebuild_index, (left, right))  # xref GH9212
        tm.assert_series_equal(left.sort_index(), right.sort_index())
