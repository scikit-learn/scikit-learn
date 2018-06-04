import numpy as np
import pandas as pd
import pandas.util.testing as tm

import pytest

from dask.dataframe.hashing import hash_pandas_object
from dask.dataframe.utils import assert_eq


@pytest.mark.parametrize('obj', [
    pd.Series([1, 2, 3]),
    pd.Series([1.0, 1.5, 3.2]),
    pd.Series([1.0, 1.5, 3.2], index=[1.5, 1.1, 3.3]),
    pd.Series(['a', 'b', 'c']),
    pd.Series([True, False, True]),
    pd.Index([1, 2, 3]),
    pd.Index([True, False, True]),
    pd.DataFrame({'x': ['a', 'b', 'c'], 'y': [1, 2, 3]}),
    pd.util.testing.makeMissingDataframe(),
    pd.util.testing.makeMixedDataFrame(),
    pd.util.testing.makeTimeDataFrame(),
    pd.util.testing.makeTimeSeries(),
    pd.util.testing.makeTimedeltaIndex()])
def test_hash_pandas_object(obj):
    a = hash_pandas_object(obj)
    b = hash_pandas_object(obj)
    if isinstance(a, np.ndarray):
        np.testing.assert_equal(a, b)
    else:
        assert_eq(a, b)


def test_categorical_consistency():
    # Check that categoricals hash consistent with their values, not codes
    # This should work for categoricals of any dtype
    for s1 in [pd.Series(['a', 'b', 'c', 'd']),
               pd.Series([1000, 2000, 3000, 4000]),
               pd.Series(pd.date_range(0, periods=4))]:
        s2 = s1.astype('category').cat.set_categories(s1)
        s3 = s2.cat.set_categories(list(reversed(s1)))
        for categorize in [True, False]:
            # These should all hash identically
            h1 = hash_pandas_object(s1, categorize=categorize)
            h2 = hash_pandas_object(s2, categorize=categorize)
            h3 = hash_pandas_object(s3, categorize=categorize)
            tm.assert_series_equal(h1, h2)
            tm.assert_series_equal(h1, h3)


def test_object_missing_values():
    # Check that the presence of missing values doesn't change how object dtype
    # is hashed.
    s = pd.Series(['a', 'b', 'c', None])
    h1 = hash_pandas_object(s).iloc[:3]
    h2 = hash_pandas_object(s.iloc[:3])
    tm.assert_series_equal(h1, h2)
