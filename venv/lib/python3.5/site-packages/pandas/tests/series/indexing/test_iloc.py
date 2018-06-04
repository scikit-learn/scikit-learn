# coding=utf-8
# pylint: disable-msg=E1101,W0612

import numpy as np

from pandas import Series

from pandas.compat import lrange, range
from pandas.util.testing import (assert_series_equal,
                                 assert_almost_equal)


def test_iloc():
    s = Series(np.random.randn(10), index=lrange(0, 20, 2))

    for i in range(len(s)):
        result = s.iloc[i]
        exp = s[s.index[i]]
        assert_almost_equal(result, exp)

    # pass a slice
    result = s.iloc[slice(1, 3)]
    expected = s.loc[2:4]
    assert_series_equal(result, expected)

    # test slice is a view
    result[:] = 0
    assert (s[1:3] == 0).all()

    # list of integers
    result = s.iloc[[0, 2, 3, 4, 5]]
    expected = s.reindex(s.index[[0, 2, 3, 4, 5]])
    assert_series_equal(result, expected)


def test_iloc_nonunique():
    s = Series([0, 1, 2], index=[0, 1, 0])
    assert s.iloc[2] == 2
