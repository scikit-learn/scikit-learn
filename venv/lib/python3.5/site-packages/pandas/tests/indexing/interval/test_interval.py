import pytest
import numpy as np
import pandas as pd

from pandas import Series, DataFrame, IntervalIndex, Interval
from pandas.compat import product
import pandas.util.testing as tm


class TestIntervalIndex(object):

    def setup_method(self, method):
        self.s = Series(np.arange(5), IntervalIndex.from_breaks(np.arange(6)))

    # To be removed, replaced by test_interval_new.py (see #16316, #16386)
    def test_loc_with_scalar(self):

        s = self.s

        expected = s.iloc[:3]
        tm.assert_series_equal(expected, s.loc[:3])
        tm.assert_series_equal(expected, s.loc[:2.5])
        tm.assert_series_equal(expected, s.loc[0.1:2.5])
        tm.assert_series_equal(expected, s.loc[-1:3])

        expected = s.iloc[1:4]
        tm.assert_series_equal(expected, s.loc[[1.5, 2.5, 3.5]])
        tm.assert_series_equal(expected, s.loc[[2, 3, 4]])
        tm.assert_series_equal(expected, s.loc[[1.5, 3, 4]])

        expected = s.iloc[2:5]
        tm.assert_series_equal(expected, s.loc[s >= 2])

    # TODO: check this behavior is consistent with test_interval_new.py
    def test_getitem_with_scalar(self):

        s = self.s

        expected = s.iloc[:3]
        tm.assert_series_equal(expected, s[:3])
        tm.assert_series_equal(expected, s[:2.5])
        tm.assert_series_equal(expected, s[0.1:2.5])
        tm.assert_series_equal(expected, s[-1:3])

        expected = s.iloc[1:4]
        tm.assert_series_equal(expected, s[[1.5, 2.5, 3.5]])
        tm.assert_series_equal(expected, s[[2, 3, 4]])
        tm.assert_series_equal(expected, s[[1.5, 3, 4]])

        expected = s.iloc[2:5]
        tm.assert_series_equal(expected, s[s >= 2])

    # TODO: check this behavior is consistent with test_interval_new.py
    @pytest.mark.parametrize('direction, closed',
                             product(('increasing', 'decreasing'),
                                     ('left', 'right', 'neither', 'both')))
    def test_nonoverlapping_monotonic(self, direction, closed):
        tpls = [(0, 1), (2, 3), (4, 5)]
        if direction == 'decreasing':
            tpls = tpls[::-1]

        idx = IntervalIndex.from_tuples(tpls, closed=closed)
        s = Series(list('abc'), idx)

        for key, expected in zip(idx.left, s):
            if idx.closed_left:
                assert s[key] == expected
                assert s.loc[key] == expected
            else:
                with pytest.raises(KeyError):
                    s[key]
                with pytest.raises(KeyError):
                    s.loc[key]

        for key, expected in zip(idx.right, s):
            if idx.closed_right:
                assert s[key] == expected
                assert s.loc[key] == expected
            else:
                with pytest.raises(KeyError):
                    s[key]
                with pytest.raises(KeyError):
                    s.loc[key]

        for key, expected in zip(idx.mid, s):
            assert s[key] == expected
            assert s.loc[key] == expected

    # To be removed, replaced by test_interval_new.py (see #16316, #16386)
    def test_with_interval(self):

        s = self.s
        expected = 0

        result = s.loc[Interval(0, 1)]
        assert result == expected

        result = s[Interval(0, 1)]
        assert result == expected

        expected = s.iloc[3:5]
        result = s.loc[Interval(3, 6)]
        tm.assert_series_equal(expected, result)

        expected = s.iloc[3:5]
        result = s.loc[[Interval(3, 6)]]
        tm.assert_series_equal(expected, result)

        expected = s.iloc[3:5]
        result = s.loc[[Interval(3, 5)]]
        tm.assert_series_equal(expected, result)

        # missing
        with pytest.raises(KeyError):
            s.loc[Interval(-2, 0)]

        with pytest.raises(KeyError):
            s[Interval(-2, 0)]

        with pytest.raises(KeyError):
            s.loc[Interval(5, 6)]

        with pytest.raises(KeyError):
            s[Interval(5, 6)]

    # To be removed, replaced by test_interval_new.py (see #16316, #16386)
    def test_with_slices(self):

        s = self.s

        # slice of interval
        with pytest.raises(NotImplementedError):
            s.loc[Interval(3, 6):]

        with pytest.raises(NotImplementedError):
            s[Interval(3, 6):]

        expected = s.iloc[3:5]
        result = s[[Interval(3, 6)]]
        tm.assert_series_equal(expected, result)

        # slice of scalar with step != 1
        with pytest.raises(ValueError):
            s[0:4:2]

    # To be removed, replaced by test_interval_new.py (see #16316, #16386)
    def test_with_overlaps(self):

        s = self.s
        expected = s.iloc[[3, 4, 3, 4]]
        result = s.loc[[Interval(3, 6), Interval(3, 6)]]
        tm.assert_series_equal(expected, result)

        idx = IntervalIndex.from_tuples([(1, 5), (3, 7)])
        s = Series(range(len(idx)), index=idx)

        result = s[4]
        expected = s
        tm.assert_series_equal(expected, result)

        result = s[[4]]
        expected = s
        tm.assert_series_equal(expected, result)

        result = s.loc[[4]]
        expected = s
        tm.assert_series_equal(expected, result)

        result = s[Interval(3, 5)]
        expected = s
        tm.assert_series_equal(expected, result)

        result = s.loc[Interval(3, 5)]
        expected = s
        tm.assert_series_equal(expected, result)

        # doesn't intersect unique set of intervals
        with pytest.raises(KeyError):
            s[[Interval(3, 5)]]

        with pytest.raises(KeyError):
            s.loc[[Interval(3, 5)]]

    # To be removed, replaced by test_interval_new.py (see #16316, #16386)
    def test_non_unique(self):

        idx = IntervalIndex.from_tuples([(1, 3), (3, 7)])

        s = Series(range(len(idx)), index=idx)

        result = s.loc[Interval(1, 3)]
        assert result == 0

        result = s.loc[[Interval(1, 3)]]
        expected = s.iloc[0:1]
        tm.assert_series_equal(expected, result)

    # To be removed, replaced by test_interval_new.py (see #16316, #16386)
    def test_non_unique_moar(self):

        idx = IntervalIndex.from_tuples([(1, 3), (1, 3), (3, 7)])
        s = Series(range(len(idx)), index=idx)

        result = s.loc[Interval(1, 3)]
        expected = s.iloc[[0, 1]]
        tm.assert_series_equal(expected, result)

        # non-unique index and slices not allowed
        with pytest.raises(ValueError):
            s.loc[Interval(1, 3):]

        with pytest.raises(ValueError):
            s[Interval(1, 3):]

        # non-unique
        with pytest.raises(ValueError):
            s[[Interval(1, 3)]]

    # TODO: check this behavior is consistent with test_interval_new.py
    def test_non_matching(self):
        s = self.s

        # this is a departure from our current
        # indexin scheme, but simpler
        with pytest.raises(KeyError):
            s.loc[[-1, 3, 4, 5]]

        with pytest.raises(KeyError):
            s.loc[[-1, 3]]

    def test_large_series(self):
        s = Series(np.arange(1000000),
                   index=IntervalIndex.from_breaks(np.arange(1000001)))

        result1 = s.loc[:80000]
        result2 = s.loc[0:80000]
        result3 = s.loc[0:80000:1]
        tm.assert_series_equal(result1, result2)
        tm.assert_series_equal(result1, result3)

    def test_loc_getitem_frame(self):

        df = DataFrame({'A': range(10)})
        s = pd.cut(df.A, 5)
        df['B'] = s
        df = df.set_index('B')

        result = df.loc[4]
        expected = df.iloc[4:6]
        tm.assert_frame_equal(result, expected)

        with pytest.raises(KeyError):
            df.loc[10]

        # single list-like
        result = df.loc[[4]]
        expected = df.iloc[4:6]
        tm.assert_frame_equal(result, expected)

        # non-unique
        result = df.loc[[4, 5]]
        expected = df.take([4, 5, 4, 5])
        tm.assert_frame_equal(result, expected)

        with pytest.raises(KeyError):
            df.loc[[10]]

        # partial missing
        with pytest.raises(KeyError):
            df.loc[[10, 4]]
