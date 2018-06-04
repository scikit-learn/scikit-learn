import numpy as np

import pandas as pd
import pandas.util.testing as tm
import pandas.core.indexes.period as period
from pandas import Series, period_range, DataFrame


def _permute(obj):
    return obj.take(np.random.permutation(len(obj)))


class TestSeriesPeriod(object):

    def setup_method(self, method):
        self.series = Series(period_range('2000-01-01', periods=10, freq='D'))

    def test_auto_conversion(self):
        series = Series(list(period_range('2000-01-01', periods=10, freq='D')))
        assert series.dtype == 'object'

        series = pd.Series([pd.Period('2011-01-01', freq='D'),
                            pd.Period('2011-02-01', freq='D')])
        assert series.dtype == 'object'

    def test_getitem(self):
        assert self.series[1] == pd.Period('2000-01-02', freq='D')

        result = self.series[[2, 4]]
        exp = pd.Series([pd.Period('2000-01-03', freq='D'),
                         pd.Period('2000-01-05', freq='D')],
                        index=[2, 4])
        tm.assert_series_equal(result, exp)
        assert result.dtype == 'object'

    def test_isna(self):
        # GH 13737
        s = Series([pd.Period('2011-01', freq='M'),
                    pd.Period('NaT', freq='M')])
        tm.assert_series_equal(s.isna(), Series([False, True]))
        tm.assert_series_equal(s.notna(), Series([True, False]))

    def test_fillna(self):
        # GH 13737
        s = Series([pd.Period('2011-01', freq='M'),
                    pd.Period('NaT', freq='M')])

        res = s.fillna(pd.Period('2012-01', freq='M'))
        exp = Series([pd.Period('2011-01', freq='M'),
                      pd.Period('2012-01', freq='M')])
        tm.assert_series_equal(res, exp)
        assert res.dtype == 'object'

        res = s.fillna('XXX')
        exp = Series([pd.Period('2011-01', freq='M'), 'XXX'])
        tm.assert_series_equal(res, exp)
        assert res.dtype == 'object'

    def test_dropna(self):
        # GH 13737
        s = Series([pd.Period('2011-01', freq='M'),
                    pd.Period('NaT', freq='M')])
        tm.assert_series_equal(s.dropna(),
                               Series([pd.Period('2011-01', freq='M')]))

    def test_between(self):
        left, right = self.series[[2, 7]]
        result = self.series.between(left, right)
        expected = (self.series >= left) & (self.series <= right)
        tm.assert_series_equal(result, expected)

    # ---------------------------------------------------------------------
    # NaT support

    """
    # ToDo: Enable when support period dtype
    def test_NaT_scalar(self):
        series = Series([0, 1000, 2000, iNaT], dtype='period[D]')

        val = series[3]
        assert isna(val)

        series[2] = val
        assert isna(series[2])

    def test_NaT_cast(self):
        result = Series([np.nan]).astype('period[D]')
        expected = Series([NaT])
        tm.assert_series_equal(result, expected)
    """

    def test_set_none_nan(self):
        # currently Period is stored as object dtype, not as NaT
        self.series[3] = None
        assert self.series[3] is None

        self.series[3:5] = None
        assert self.series[4] is None

        self.series[5] = np.nan
        assert np.isnan(self.series[5])

        self.series[5:7] = np.nan
        assert np.isnan(self.series[6])

    def test_intercept_astype_object(self):
        expected = self.series.astype('object')

        df = DataFrame({'a': self.series,
                        'b': np.random.randn(len(self.series))})

        result = df.values.squeeze()
        assert (result[:, 0] == expected.values).all()

        df = DataFrame({'a': self.series, 'b': ['foo'] * len(self.series)})

        result = df.values.squeeze()
        assert (result[:, 0] == expected.values).all()

    def test_add_series(self):
        rng = period_range('1/1/2000', '1/1/2010', freq='A')
        ts = Series(np.random.randn(len(rng)), index=rng)

        result = ts + ts[::2]
        expected = ts + ts
        expected[1::2] = np.nan
        tm.assert_series_equal(result, expected)

        result = ts + _permute(ts[::2])
        tm.assert_series_equal(result, expected)

        msg = "Input has different freq=D from PeriodIndex\\(freq=A-DEC\\)"
        with tm.assert_raises_regex(period.IncompatibleFrequency, msg):
            ts + ts.asfreq('D', how="end")

    def test_align_series(self, join_type):
        rng = period_range('1/1/2000', '1/1/2010', freq='A')
        ts = Series(np.random.randn(len(rng)), index=rng)

        ts.align(ts[::2], join=join_type)

    def test_truncate(self):
        # GH 17717
        idx1 = pd.PeriodIndex([
            pd.Period('2017-09-02'),
            pd.Period('2017-09-02'),
            pd.Period('2017-09-03')
        ])
        series1 = pd.Series([1, 2, 3], index=idx1)
        result1 = series1.truncate(after='2017-09-02')

        expected_idx1 = pd.PeriodIndex([
            pd.Period('2017-09-02'),
            pd.Period('2017-09-02')
        ])
        tm.assert_series_equal(result1, pd.Series([1, 2], index=expected_idx1))

        idx2 = pd.PeriodIndex([
            pd.Period('2017-09-03'),
            pd.Period('2017-09-02'),
            pd.Period('2017-09-03')
        ])
        series2 = pd.Series([1, 2, 3], index=idx2)
        result2 = series2.sort_index().truncate(after='2017-09-02')

        expected_idx2 = pd.PeriodIndex([
            pd.Period('2017-09-02')
        ])
        tm.assert_series_equal(result2, pd.Series([2], index=expected_idx2))
