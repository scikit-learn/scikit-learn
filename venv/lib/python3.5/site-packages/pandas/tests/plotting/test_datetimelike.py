""" Test cases for time series specific (freq conversion, etc) """

from datetime import datetime, timedelta, date, time
import pickle

import pytest
from pandas.compat import lrange, zip

import numpy as np
from pandas import Index, Series, DataFrame, NaT
from pandas.compat import PY3
from pandas.core.indexes.datetimes import date_range, bdate_range
from pandas.core.indexes.timedeltas import timedelta_range
from pandas.tseries.offsets import DateOffset
from pandas.core.indexes.period import period_range, Period, PeriodIndex
from pandas.core.resample import DatetimeIndex

from pandas.util.testing import assert_series_equal, ensure_clean
import pandas.util.testing as tm
import pandas.util._test_decorators as td

from pandas.tests.plotting.common import (TestPlotBase,
                                          _skip_if_no_scipy_gaussian_kde)


@td.skip_if_no_mpl
class TestTSPlot(TestPlotBase):

    def setup_method(self, method):
        TestPlotBase.setup_method(self, method)

        freq = ['S', 'T', 'H', 'D', 'W', 'M', 'Q', 'A']
        idx = [period_range('12/31/1999', freq=x, periods=100) for x in freq]
        self.period_ser = [Series(np.random.randn(len(x)), x) for x in idx]
        self.period_df = [DataFrame(np.random.randn(len(x), 3), index=x,
                                    columns=['A', 'B', 'C'])
                          for x in idx]

        freq = ['S', 'T', 'H', 'D', 'W', 'M', 'Q-DEC', 'A', '1B30Min']
        idx = [date_range('12/31/1999', freq=x, periods=100) for x in freq]
        self.datetime_ser = [Series(np.random.randn(len(x)), x) for x in idx]
        self.datetime_df = [DataFrame(np.random.randn(len(x), 3), index=x,
                                      columns=['A', 'B', 'C'])
                            for x in idx]

    def teardown_method(self, method):
        tm.close()

    @pytest.mark.slow
    def test_ts_plot_with_tz(self):
        # GH2877
        index = date_range('1/1/2011', periods=2, freq='H',
                           tz='Europe/Brussels')
        ts = Series([188.5, 328.25], index=index)
        _check_plot_works(ts.plot)

    def test_fontsize_set_correctly(self):
        # For issue #8765
        df = DataFrame(np.random.randn(10, 9), index=range(10))
        fig, ax = self.plt.subplots()
        df.plot(fontsize=2, ax=ax)
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            assert label.get_fontsize() == 2

    @pytest.mark.slow
    def test_frame_inferred(self):
        # inferred freq
        idx = date_range('1/1/1987', freq='MS', periods=100)
        idx = DatetimeIndex(idx.values, freq=None)

        df = DataFrame(np.random.randn(len(idx), 3), index=idx)
        _check_plot_works(df.plot)

        # axes freq
        idx = idx[0:40].union(idx[45:99])
        df2 = DataFrame(np.random.randn(len(idx), 3), index=idx)
        _check_plot_works(df2.plot)

        # N > 1
        idx = date_range('2008-1-1 00:15:00', freq='15T', periods=10)
        idx = DatetimeIndex(idx.values, freq=None)
        df = DataFrame(np.random.randn(len(idx), 3), index=idx)
        _check_plot_works(df.plot)

    def test_is_error_nozeroindex(self):
        # GH11858
        i = np.array([1, 2, 3])
        a = DataFrame(i, index=i)
        _check_plot_works(a.plot, xerr=a)
        _check_plot_works(a.plot, yerr=a)

    def test_nonnumeric_exclude(self):
        idx = date_range('1/1/1987', freq='A', periods=3)
        df = DataFrame({'A': ["x", "y", "z"], 'B': [1, 2, 3]}, idx)

        fig, ax = self.plt.subplots()
        df.plot(ax=ax)  # it works
        assert len(ax.get_lines()) == 1  # B was plotted
        self.plt.close(fig)

        pytest.raises(TypeError, df['A'].plot)

    def test_tsplot_deprecated(self):
        from pandas.tseries.plotting import tsplot

        _, ax = self.plt.subplots()
        ts = tm.makeTimeSeries()

        with tm.assert_produces_warning(FutureWarning):
            tsplot(ts, self.plt.Axes.plot, ax=ax)

    @pytest.mark.slow
    def test_tsplot(self):

        from pandas.tseries.plotting import tsplot

        _, ax = self.plt.subplots()
        ts = tm.makeTimeSeries()

        def f(*args, **kwds):
            with tm.assert_produces_warning(FutureWarning):
                return tsplot(s, self.plt.Axes.plot, *args, **kwds)

        for s in self.period_ser:
            _check_plot_works(f, s.index.freq, ax=ax, series=s)

        for s in self.datetime_ser:
            _check_plot_works(f, s.index.freq.rule_code, ax=ax, series=s)

        for s in self.period_ser:
            _check_plot_works(s.plot, ax=ax)

        for s in self.datetime_ser:
            _check_plot_works(s.plot, ax=ax)

        _, ax = self.plt.subplots()
        ts.plot(style='k', ax=ax)
        color = (0., 0., 0., 1) if self.mpl_ge_2_0_0 else (0., 0., 0.)
        assert color == ax.get_lines()[0].get_color()

    def test_both_style_and_color(self):

        ts = tm.makeTimeSeries()
        pytest.raises(ValueError, ts.plot, style='b-', color='#000099')

        s = ts.reset_index(drop=True)
        pytest.raises(ValueError, s.plot, style='b-', color='#000099')

    @pytest.mark.slow
    def test_high_freq(self):
        freaks = ['ms', 'us']
        for freq in freaks:
            _, ax = self.plt.subplots()
            rng = date_range('1/1/2012', periods=100000, freq=freq)
            ser = Series(np.random.randn(len(rng)), rng)
            _check_plot_works(ser.plot, ax=ax)

    def test_get_datevalue(self):
        from pandas.plotting._converter import get_datevalue
        assert get_datevalue(None, 'D') is None
        assert get_datevalue(1987, 'A') == 1987
        assert (get_datevalue(Period(1987, 'A'), 'M') ==
                Period('1987-12', 'M').ordinal)
        assert (get_datevalue('1/1/1987', 'D') ==
                Period('1987-1-1', 'D').ordinal)

    @pytest.mark.slow
    def test_ts_plot_format_coord(self):
        def check_format_of_first_point(ax, expected_string):
            first_line = ax.get_lines()[0]
            first_x = first_line.get_xdata()[0].ordinal
            first_y = first_line.get_ydata()[0]
            try:
                assert expected_string == ax.format_coord(first_x, first_y)
            except (ValueError):
                pytest.skip("skipping test because issue forming "
                            "test comparison GH7664")

        annual = Series(1, index=date_range('2014-01-01', periods=3,
                                            freq='A-DEC'))
        _, ax = self.plt.subplots()
        annual.plot(ax=ax)
        check_format_of_first_point(ax, 't = 2014  y = 1.000000')

        # note this is added to the annual plot already in existence, and
        # changes its freq field
        daily = Series(1, index=date_range('2014-01-01', periods=3, freq='D'))
        daily.plot(ax=ax)
        check_format_of_first_point(ax,
                                    't = 2014-01-01  y = 1.000000')
        tm.close()

        # tsplot
        from pandas.tseries.plotting import tsplot
        _, ax = self.plt.subplots()
        with tm.assert_produces_warning(FutureWarning):
            tsplot(annual, self.plt.Axes.plot, ax=ax)
        check_format_of_first_point(ax, 't = 2014  y = 1.000000')
        with tm.assert_produces_warning(FutureWarning):
            tsplot(daily, self.plt.Axes.plot, ax=ax)
        check_format_of_first_point(ax, 't = 2014-01-01  y = 1.000000')

    @pytest.mark.slow
    def test_line_plot_period_series(self):
        for s in self.period_ser:
            _check_plot_works(s.plot, s.index.freq)

    @pytest.mark.slow
    def test_line_plot_datetime_series(self):
        for s in self.datetime_ser:
            _check_plot_works(s.plot, s.index.freq.rule_code)

    @pytest.mark.slow
    def test_line_plot_period_frame(self):
        for df in self.period_df:
            _check_plot_works(df.plot, df.index.freq)

    @pytest.mark.slow
    def test_line_plot_datetime_frame(self):
        for df in self.datetime_df:
            freq = df.index.to_period(df.index.freq.rule_code).freq
            _check_plot_works(df.plot, freq)

    @pytest.mark.slow
    def test_line_plot_inferred_freq(self):
        for ser in self.datetime_ser:
            ser = Series(ser.values, Index(np.asarray(ser.index)))
            _check_plot_works(ser.plot, ser.index.inferred_freq)

            ser = ser[[0, 3, 5, 6]]
            _check_plot_works(ser.plot)

    def test_fake_inferred_business(self):
        _, ax = self.plt.subplots()
        rng = date_range('2001-1-1', '2001-1-10')
        ts = Series(lrange(len(rng)), rng)
        ts = ts[:3].append(ts[5:])
        ts.plot(ax=ax)
        assert not hasattr(ax, 'freq')

    @pytest.mark.slow
    def test_plot_offset_freq(self):
        ser = tm.makeTimeSeries()
        _check_plot_works(ser.plot)

        dr = date_range(ser.index[0], freq='BQS', periods=10)
        ser = Series(np.random.randn(len(dr)), dr)
        _check_plot_works(ser.plot)

    @pytest.mark.slow
    def test_plot_multiple_inferred_freq(self):
        dr = Index([datetime(2000, 1, 1), datetime(2000, 1, 6), datetime(
            2000, 1, 11)])
        ser = Series(np.random.randn(len(dr)), dr)
        _check_plot_works(ser.plot)

    @pytest.mark.slow
    def test_uhf(self):
        import pandas.plotting._converter as conv
        idx = date_range('2012-6-22 21:59:51.960928', freq='L', periods=500)
        df = DataFrame(np.random.randn(len(idx), 2), idx)

        _, ax = self.plt.subplots()
        df.plot(ax=ax)
        axis = ax.get_xaxis()

        tlocs = axis.get_ticklocs()
        tlabels = axis.get_ticklabels()
        for loc, label in zip(tlocs, tlabels):
            xp = conv._from_ordinal(loc).strftime('%H:%M:%S.%f')
            rs = str(label.get_text())
            if len(rs):
                assert xp == rs

    @pytest.mark.slow
    def test_irreg_hf(self):
        idx = date_range('2012-6-22 21:59:51', freq='S', periods=100)
        df = DataFrame(np.random.randn(len(idx), 2), idx)

        irreg = df.iloc[[0, 1, 3, 4]]
        _, ax = self.plt.subplots()
        irreg.plot(ax=ax)
        diffs = Series(ax.get_lines()[0].get_xydata()[:, 0]).diff()

        sec = 1. / 24 / 60 / 60
        assert (np.fabs(diffs[1:] - [sec, sec * 2, sec]) < 1e-8).all()

        _, ax = self.plt.subplots()
        df2 = df.copy()
        df2.index = df.index.astype(object)
        df2.plot(ax=ax)
        diffs = Series(ax.get_lines()[0].get_xydata()[:, 0]).diff()
        assert (np.fabs(diffs[1:] - sec) < 1e-8).all()

    def test_irregular_datetime64_repr_bug(self):
        ser = tm.makeTimeSeries()
        ser = ser[[0, 1, 2, 7]]

        _, ax = self.plt.subplots()

        ret = ser.plot(ax=ax)
        assert ret is not None

        for rs, xp in zip(ax.get_lines()[0].get_xdata(), ser.index):
            assert rs == xp

    def test_business_freq(self):
        bts = tm.makePeriodSeries()
        _, ax = self.plt.subplots()
        bts.plot(ax=ax)
        assert ax.get_lines()[0].get_xydata()[0, 0] == bts.index[0].ordinal
        idx = ax.get_lines()[0].get_xdata()
        assert PeriodIndex(data=idx).freqstr == 'B'

    @pytest.mark.slow
    def test_business_freq_convert(self):
        n = tm.N
        tm.N = 300
        bts = tm.makeTimeSeries().asfreq('BM')
        tm.N = n
        ts = bts.to_period('M')
        _, ax = self.plt.subplots()
        bts.plot(ax=ax)
        assert ax.get_lines()[0].get_xydata()[0, 0] == ts.index[0].ordinal
        idx = ax.get_lines()[0].get_xdata()
        assert PeriodIndex(data=idx).freqstr == 'M'

    def test_nonzero_base(self):
        # GH2571
        idx = (date_range('2012-12-20', periods=24, freq='H') + timedelta(
            minutes=30))
        df = DataFrame(np.arange(24), index=idx)
        _, ax = self.plt.subplots()
        df.plot(ax=ax)
        rs = ax.get_lines()[0].get_xdata()
        assert not Index(rs).is_normalized

    def test_dataframe(self):
        bts = DataFrame({'a': tm.makeTimeSeries()})
        _, ax = self.plt.subplots()
        bts.plot(ax=ax)
        idx = ax.get_lines()[0].get_xdata()
        tm.assert_index_equal(bts.index.to_period(), PeriodIndex(idx))

    @pytest.mark.slow
    def test_axis_limits(self):

        def _test(ax):
            xlim = ax.get_xlim()
            ax.set_xlim(xlim[0] - 5, xlim[1] + 10)
            ax.get_figure().canvas.draw()
            result = ax.get_xlim()
            assert result[0] == xlim[0] - 5
            assert result[1] == xlim[1] + 10

            # string
            expected = (Period('1/1/2000', ax.freq),
                        Period('4/1/2000', ax.freq))
            ax.set_xlim('1/1/2000', '4/1/2000')
            ax.get_figure().canvas.draw()
            result = ax.get_xlim()
            assert int(result[0]) == expected[0].ordinal
            assert int(result[1]) == expected[1].ordinal

            # datetime
            expected = (Period('1/1/2000', ax.freq),
                        Period('4/1/2000', ax.freq))
            ax.set_xlim(datetime(2000, 1, 1), datetime(2000, 4, 1))
            ax.get_figure().canvas.draw()
            result = ax.get_xlim()
            assert int(result[0]) == expected[0].ordinal
            assert int(result[1]) == expected[1].ordinal
            fig = ax.get_figure()
            self.plt.close(fig)

        ser = tm.makeTimeSeries()
        _, ax = self.plt.subplots()
        ser.plot(ax=ax)
        _test(ax)

        _, ax = self.plt.subplots()
        df = DataFrame({'a': ser, 'b': ser + 1})
        df.plot(ax=ax)
        _test(ax)

        df = DataFrame({'a': ser, 'b': ser + 1})
        axes = df.plot(subplots=True)

        for ax in axes:
            _test(ax)

    def test_get_finder(self):
        import pandas.plotting._converter as conv

        assert conv.get_finder('B') == conv._daily_finder
        assert conv.get_finder('D') == conv._daily_finder
        assert conv.get_finder('M') == conv._monthly_finder
        assert conv.get_finder('Q') == conv._quarterly_finder
        assert conv.get_finder('A') == conv._annual_finder
        assert conv.get_finder('W') == conv._daily_finder

    @pytest.mark.slow
    def test_finder_daily(self):
        day_lst = [10, 40, 252, 400, 950, 2750, 10000]

        if self.mpl_ge_2_0_0:
            xpl1 = [7565, 7564, 7553, 7546, 7518, 7428, 7066]
            xpl2 = [7566, 7564, 7554, 7546, 7519, 7429, 7066]
        else:
            xpl1 = xpl2 = [Period('1999-1-1', freq='B').ordinal] * len(day_lst)

        for i, n in enumerate(day_lst):
            xp = xpl1[i]
            rng = bdate_range('1999-1-1', periods=n)
            ser = Series(np.random.randn(len(rng)), rng)
            _, ax = self.plt.subplots()
            ser.plot(ax=ax)
            xaxis = ax.get_xaxis()
            rs = xaxis.get_majorticklocs()[0]
            assert xp == rs
            xp = xpl2[i]
            vmin, vmax = ax.get_xlim()
            ax.set_xlim(vmin + 0.9, vmax)
            rs = xaxis.get_majorticklocs()[0]
            assert xp == rs
            self.plt.close(ax.get_figure())

    @pytest.mark.slow
    def test_finder_quarterly(self):
        yrs = [3.5, 11]

        if self.mpl_ge_2_0_0:
            xpl1 = [68, 68]
            xpl2 = [72, 68]
        else:
            xpl1 = xpl2 = [Period('1988Q1').ordinal] * len(yrs)

        for i, n in enumerate(yrs):
            xp = xpl1[i]
            rng = period_range('1987Q2', periods=int(n * 4), freq='Q')
            ser = Series(np.random.randn(len(rng)), rng)
            _, ax = self.plt.subplots()
            ser.plot(ax=ax)
            xaxis = ax.get_xaxis()
            rs = xaxis.get_majorticklocs()[0]
            assert rs == xp
            xp = xpl2[i]
            (vmin, vmax) = ax.get_xlim()
            ax.set_xlim(vmin + 0.9, vmax)
            rs = xaxis.get_majorticklocs()[0]
            assert xp == rs
            self.plt.close(ax.get_figure())

    @pytest.mark.slow
    def test_finder_monthly(self):
        yrs = [1.15, 2.5, 4, 11]

        if self.mpl_ge_2_0_0:
            xpl1 = [216, 216, 204, 204]
            xpl2 = [216, 216, 216, 204]
        else:
            xpl1 = xpl2 = [Period('Jan 1988').ordinal] * len(yrs)

        for i, n in enumerate(yrs):
            xp = xpl1[i]
            rng = period_range('1987Q2', periods=int(n * 12), freq='M')
            ser = Series(np.random.randn(len(rng)), rng)
            _, ax = self.plt.subplots()
            ser.plot(ax=ax)
            xaxis = ax.get_xaxis()
            rs = xaxis.get_majorticklocs()[0]
            assert rs == xp
            xp = xpl2[i]
            vmin, vmax = ax.get_xlim()
            ax.set_xlim(vmin + 0.9, vmax)
            rs = xaxis.get_majorticklocs()[0]
            assert xp == rs
            self.plt.close(ax.get_figure())

    def test_finder_monthly_long(self):
        rng = period_range('1988Q1', periods=24 * 12, freq='M')
        ser = Series(np.random.randn(len(rng)), rng)
        _, ax = self.plt.subplots()
        ser.plot(ax=ax)
        xaxis = ax.get_xaxis()
        rs = xaxis.get_majorticklocs()[0]
        xp = Period('1989Q1', 'M').ordinal
        assert rs == xp

    @pytest.mark.slow
    def test_finder_annual(self):
        if self.mpl_ge_2_0_0:
            xp = [1986, 1986, 1990, 1990, 1995, 2020, 1970, 1970]
        else:
            xp = [1987, 1988, 1990, 1990, 1995, 2020, 2070, 2170]

        for i, nyears in enumerate([5, 10, 19, 49, 99, 199, 599, 1001]):
            rng = period_range('1987', periods=nyears, freq='A')
            ser = Series(np.random.randn(len(rng)), rng)
            _, ax = self.plt.subplots()
            ser.plot(ax=ax)
            xaxis = ax.get_xaxis()
            rs = xaxis.get_majorticklocs()[0]
            assert rs == Period(xp[i], freq='A').ordinal
            self.plt.close(ax.get_figure())

    @pytest.mark.slow
    def test_finder_minutely(self):
        nminutes = 50 * 24 * 60
        rng = date_range('1/1/1999', freq='Min', periods=nminutes)
        ser = Series(np.random.randn(len(rng)), rng)
        _, ax = self.plt.subplots()
        ser.plot(ax=ax)
        xaxis = ax.get_xaxis()
        rs = xaxis.get_majorticklocs()[0]
        if self.mpl_ge_2_0_0:
            xp = Period('1998-12-29 12:00', freq='Min').ordinal
        else:
            xp = Period('1/1/1999', freq='Min').ordinal
        assert rs == xp

    def test_finder_hourly(self):
        nhours = 23
        rng = date_range('1/1/1999', freq='H', periods=nhours)
        ser = Series(np.random.randn(len(rng)), rng)
        _, ax = self.plt.subplots()
        ser.plot(ax=ax)
        xaxis = ax.get_xaxis()
        rs = xaxis.get_majorticklocs()[0]
        if self.mpl_ge_2_0_0:
            xp = Period('1998-12-31 22:00', freq='H').ordinal
        else:
            xp = Period('1/1/1999', freq='H').ordinal
        assert rs == xp

    @td.skip_if_mpl_1_5
    @pytest.mark.slow
    def test_gaps(self):
        ts = tm.makeTimeSeries()
        ts[5:25] = np.nan
        _, ax = self.plt.subplots()
        ts.plot(ax=ax)
        lines = ax.get_lines()
        assert len(lines) == 1
        l = lines[0]
        data = l.get_xydata()
        assert isinstance(data, np.ma.core.MaskedArray)
        mask = data.mask
        assert mask[5:25, 1].all()
        self.plt.close(ax.get_figure())

        # irregular
        ts = tm.makeTimeSeries()
        ts = ts[[0, 1, 2, 5, 7, 9, 12, 15, 20]]
        ts[2:5] = np.nan
        _, ax = self.plt.subplots()
        ax = ts.plot(ax=ax)
        lines = ax.get_lines()
        assert len(lines) == 1
        l = lines[0]
        data = l.get_xydata()
        assert isinstance(data, np.ma.core.MaskedArray)
        mask = data.mask
        assert mask[2:5, 1].all()
        self.plt.close(ax.get_figure())

        # non-ts
        idx = [0, 1, 2, 5, 7, 9, 12, 15, 20]
        ser = Series(np.random.randn(len(idx)), idx)
        ser[2:5] = np.nan
        _, ax = self.plt.subplots()
        ser.plot(ax=ax)
        lines = ax.get_lines()
        assert len(lines) == 1
        l = lines[0]
        data = l.get_xydata()
        assert isinstance(data, np.ma.core.MaskedArray)
        mask = data.mask
        assert mask[2:5, 1].all()

    @td.skip_if_mpl_1_5
    @pytest.mark.slow
    def test_gap_upsample(self):
        low = tm.makeTimeSeries()
        low[5:25] = np.nan
        _, ax = self.plt.subplots()
        low.plot(ax=ax)

        idxh = date_range(low.index[0], low.index[-1], freq='12h')
        s = Series(np.random.randn(len(idxh)), idxh)
        s.plot(secondary_y=True)
        lines = ax.get_lines()
        assert len(lines) == 1
        assert len(ax.right_ax.get_lines()) == 1
        l = lines[0]
        data = l.get_xydata()

        assert isinstance(data, np.ma.core.MaskedArray)
        mask = data.mask
        assert mask[5:25, 1].all()

    @pytest.mark.slow
    def test_secondary_y(self):
        ser = Series(np.random.randn(10))
        ser2 = Series(np.random.randn(10))
        fig, _ = self.plt.subplots()
        ax = ser.plot(secondary_y=True)
        assert hasattr(ax, 'left_ax')
        assert not hasattr(ax, 'right_ax')
        axes = fig.get_axes()
        l = ax.get_lines()[0]
        xp = Series(l.get_ydata(), l.get_xdata())
        assert_series_equal(ser, xp)
        assert ax.get_yaxis().get_ticks_position() == 'right'
        assert not axes[0].get_yaxis().get_visible()
        self.plt.close(fig)

        _, ax2 = self.plt.subplots()
        ser2.plot(ax=ax2)
        assert (ax2.get_yaxis().get_ticks_position() ==
                self.default_tick_position)
        self.plt.close(ax2.get_figure())

        ax = ser2.plot()
        ax2 = ser.plot(secondary_y=True)
        assert ax.get_yaxis().get_visible()
        assert not hasattr(ax, 'left_ax')
        assert hasattr(ax, 'right_ax')
        assert hasattr(ax2, 'left_ax')
        assert not hasattr(ax2, 'right_ax')

    @pytest.mark.slow
    def test_secondary_y_ts(self):
        idx = date_range('1/1/2000', periods=10)
        ser = Series(np.random.randn(10), idx)
        ser2 = Series(np.random.randn(10), idx)
        fig, _ = self.plt.subplots()
        ax = ser.plot(secondary_y=True)
        assert hasattr(ax, 'left_ax')
        assert not hasattr(ax, 'right_ax')
        axes = fig.get_axes()
        l = ax.get_lines()[0]
        xp = Series(l.get_ydata(), l.get_xdata()).to_timestamp()
        assert_series_equal(ser, xp)
        assert ax.get_yaxis().get_ticks_position() == 'right'
        assert not axes[0].get_yaxis().get_visible()
        self.plt.close(fig)

        _, ax2 = self.plt.subplots()
        ser2.plot(ax=ax2)
        assert (ax2.get_yaxis().get_ticks_position() ==
                self.default_tick_position)
        self.plt.close(ax2.get_figure())

        ax = ser2.plot()
        ax2 = ser.plot(secondary_y=True)
        assert ax.get_yaxis().get_visible()

    @pytest.mark.slow
    @td.skip_if_no_scipy
    def test_secondary_kde(self):
        if not self.mpl_ge_1_5_0:
            pytest.skip("mpl is not supported")
        _skip_if_no_scipy_gaussian_kde()

        ser = Series(np.random.randn(10))
        fig, ax = self.plt.subplots()
        ax = ser.plot(secondary_y=True, kind='density', ax=ax)
        assert hasattr(ax, 'left_ax')
        assert not hasattr(ax, 'right_ax')
        axes = fig.get_axes()
        assert axes[1].get_yaxis().get_ticks_position() == 'right'

    @pytest.mark.slow
    def test_secondary_bar(self):
        ser = Series(np.random.randn(10))
        fig, ax = self.plt.subplots()
        ser.plot(secondary_y=True, kind='bar', ax=ax)
        axes = fig.get_axes()
        assert axes[1].get_yaxis().get_ticks_position() == 'right'

    @pytest.mark.slow
    def test_secondary_frame(self):
        df = DataFrame(np.random.randn(5, 3), columns=['a', 'b', 'c'])
        axes = df.plot(secondary_y=['a', 'c'], subplots=True)
        assert axes[0].get_yaxis().get_ticks_position() == 'right'
        assert (axes[1].get_yaxis().get_ticks_position() ==
                self.default_tick_position)
        assert axes[2].get_yaxis().get_ticks_position() == 'right'

    @pytest.mark.slow
    def test_secondary_bar_frame(self):
        df = DataFrame(np.random.randn(5, 3), columns=['a', 'b', 'c'])
        axes = df.plot(kind='bar', secondary_y=['a', 'c'], subplots=True)
        assert axes[0].get_yaxis().get_ticks_position() == 'right'
        assert (axes[1].get_yaxis().get_ticks_position() ==
                self.default_tick_position)
        assert axes[2].get_yaxis().get_ticks_position() == 'right'

    def test_mixed_freq_regular_first(self):
        # TODO
        s1 = tm.makeTimeSeries()
        s2 = s1[[0, 5, 10, 11, 12, 13, 14, 15]]

        # it works!
        _, ax = self.plt.subplots()
        s1.plot(ax=ax)

        ax2 = s2.plot(style='g', ax=ax)
        lines = ax2.get_lines()
        idx1 = PeriodIndex(lines[0].get_xdata())
        idx2 = PeriodIndex(lines[1].get_xdata())

        tm.assert_index_equal(idx1, s1.index.to_period('B'))
        tm.assert_index_equal(idx2, s2.index.to_period('B'))

        left, right = ax2.get_xlim()
        pidx = s1.index.to_period()
        assert left <= pidx[0].ordinal
        assert right >= pidx[-1].ordinal

    @pytest.mark.slow
    def test_mixed_freq_irregular_first(self):
        s1 = tm.makeTimeSeries()
        s2 = s1[[0, 5, 10, 11, 12, 13, 14, 15]]
        _, ax = self.plt.subplots()
        s2.plot(style='g', ax=ax)
        s1.plot(ax=ax)
        assert not hasattr(ax, 'freq')
        lines = ax.get_lines()
        x1 = lines[0].get_xdata()
        tm.assert_numpy_array_equal(x1, s2.index.astype(object).values)
        x2 = lines[1].get_xdata()
        tm.assert_numpy_array_equal(x2, s1.index.astype(object).values)

    def test_mixed_freq_regular_first_df(self):
        # GH 9852
        s1 = tm.makeTimeSeries().to_frame()
        s2 = s1.iloc[[0, 5, 10, 11, 12, 13, 14, 15], :]
        _, ax = self.plt.subplots()
        s1.plot(ax=ax)
        ax2 = s2.plot(style='g', ax=ax)
        lines = ax2.get_lines()
        idx1 = PeriodIndex(lines[0].get_xdata())
        idx2 = PeriodIndex(lines[1].get_xdata())
        assert idx1.equals(s1.index.to_period('B'))
        assert idx2.equals(s2.index.to_period('B'))
        left, right = ax2.get_xlim()
        pidx = s1.index.to_period()
        assert left <= pidx[0].ordinal
        assert right >= pidx[-1].ordinal

    @pytest.mark.slow
    def test_mixed_freq_irregular_first_df(self):
        # GH 9852
        s1 = tm.makeTimeSeries().to_frame()
        s2 = s1.iloc[[0, 5, 10, 11, 12, 13, 14, 15], :]
        _, ax = self.plt.subplots()
        s2.plot(style='g', ax=ax)
        s1.plot(ax=ax)
        assert not hasattr(ax, 'freq')
        lines = ax.get_lines()
        x1 = lines[0].get_xdata()
        tm.assert_numpy_array_equal(x1, s2.index.astype(object).values)
        x2 = lines[1].get_xdata()
        tm.assert_numpy_array_equal(x2, s1.index.astype(object).values)

    def test_mixed_freq_hf_first(self):
        idxh = date_range('1/1/1999', periods=365, freq='D')
        idxl = date_range('1/1/1999', periods=12, freq='M')
        high = Series(np.random.randn(len(idxh)), idxh)
        low = Series(np.random.randn(len(idxl)), idxl)
        _, ax = self.plt.subplots()
        high.plot(ax=ax)
        low.plot(ax=ax)
        for l in ax.get_lines():
            assert PeriodIndex(data=l.get_xdata()).freq == 'D'

    @pytest.mark.slow
    def test_mixed_freq_alignment(self):
        ts_ind = date_range('2012-01-01 13:00', '2012-01-02', freq='H')
        ts_data = np.random.randn(12)

        ts = Series(ts_data, index=ts_ind)
        ts2 = ts.asfreq('T').interpolate()

        _, ax = self.plt.subplots()
        ax = ts.plot(ax=ax)
        ts2.plot(style='r', ax=ax)

        assert ax.lines[0].get_xdata()[0] == ax.lines[1].get_xdata()[0]

    @pytest.mark.slow
    def test_mixed_freq_lf_first(self):

        idxh = date_range('1/1/1999', periods=365, freq='D')
        idxl = date_range('1/1/1999', periods=12, freq='M')
        high = Series(np.random.randn(len(idxh)), idxh)
        low = Series(np.random.randn(len(idxl)), idxl)
        _, ax = self.plt.subplots()
        low.plot(legend=True, ax=ax)
        high.plot(legend=True, ax=ax)
        for l in ax.get_lines():
            assert PeriodIndex(data=l.get_xdata()).freq == 'D'
        leg = ax.get_legend()
        assert len(leg.texts) == 2
        self.plt.close(ax.get_figure())

        idxh = date_range('1/1/1999', periods=240, freq='T')
        idxl = date_range('1/1/1999', periods=4, freq='H')
        high = Series(np.random.randn(len(idxh)), idxh)
        low = Series(np.random.randn(len(idxl)), idxl)
        _, ax = self.plt.subplots()
        low.plot(ax=ax)
        high.plot(ax=ax)
        for l in ax.get_lines():
            assert PeriodIndex(data=l.get_xdata()).freq == 'T'

    def test_mixed_freq_irreg_period(self):
        ts = tm.makeTimeSeries()
        irreg = ts[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 16, 17, 18, 29]]
        rng = period_range('1/3/2000', periods=30, freq='B')
        ps = Series(np.random.randn(len(rng)), rng)
        _, ax = self.plt.subplots()
        irreg.plot(ax=ax)
        ps.plot(ax=ax)

    def test_mixed_freq_shared_ax(self):

        # GH13341, using sharex=True
        idx1 = date_range('2015-01-01', periods=3, freq='M')
        idx2 = idx1[:1].union(idx1[2:])
        s1 = Series(range(len(idx1)), idx1)
        s2 = Series(range(len(idx2)), idx2)

        fig, (ax1, ax2) = self.plt.subplots(nrows=2, sharex=True)
        s1.plot(ax=ax1)
        s2.plot(ax=ax2)

        assert ax1.freq == 'M'
        assert ax2.freq == 'M'
        assert (ax1.lines[0].get_xydata()[0, 0] ==
                ax2.lines[0].get_xydata()[0, 0])

        # using twinx
        fig, ax1 = self.plt.subplots()
        ax2 = ax1.twinx()
        s1.plot(ax=ax1)
        s2.plot(ax=ax2)

        assert (ax1.lines[0].get_xydata()[0, 0] ==
                ax2.lines[0].get_xydata()[0, 0])

        # TODO (GH14330, GH14322)
        # plotting the irregular first does not yet work
        # fig, ax1 = plt.subplots()
        # ax2 = ax1.twinx()
        # s2.plot(ax=ax1)
        # s1.plot(ax=ax2)
        # assert (ax1.lines[0].get_xydata()[0, 0] ==
        #         ax2.lines[0].get_xydata()[0, 0])

    def test_nat_handling(self):

        _, ax = self.plt.subplots()

        dti = DatetimeIndex(['2015-01-01', NaT, '2015-01-03'])
        s = Series(range(len(dti)), dti)
        s.plot(ax=ax)
        xdata = ax.get_lines()[0].get_xdata()
        # plot x data is bounded by index values
        assert s.index.min() <= Series(xdata).min()
        assert Series(xdata).max() <= s.index.max()

    @pytest.mark.slow
    def test_to_weekly_resampling(self):
        idxh = date_range('1/1/1999', periods=52, freq='W')
        idxl = date_range('1/1/1999', periods=12, freq='M')
        high = Series(np.random.randn(len(idxh)), idxh)
        low = Series(np.random.randn(len(idxl)), idxl)
        _, ax = self.plt.subplots()
        high.plot(ax=ax)
        low.plot(ax=ax)
        for l in ax.get_lines():
            assert PeriodIndex(data=l.get_xdata()).freq == idxh.freq

        _, ax = self.plt.subplots()
        from pandas.tseries.plotting import tsplot
        with tm.assert_produces_warning(FutureWarning):
            tsplot(high, self.plt.Axes.plot, ax=ax)
        with tm.assert_produces_warning(FutureWarning):
            lines = tsplot(low, self.plt.Axes.plot, ax=ax)
        for l in lines:
            assert PeriodIndex(data=l.get_xdata()).freq == idxh.freq

    @pytest.mark.slow
    def test_from_weekly_resampling(self):
        idxh = date_range('1/1/1999', periods=52, freq='W')
        idxl = date_range('1/1/1999', periods=12, freq='M')
        high = Series(np.random.randn(len(idxh)), idxh)
        low = Series(np.random.randn(len(idxl)), idxl)
        _, ax = self.plt.subplots()
        low.plot(ax=ax)
        high.plot(ax=ax)

        expected_h = idxh.to_period().asi8.astype(np.float64)
        expected_l = np.array([1514, 1519, 1523, 1527, 1531, 1536, 1540, 1544,
                               1549, 1553, 1558, 1562], dtype=np.float64)
        for l in ax.get_lines():
            assert PeriodIndex(data=l.get_xdata()).freq == idxh.freq
            xdata = l.get_xdata(orig=False)
            if len(xdata) == 12:  # idxl lines
                tm.assert_numpy_array_equal(xdata, expected_l)
            else:
                tm.assert_numpy_array_equal(xdata, expected_h)
        tm.close()

        _, ax = self.plt.subplots()
        from pandas.tseries.plotting import tsplot
        with tm.assert_produces_warning(FutureWarning):
            tsplot(low, self.plt.Axes.plot, ax=ax)
        with tm.assert_produces_warning(FutureWarning):
            lines = tsplot(high, self.plt.Axes.plot, ax=ax)
        for l in lines:
            assert PeriodIndex(data=l.get_xdata()).freq == idxh.freq
            xdata = l.get_xdata(orig=False)
            if len(xdata) == 12:  # idxl lines
                tm.assert_numpy_array_equal(xdata, expected_l)
            else:
                tm.assert_numpy_array_equal(xdata, expected_h)

    @pytest.mark.slow
    def test_from_resampling_area_line_mixed(self):
        idxh = date_range('1/1/1999', periods=52, freq='W')
        idxl = date_range('1/1/1999', periods=12, freq='M')
        high = DataFrame(np.random.rand(len(idxh), 3),
                         index=idxh, columns=[0, 1, 2])
        low = DataFrame(np.random.rand(len(idxl), 3),
                        index=idxl, columns=[0, 1, 2])

        # low to high
        for kind1, kind2 in [('line', 'area'), ('area', 'line')]:
            _, ax = self.plt.subplots()
            low.plot(kind=kind1, stacked=True, ax=ax)
            high.plot(kind=kind2, stacked=True, ax=ax)

            # check low dataframe result
            expected_x = np.array([1514, 1519, 1523, 1527, 1531, 1536, 1540,
                                   1544, 1549, 1553, 1558, 1562],
                                  dtype=np.float64)
            expected_y = np.zeros(len(expected_x), dtype=np.float64)
            for i in range(3):
                l = ax.lines[i]
                assert PeriodIndex(l.get_xdata()).freq == idxh.freq
                tm.assert_numpy_array_equal(l.get_xdata(orig=False),
                                            expected_x)
                # check stacked values are correct
                expected_y += low[i].values
                tm.assert_numpy_array_equal(l.get_ydata(orig=False),
                                            expected_y)

            # check high dataframe result
            expected_x = idxh.to_period().asi8.astype(np.float64)
            expected_y = np.zeros(len(expected_x), dtype=np.float64)
            for i in range(3):
                l = ax.lines[3 + i]
                assert PeriodIndex(data=l.get_xdata()).freq == idxh.freq
                tm.assert_numpy_array_equal(l.get_xdata(orig=False),
                                            expected_x)
                expected_y += high[i].values
                tm.assert_numpy_array_equal(l.get_ydata(orig=False),
                                            expected_y)

        # high to low
        for kind1, kind2 in [('line', 'area'), ('area', 'line')]:
            _, ax = self.plt.subplots()
            high.plot(kind=kind1, stacked=True, ax=ax)
            low.plot(kind=kind2, stacked=True, ax=ax)

            # check high dataframe result
            expected_x = idxh.to_period().asi8.astype(np.float64)
            expected_y = np.zeros(len(expected_x), dtype=np.float64)
            for i in range(3):
                l = ax.lines[i]
                assert PeriodIndex(data=l.get_xdata()).freq == idxh.freq
                tm.assert_numpy_array_equal(l.get_xdata(orig=False),
                                            expected_x)
                expected_y += high[i].values
                tm.assert_numpy_array_equal(l.get_ydata(orig=False),
                                            expected_y)

            # check low dataframe result
            expected_x = np.array([1514, 1519, 1523, 1527, 1531, 1536, 1540,
                                   1544, 1549, 1553, 1558, 1562],
                                  dtype=np.float64)
            expected_y = np.zeros(len(expected_x), dtype=np.float64)
            for i in range(3):
                l = ax.lines[3 + i]
                assert PeriodIndex(data=l.get_xdata()).freq == idxh.freq
                tm.assert_numpy_array_equal(l.get_xdata(orig=False),
                                            expected_x)
                expected_y += low[i].values
                tm.assert_numpy_array_equal(l.get_ydata(orig=False),
                                            expected_y)

    @pytest.mark.slow
    def test_mixed_freq_second_millisecond(self):
        # GH 7772, GH 7760
        idxh = date_range('2014-07-01 09:00', freq='S', periods=50)
        idxl = date_range('2014-07-01 09:00', freq='100L', periods=500)
        high = Series(np.random.randn(len(idxh)), idxh)
        low = Series(np.random.randn(len(idxl)), idxl)
        # high to low
        _, ax = self.plt.subplots()
        high.plot(ax=ax)
        low.plot(ax=ax)
        assert len(ax.get_lines()) == 2
        for l in ax.get_lines():
            assert PeriodIndex(data=l.get_xdata()).freq == 'L'
        tm.close()

        # low to high
        _, ax = self.plt.subplots()
        low.plot(ax=ax)
        high.plot(ax=ax)
        assert len(ax.get_lines()) == 2
        for l in ax.get_lines():
            assert PeriodIndex(data=l.get_xdata()).freq == 'L'

    @pytest.mark.slow
    def test_irreg_dtypes(self):
        # date
        idx = [date(2000, 1, 1), date(2000, 1, 5), date(2000, 1, 20)]
        df = DataFrame(np.random.randn(len(idx), 3), Index(idx, dtype=object))
        _check_plot_works(df.plot)

        # np.datetime64
        idx = date_range('1/1/2000', periods=10)
        idx = idx[[0, 2, 5, 9]].astype(object)
        df = DataFrame(np.random.randn(len(idx), 3), idx)
        _, ax = self.plt.subplots()
        _check_plot_works(df.plot, ax=ax)

    @pytest.mark.xfail(not PY3, reason="failing on mpl 1.4.3 on PY2")
    @pytest.mark.slow
    def test_time(self):
        t = datetime(1, 1, 1, 3, 30, 0)
        deltas = np.random.randint(1, 20, 3).cumsum()
        ts = np.array([(t + timedelta(minutes=int(x))).time() for x in deltas])
        df = DataFrame({'a': np.random.randn(len(ts)),
                        'b': np.random.randn(len(ts))},
                       index=ts)
        fig, ax = self.plt.subplots()
        df.plot(ax=ax)

        # verify tick labels
        fig.canvas.draw()
        ticks = ax.get_xticks()
        labels = ax.get_xticklabels()
        for t, l in zip(ticks, labels):
            m, s = divmod(int(t), 60)
            h, m = divmod(m, 60)
            rs = l.get_text()
            if len(rs) > 0:
                if s != 0:
                    xp = time(h, m, s).strftime('%H:%M:%S')
                else:
                    xp = time(h, m, s).strftime('%H:%M')
                assert xp == rs

        # change xlim
        ax.set_xlim('1:30', '5:00')

        # check tick labels again
        fig.canvas.draw()
        ticks = ax.get_xticks()
        labels = ax.get_xticklabels()
        for t, l in zip(ticks, labels):
            m, s = divmod(int(t), 60)
            h, m = divmod(m, 60)
            rs = l.get_text()
            if len(rs) > 0:
                if s != 0:
                    xp = time(h, m, s).strftime('%H:%M:%S')
                else:
                    xp = time(h, m, s).strftime('%H:%M')
                assert xp == rs

    @pytest.mark.slow
    def test_time_musec(self):
        t = datetime(1, 1, 1, 3, 30, 0)
        deltas = np.random.randint(1, 20, 3).cumsum()
        ts = np.array([(t + timedelta(microseconds=int(x))).time()
                       for x in deltas])
        df = DataFrame({'a': np.random.randn(len(ts)),
                        'b': np.random.randn(len(ts))},
                       index=ts)
        fig, ax = self.plt.subplots()
        ax = df.plot(ax=ax)

        # verify tick labels
        fig.canvas.draw()
        ticks = ax.get_xticks()
        labels = ax.get_xticklabels()
        for t, l in zip(ticks, labels):
            m, s = divmod(int(t), 60)

            us = int(round((t - int(t)) * 1e6))

            h, m = divmod(m, 60)
            rs = l.get_text()
            if len(rs) > 0:
                if (us % 1000) != 0:
                    xp = time(h, m, s, us).strftime('%H:%M:%S.%f')
                elif (us // 1000) != 0:
                    xp = time(h, m, s, us).strftime('%H:%M:%S.%f')[:-3]
                elif s != 0:
                    xp = time(h, m, s, us).strftime('%H:%M:%S')
                else:
                    xp = time(h, m, s, us).strftime('%H:%M')
                assert xp == rs

    @pytest.mark.slow
    def test_secondary_upsample(self):
        idxh = date_range('1/1/1999', periods=365, freq='D')
        idxl = date_range('1/1/1999', periods=12, freq='M')
        high = Series(np.random.randn(len(idxh)), idxh)
        low = Series(np.random.randn(len(idxl)), idxl)
        _, ax = self.plt.subplots()
        low.plot(ax=ax)
        ax = high.plot(secondary_y=True, ax=ax)
        for l in ax.get_lines():
            assert PeriodIndex(l.get_xdata()).freq == 'D'
        assert hasattr(ax, 'left_ax')
        assert not hasattr(ax, 'right_ax')
        for l in ax.left_ax.get_lines():
            assert PeriodIndex(l.get_xdata()).freq == 'D'

    @pytest.mark.slow
    def test_secondary_legend(self):
        fig = self.plt.figure()
        ax = fig.add_subplot(211)

        # ts
        df = tm.makeTimeDataFrame()
        df.plot(secondary_y=['A', 'B'], ax=ax)
        leg = ax.get_legend()
        assert len(leg.get_lines()) == 4
        assert leg.get_texts()[0].get_text() == 'A (right)'
        assert leg.get_texts()[1].get_text() == 'B (right)'
        assert leg.get_texts()[2].get_text() == 'C'
        assert leg.get_texts()[3].get_text() == 'D'
        assert ax.right_ax.get_legend() is None
        colors = set()
        for line in leg.get_lines():
            colors.add(line.get_color())

        # TODO: color cycle problems
        assert len(colors) == 4
        self.plt.close(fig)

        fig = self.plt.figure()
        ax = fig.add_subplot(211)
        df.plot(secondary_y=['A', 'C'], mark_right=False, ax=ax)
        leg = ax.get_legend()
        assert len(leg.get_lines()) == 4
        assert leg.get_texts()[0].get_text() == 'A'
        assert leg.get_texts()[1].get_text() == 'B'
        assert leg.get_texts()[2].get_text() == 'C'
        assert leg.get_texts()[3].get_text() == 'D'
        self.plt.close(fig)

        fig, ax = self.plt.subplots()
        df.plot(kind='bar', secondary_y=['A'], ax=ax)
        leg = ax.get_legend()
        assert leg.get_texts()[0].get_text() == 'A (right)'
        assert leg.get_texts()[1].get_text() == 'B'
        self.plt.close(fig)

        fig, ax = self.plt.subplots()
        df.plot(kind='bar', secondary_y=['A'], mark_right=False, ax=ax)
        leg = ax.get_legend()
        assert leg.get_texts()[0].get_text() == 'A'
        assert leg.get_texts()[1].get_text() == 'B'
        self.plt.close(fig)

        fig = self.plt.figure()
        ax = fig.add_subplot(211)
        df = tm.makeTimeDataFrame()
        ax = df.plot(secondary_y=['C', 'D'], ax=ax)
        leg = ax.get_legend()
        assert len(leg.get_lines()) == 4
        assert ax.right_ax.get_legend() is None
        colors = set()
        for line in leg.get_lines():
            colors.add(line.get_color())

        # TODO: color cycle problems
        assert len(colors) == 4
        self.plt.close(fig)

        # non-ts
        df = tm.makeDataFrame()
        fig = self.plt.figure()
        ax = fig.add_subplot(211)
        ax = df.plot(secondary_y=['A', 'B'], ax=ax)
        leg = ax.get_legend()
        assert len(leg.get_lines()) == 4
        assert ax.right_ax.get_legend() is None
        colors = set()
        for line in leg.get_lines():
            colors.add(line.get_color())

        # TODO: color cycle problems
        assert len(colors) == 4
        self.plt.close()

        fig = self.plt.figure()
        ax = fig.add_subplot(211)
        ax = df.plot(secondary_y=['C', 'D'], ax=ax)
        leg = ax.get_legend()
        assert len(leg.get_lines()) == 4
        assert ax.right_ax.get_legend() is None
        colors = set()
        for line in leg.get_lines():
            colors.add(line.get_color())

        # TODO: color cycle problems
        assert len(colors) == 4

    def test_format_date_axis(self):
        rng = date_range('1/1/2012', periods=12, freq='M')
        df = DataFrame(np.random.randn(len(rng), 3), rng)
        _, ax = self.plt.subplots()
        ax = df.plot(ax=ax)
        xaxis = ax.get_xaxis()
        for l in xaxis.get_ticklabels():
            if len(l.get_text()) > 0:
                assert l.get_rotation() == 30

    @pytest.mark.slow
    def test_ax_plot(self):
        x = DatetimeIndex(start='2012-01-02', periods=10, freq='D')
        y = lrange(len(x))
        _, ax = self.plt.subplots()
        lines = ax.plot(x, y, label='Y')
        tm.assert_index_equal(DatetimeIndex(lines[0].get_xdata()), x)

    @pytest.mark.slow
    def test_mpl_nopandas(self):
        dates = [date(2008, 12, 31), date(2009, 1, 31)]
        values1 = np.arange(10.0, 11.0, 0.5)
        values2 = np.arange(11.0, 12.0, 0.5)

        kw = dict(fmt='-', lw=4)

        _, ax = self.plt.subplots()
        ax.plot_date([x.toordinal() for x in dates], values1, **kw)
        ax.plot_date([x.toordinal() for x in dates], values2, **kw)

        line1, line2 = ax.get_lines()

        exp = np.array([x.toordinal() for x in dates], dtype=np.float64)
        tm.assert_numpy_array_equal(line1.get_xydata()[:, 0], exp)
        exp = np.array([x.toordinal() for x in dates], dtype=np.float64)
        tm.assert_numpy_array_equal(line2.get_xydata()[:, 0], exp)

    @pytest.mark.slow
    def test_irregular_ts_shared_ax_xlim(self):
        # GH 2960
        ts = tm.makeTimeSeries()[:20]
        ts_irregular = ts[[1, 4, 5, 6, 8, 9, 10, 12, 13, 14, 15, 17, 18]]

        # plot the left section of the irregular series, then the right section
        _, ax = self.plt.subplots()
        ts_irregular[:5].plot(ax=ax)
        ts_irregular[5:].plot(ax=ax)

        # check that axis limits are correct
        left, right = ax.get_xlim()
        assert left <= ts_irregular.index.min().toordinal()
        assert right >= ts_irregular.index.max().toordinal()

    @pytest.mark.slow
    def test_secondary_y_non_ts_xlim(self):
        # GH 3490 - non-timeseries with secondary y
        index_1 = [1, 2, 3, 4]
        index_2 = [5, 6, 7, 8]
        s1 = Series(1, index=index_1)
        s2 = Series(2, index=index_2)

        _, ax = self.plt.subplots()
        s1.plot(ax=ax)
        left_before, right_before = ax.get_xlim()
        s2.plot(secondary_y=True, ax=ax)
        left_after, right_after = ax.get_xlim()

        assert left_before >= left_after
        assert right_before < right_after

    @pytest.mark.slow
    def test_secondary_y_regular_ts_xlim(self):
        # GH 3490 - regular-timeseries with secondary y
        index_1 = date_range(start='2000-01-01', periods=4, freq='D')
        index_2 = date_range(start='2000-01-05', periods=4, freq='D')
        s1 = Series(1, index=index_1)
        s2 = Series(2, index=index_2)

        _, ax = self.plt.subplots()
        s1.plot(ax=ax)
        left_before, right_before = ax.get_xlim()
        s2.plot(secondary_y=True, ax=ax)
        left_after, right_after = ax.get_xlim()

        assert left_before >= left_after
        assert right_before < right_after

    @pytest.mark.slow
    def test_secondary_y_mixed_freq_ts_xlim(self):
        # GH 3490 - mixed frequency timeseries with secondary y
        rng = date_range('2000-01-01', periods=10000, freq='min')
        ts = Series(1, index=rng)

        _, ax = self.plt.subplots()
        ts.plot(ax=ax)
        left_before, right_before = ax.get_xlim()
        ts.resample('D').mean().plot(secondary_y=True, ax=ax)
        left_after, right_after = ax.get_xlim()

        # a downsample should not have changed either limit
        assert left_before == left_after
        assert right_before == right_after

    @pytest.mark.slow
    def test_secondary_y_irregular_ts_xlim(self):
        # GH 3490 - irregular-timeseries with secondary y
        ts = tm.makeTimeSeries()[:20]
        ts_irregular = ts[[1, 4, 5, 6, 8, 9, 10, 12, 13, 14, 15, 17, 18]]

        _, ax = self.plt.subplots()
        ts_irregular[:5].plot(ax=ax)
        # plot higher-x values on secondary axis
        ts_irregular[5:].plot(secondary_y=True, ax=ax)
        # ensure secondary limits aren't overwritten by plot on primary
        ts_irregular[:5].plot(ax=ax)

        left, right = ax.get_xlim()
        assert left <= ts_irregular.index.min().toordinal()
        assert right >= ts_irregular.index.max().toordinal()

    def test_plot_outofbounds_datetime(self):
        # 2579 - checking this does not raise
        values = [date(1677, 1, 1), date(1677, 1, 2)]
        _, ax = self.plt.subplots()
        ax.plot(values)

        values = [datetime(1677, 1, 1, 12), datetime(1677, 1, 2, 12)]
        ax.plot(values)

    def test_format_timedelta_ticks_narrow(self):

        if self.mpl_ge_2_2_0:
            expected_labels = (['-1 days 23:59:59.999999998'] +
                               ['00:00:00.0000000{:0>2d}'.format(2 * i)
                                for i in range(6)])
        elif self.mpl_ge_2_0_0:
            expected_labels = [''] + [
                '00:00:00.00000000{:d}'.format(2 * i)
                for i in range(5)] + ['']
        else:
            expected_labels = [
                '00:00:00.00000000{:d}'.format(i)
                for i in range(10)]

        rng = timedelta_range('0', periods=10, freq='ns')
        df = DataFrame(np.random.randn(len(rng), 3), rng)
        fig, ax = self.plt.subplots()
        df.plot(fontsize=2, ax=ax)
        fig.canvas.draw()
        labels = ax.get_xticklabels()
        assert len(labels) == len(expected_labels)
        for l, l_expected in zip(labels, expected_labels):
            assert l.get_text() == l_expected

    def test_format_timedelta_ticks_wide(self):

        if self.mpl_ge_2_0_0:
            expected_labels = [
                '',
                '00:00:00',
                '1 days 03:46:40',
                '2 days 07:33:20',
                '3 days 11:20:00',
                '4 days 15:06:40',
                '5 days 18:53:20',
                '6 days 22:40:00',
                '8 days 02:26:40',
                '9 days 06:13:20',
                ''
            ]
            if self.mpl_ge_2_2_0:
                expected_labels[0] = '-2 days 20:13:20'
                expected_labels[-1] = '10 days 10:00:00'
        else:
            expected_labels = [
                '00:00:00',
                '1 days 03:46:40',
                '2 days 07:33:20',
                '3 days 11:20:00',
                '4 days 15:06:40',
                '5 days 18:53:20',
                '6 days 22:40:00',
                '8 days 02:26:40',
                ''
            ]

        rng = timedelta_range('0', periods=10, freq='1 d')
        df = DataFrame(np.random.randn(len(rng), 3), rng)
        fig, ax = self.plt.subplots()
        ax = df.plot(fontsize=2, ax=ax)
        fig.canvas.draw()
        labels = ax.get_xticklabels()
        assert len(labels) == len(expected_labels)
        for l, l_expected in zip(labels, expected_labels):
            assert l.get_text() == l_expected

    def test_timedelta_plot(self):
        # test issue #8711
        s = Series(range(5), timedelta_range('1day', periods=5))
        _, ax = self.plt.subplots()
        _check_plot_works(s.plot, ax=ax)

        # test long period
        index = timedelta_range('1 day 2 hr 30 min 10 s',
                                periods=10, freq='1 d')
        s = Series(np.random.randn(len(index)), index)
        _, ax = self.plt.subplots()
        _check_plot_works(s.plot, ax=ax)

        # test short period
        index = timedelta_range('1 day 2 hr 30 min 10 s',
                                periods=10, freq='1 ns')
        s = Series(np.random.randn(len(index)), index)
        _, ax = self.plt.subplots()
        _check_plot_works(s.plot, ax=ax)

    def test_hist(self):
        # https://github.com/matplotlib/matplotlib/issues/8459
        rng = date_range('1/1/2011', periods=10, freq='H')
        x = rng
        w1 = np.arange(0, 1, .1)
        w2 = np.arange(0, 1, .1)[::-1]
        _, ax = self.plt.subplots()
        ax.hist([x, x], weights=[w1, w2])

    @pytest.mark.slow
    def test_overlapping_datetime(self):
        # GB 6608
        s1 = Series([1, 2, 3], index=[datetime(1995, 12, 31),
                                      datetime(2000, 12, 31),
                                      datetime(2005, 12, 31)])
        s2 = Series([1, 2, 3], index=[datetime(1997, 12, 31),
                                      datetime(2003, 12, 31),
                                      datetime(2008, 12, 31)])

        # plot first series, then add the second series to those axes,
        # then try adding the first series again
        _, ax = self.plt.subplots()
        s1.plot(ax=ax)
        s2.plot(ax=ax)
        s1.plot(ax=ax)

    @pytest.mark.xfail(reason="GH9053 matplotlib does not use"
                              " ax.xaxis.converter")
    def test_add_matplotlib_datetime64(self):
        # GH9053 - ensure that a plot with PeriodConverter still understands
        # datetime64 data. This still fails because matplotlib overrides the
        # ax.xaxis.converter with a DatetimeConverter
        s = Series(np.random.randn(10),
                   index=date_range('1970-01-02', periods=10))
        ax = s.plot()
        ax.plot(s.index, s.values, color='g')
        l1, l2 = ax.lines
        tm.assert_numpy_array_equal(l1.get_xydata(), l2.get_xydata())


def _check_plot_works(f, freq=None, series=None, *args, **kwargs):
    import matplotlib.pyplot as plt

    fig = plt.gcf()

    try:
        plt.clf()
        ax = fig.add_subplot(211)
        orig_ax = kwargs.pop('ax', plt.gca())
        orig_axfreq = getattr(orig_ax, 'freq', None)

        ret = f(*args, **kwargs)
        assert ret is not None  # do something more intelligent

        ax = kwargs.pop('ax', plt.gca())
        if series is not None:
            dfreq = series.index.freq
            if isinstance(dfreq, DateOffset):
                dfreq = dfreq.rule_code
            if orig_axfreq is None:
                assert ax.freq == dfreq

        if freq is not None and orig_axfreq is None:
            assert ax.freq == freq

        ax = fig.add_subplot(212)
        try:
            kwargs['ax'] = ax
            ret = f(*args, **kwargs)
            assert ret is not None  # do something more intelligent
        except Exception:
            pass

        with ensure_clean(return_filelike=True) as path:
            plt.savefig(path)

        # GH18439
        # this is supported only in Python 3 pickle since
        # pickle in Python2 doesn't support instancemethod pickling
        if PY3:
            with ensure_clean(return_filelike=True) as path:
                pickle.dump(fig, path)
    finally:
        plt.close(fig)
