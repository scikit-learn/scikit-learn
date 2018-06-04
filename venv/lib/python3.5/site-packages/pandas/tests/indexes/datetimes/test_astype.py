import pytest

import pytz
import dateutil
import numpy as np

from datetime import datetime
from dateutil.tz import tzlocal

import pandas as pd
import pandas.util.testing as tm
from pandas import (DatetimeIndex, date_range, Series, NaT, Index, Timestamp,
                    Int64Index, Period)


class TestDatetimeIndex(object):

    def test_astype(self):
        # GH 13149, GH 13209
        idx = DatetimeIndex(['2016-05-16', 'NaT', NaT, np.NaN])

        result = idx.astype(object)
        expected = Index([Timestamp('2016-05-16')] + [NaT] * 3, dtype=object)
        tm.assert_index_equal(result, expected)

        result = idx.astype(int)
        expected = Int64Index([1463356800000000000] +
                              [-9223372036854775808] * 3, dtype=np.int64)
        tm.assert_index_equal(result, expected)

        rng = date_range('1/1/2000', periods=10)
        result = rng.astype('i8')
        tm.assert_index_equal(result, Index(rng.asi8))
        tm.assert_numpy_array_equal(result.values, rng.asi8)

    def test_astype_with_tz(self):

        # with tz
        rng = date_range('1/1/2000', periods=10, tz='US/Eastern')
        result = rng.astype('datetime64[ns]')
        expected = (date_range('1/1/2000', periods=10,
                               tz='US/Eastern')
                    .tz_convert('UTC').tz_localize(None))
        tm.assert_index_equal(result, expected)

        # BUG#10442 : testing astype(str) is correct for Series/DatetimeIndex
        result = pd.Series(pd.date_range('2012-01-01', periods=3)).astype(str)
        expected = pd.Series(
            ['2012-01-01', '2012-01-02', '2012-01-03'], dtype=object)
        tm.assert_series_equal(result, expected)

        result = Series(pd.date_range('2012-01-01', periods=3,
                                      tz='US/Eastern')).astype(str)
        expected = Series(['2012-01-01 00:00:00-05:00',
                           '2012-01-02 00:00:00-05:00',
                           '2012-01-03 00:00:00-05:00'],
                          dtype=object)
        tm.assert_series_equal(result, expected)

        # GH 18951: tz-aware to tz-aware
        idx = date_range('20170101', periods=4, tz='US/Pacific')
        result = idx.astype('datetime64[ns, US/Eastern]')
        expected = date_range('20170101 03:00:00', periods=4, tz='US/Eastern')
        tm.assert_index_equal(result, expected)

        # GH 18951: tz-naive to tz-aware
        idx = date_range('20170101', periods=4)
        result = idx.astype('datetime64[ns, US/Eastern]')
        expected = date_range('20170101', periods=4, tz='US/Eastern')
        tm.assert_index_equal(result, expected)

    def test_astype_str_compat(self):
        # GH 13149, GH 13209
        # verify that we are returning NaT as a string (and not unicode)

        idx = DatetimeIndex(['2016-05-16', 'NaT', NaT, np.NaN])
        result = idx.astype(str)
        expected = Index(['2016-05-16', 'NaT', 'NaT', 'NaT'], dtype=object)
        tm.assert_index_equal(result, expected)

    def test_astype_str(self):
        # test astype string - #10442
        result = date_range('2012-01-01', periods=4,
                            name='test_name').astype(str)
        expected = Index(['2012-01-01', '2012-01-02', '2012-01-03',
                          '2012-01-04'], name='test_name', dtype=object)
        tm.assert_index_equal(result, expected)

        # test astype string with tz and name
        result = date_range('2012-01-01', periods=3, name='test_name',
                            tz='US/Eastern').astype(str)
        expected = Index(['2012-01-01 00:00:00-05:00',
                          '2012-01-02 00:00:00-05:00',
                          '2012-01-03 00:00:00-05:00'],
                         name='test_name', dtype=object)
        tm.assert_index_equal(result, expected)

        # test astype string with freqH and name
        result = date_range('1/1/2011', periods=3, freq='H',
                            name='test_name').astype(str)
        expected = Index(['2011-01-01 00:00:00', '2011-01-01 01:00:00',
                          '2011-01-01 02:00:00'],
                         name='test_name', dtype=object)
        tm.assert_index_equal(result, expected)

        # test astype string with freqH and timezone
        result = date_range('3/6/2012 00:00', periods=2, freq='H',
                            tz='Europe/London', name='test_name').astype(str)
        expected = Index(['2012-03-06 00:00:00+00:00',
                          '2012-03-06 01:00:00+00:00'],
                         dtype=object, name='test_name')
        tm.assert_index_equal(result, expected)

    def test_astype_datetime64(self):
        # GH 13149, GH 13209
        idx = DatetimeIndex(['2016-05-16', 'NaT', NaT, np.NaN])

        result = idx.astype('datetime64[ns]')
        tm.assert_index_equal(result, idx)
        assert result is not idx

        result = idx.astype('datetime64[ns]', copy=False)
        tm.assert_index_equal(result, idx)
        assert result is idx

        idx_tz = DatetimeIndex(['2016-05-16', 'NaT', NaT, np.NaN], tz='EST')
        result = idx_tz.astype('datetime64[ns]')
        expected = DatetimeIndex(['2016-05-16 05:00:00', 'NaT', 'NaT', 'NaT'],
                                 dtype='datetime64[ns]')
        tm.assert_index_equal(result, expected)

    def test_astype_object(self):
        rng = date_range('1/1/2000', periods=20)

        casted = rng.astype('O')
        exp_values = list(rng)

        tm.assert_index_equal(casted, Index(exp_values, dtype=np.object_))
        assert casted.tolist() == exp_values

    @pytest.mark.parametrize('tz', [None, 'Asia/Tokyo'])
    def test_astype_object_tz(self, tz):
        idx = pd.date_range(start='2013-01-01', periods=4, freq='M',
                            name='idx', tz=tz)
        expected_list = [Timestamp('2013-01-31', tz=tz),
                         Timestamp('2013-02-28', tz=tz),
                         Timestamp('2013-03-31', tz=tz),
                         Timestamp('2013-04-30', tz=tz)]
        expected = pd.Index(expected_list, dtype=object, name='idx')
        result = idx.astype(object)
        tm.assert_index_equal(result, expected)
        assert idx.tolist() == expected_list

    def test_astype_object_with_nat(self):
        idx = DatetimeIndex([datetime(2013, 1, 1), datetime(2013, 1, 2),
                             pd.NaT, datetime(2013, 1, 4)], name='idx')
        expected_list = [Timestamp('2013-01-01'),
                         Timestamp('2013-01-02'), pd.NaT,
                         Timestamp('2013-01-04')]
        expected = pd.Index(expected_list, dtype=object, name='idx')
        result = idx.astype(object)
        tm.assert_index_equal(result, expected)
        assert idx.tolist() == expected_list

    @pytest.mark.parametrize('dtype', [
        float, 'timedelta64', 'timedelta64[ns]', 'datetime64',
        'datetime64[D]'])
    def test_astype_raises(self, dtype):
        # GH 13149, GH 13209
        idx = DatetimeIndex(['2016-05-16', 'NaT', NaT, np.NaN])
        msg = 'Cannot cast DatetimeIndex to dtype'
        with tm.assert_raises_regex(TypeError, msg):
            idx.astype(dtype)

    def test_index_convert_to_datetime_array(self):
        def _check_rng(rng):
            converted = rng.to_pydatetime()
            assert isinstance(converted, np.ndarray)
            for x, stamp in zip(converted, rng):
                assert isinstance(x, datetime)
                assert x == stamp.to_pydatetime()
                assert x.tzinfo == stamp.tzinfo

        rng = date_range('20090415', '20090519')
        rng_eastern = date_range('20090415', '20090519', tz='US/Eastern')
        rng_utc = date_range('20090415', '20090519', tz='utc')

        _check_rng(rng)
        _check_rng(rng_eastern)
        _check_rng(rng_utc)

    def test_index_convert_to_datetime_array_explicit_pytz(self):
        def _check_rng(rng):
            converted = rng.to_pydatetime()
            assert isinstance(converted, np.ndarray)
            for x, stamp in zip(converted, rng):
                assert isinstance(x, datetime)
                assert x == stamp.to_pydatetime()
                assert x.tzinfo == stamp.tzinfo

        rng = date_range('20090415', '20090519')
        rng_eastern = date_range('20090415', '20090519',
                                 tz=pytz.timezone('US/Eastern'))
        rng_utc = date_range('20090415', '20090519', tz=pytz.utc)

        _check_rng(rng)
        _check_rng(rng_eastern)
        _check_rng(rng_utc)

    def test_index_convert_to_datetime_array_dateutil(self):
        def _check_rng(rng):
            converted = rng.to_pydatetime()
            assert isinstance(converted, np.ndarray)
            for x, stamp in zip(converted, rng):
                assert isinstance(x, datetime)
                assert x == stamp.to_pydatetime()
                assert x.tzinfo == stamp.tzinfo

        rng = date_range('20090415', '20090519')
        rng_eastern = date_range('20090415', '20090519',
                                 tz='dateutil/US/Eastern')
        rng_utc = date_range('20090415', '20090519', tz=dateutil.tz.tzutc())

        _check_rng(rng)
        _check_rng(rng_eastern)
        _check_rng(rng_utc)


class TestToPeriod(object):

    def setup_method(self, method):
        data = [Timestamp('2007-01-01 10:11:12.123456Z'),
                Timestamp('2007-01-01 10:11:13.789123Z')]
        self.index = DatetimeIndex(data)

    def test_to_period_millisecond(self):
        index = self.index

        period = index.to_period(freq='L')
        assert 2 == len(period)
        assert period[0] == Period('2007-01-01 10:11:12.123Z', 'L')
        assert period[1] == Period('2007-01-01 10:11:13.789Z', 'L')

    def test_to_period_microsecond(self):
        index = self.index

        period = index.to_period(freq='U')
        assert 2 == len(period)
        assert period[0] == Period('2007-01-01 10:11:12.123456Z', 'U')
        assert period[1] == Period('2007-01-01 10:11:13.789123Z', 'U')

    def test_to_period_tz_pytz(self):
        from pytz import utc as UTC

        xp = date_range('1/1/2000', '4/1/2000').to_period()

        ts = date_range('1/1/2000', '4/1/2000', tz='US/Eastern')

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        assert result == expected
        tm.assert_index_equal(ts.to_period(), xp)

        ts = date_range('1/1/2000', '4/1/2000', tz=UTC)

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        assert result == expected
        tm.assert_index_equal(ts.to_period(), xp)

        ts = date_range('1/1/2000', '4/1/2000', tz=tzlocal())

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        assert result == expected
        tm.assert_index_equal(ts.to_period(), xp)

    def test_to_period_tz_explicit_pytz(self):
        xp = date_range('1/1/2000', '4/1/2000').to_period()

        ts = date_range('1/1/2000', '4/1/2000', tz=pytz.timezone('US/Eastern'))

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        assert result == expected
        tm.assert_index_equal(ts.to_period(), xp)

        ts = date_range('1/1/2000', '4/1/2000', tz=pytz.utc)

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        assert result == expected
        tm.assert_index_equal(ts.to_period(), xp)

        ts = date_range('1/1/2000', '4/1/2000', tz=tzlocal())

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        assert result == expected
        tm.assert_index_equal(ts.to_period(), xp)

    def test_to_period_tz_dateutil(self):
        xp = date_range('1/1/2000', '4/1/2000').to_period()

        ts = date_range('1/1/2000', '4/1/2000', tz='dateutil/US/Eastern')

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        assert result == expected
        tm.assert_index_equal(ts.to_period(), xp)

        ts = date_range('1/1/2000', '4/1/2000', tz=dateutil.tz.tzutc())

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        assert result == expected
        tm.assert_index_equal(ts.to_period(), xp)

        ts = date_range('1/1/2000', '4/1/2000', tz=tzlocal())

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        assert result == expected
        tm.assert_index_equal(ts.to_period(), xp)

    def test_to_period_nofreq(self):
        idx = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-04'])
        pytest.raises(ValueError, idx.to_period)

        idx = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-03'],
                            freq='infer')
        assert idx.freqstr == 'D'
        expected = pd.PeriodIndex(['2000-01-01', '2000-01-02',
                                   '2000-01-03'], freq='D')
        tm.assert_index_equal(idx.to_period(), expected)

        # GH 7606
        idx = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-03'])
        assert idx.freqstr is None
        tm.assert_index_equal(idx.to_period(), expected)
