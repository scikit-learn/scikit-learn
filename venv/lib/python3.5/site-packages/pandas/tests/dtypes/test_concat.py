# -*- coding: utf-8 -*-

import pytest
import pandas.core.dtypes.concat as _concat
from pandas import (
    Index, DatetimeIndex, PeriodIndex, TimedeltaIndex, Series, Period)


@pytest.mark.parametrize('to_concat, expected', [
    # int/float/str
    ([['a'], [1, 2]], ['i', 'object']),
    ([[3, 4], [1, 2]], ['i']),
    ([[3, 4], [1, 2.1]], ['i', 'f']),

    # datetimelike
    ([DatetimeIndex(['2011-01-01']), DatetimeIndex(['2011-01-02'])],
     ['datetime']),
    ([TimedeltaIndex(['1 days']), TimedeltaIndex(['2 days'])],
     ['timedelta']),

    # datetimelike object
    ([DatetimeIndex(['2011-01-01']),
      DatetimeIndex(['2011-01-02'], tz='US/Eastern')],
     ['datetime', 'datetime64[ns, US/Eastern]']),
    ([DatetimeIndex(['2011-01-01'], tz='Asia/Tokyo'),
      DatetimeIndex(['2011-01-02'], tz='US/Eastern')],
     ['datetime64[ns, Asia/Tokyo]', 'datetime64[ns, US/Eastern]']),
    ([TimedeltaIndex(['1 days']), TimedeltaIndex(['2 hours'])],
     ['timedelta']),
    ([DatetimeIndex(['2011-01-01'], tz='Asia/Tokyo'),
      TimedeltaIndex(['1 days'])],
     ['datetime64[ns, Asia/Tokyo]', 'timedelta'])])
@pytest.mark.parametrize('klass', [Index, Series])
def test_get_dtype_kinds(klass, to_concat, expected):
    to_concat_klass = [klass(c) for c in to_concat]
    result = _concat.get_dtype_kinds(to_concat_klass)
    assert result == set(expected)


@pytest.mark.parametrize('to_concat, expected', [
    # because we don't have Period dtype (yet),
    # Series results in object dtype
    ([PeriodIndex(['2011-01'], freq='M'),
      PeriodIndex(['2011-01'], freq='M')], ['period[M]']),
    ([Series([Period('2011-01', freq='M')]),
      Series([Period('2011-02', freq='M')])], ['object']),
    ([PeriodIndex(['2011-01'], freq='M'),
      PeriodIndex(['2011-01'], freq='D')], ['period[M]', 'period[D]']),
    ([Series([Period('2011-01', freq='M')]),
      Series([Period('2011-02', freq='D')])], ['object'])])
def test_get_dtype_kinds_period(to_concat, expected):
    result = _concat.get_dtype_kinds(to_concat)
    assert result == set(expected)
