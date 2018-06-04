# -*- coding: utf-8 -*-
from datetime import datetime

import pytest
import pytz
import dateutil.tz

from pandas._libs import tslib
from pandas._libs.tslibs import timezones
from pandas import Timestamp


@pytest.mark.parametrize('tz_name', list(pytz.common_timezones))
def test_cache_keys_are_distinct_for_pytz_vs_dateutil(tz_name):
    if tz_name == 'UTC':
        # skip utc as it's a special case in dateutil
        return
    tz_p = timezones.maybe_get_tz(tz_name)
    tz_d = timezones.maybe_get_tz('dateutil/' + tz_name)
    if tz_d is None:
        # skip timezones that dateutil doesn't know about.
        return
    assert timezones._p_tz_cache_key(tz_p) != timezones._p_tz_cache_key(tz_d)


def test_tzlocal():
    # GH#13583
    ts = Timestamp('2011-01-01', tz=dateutil.tz.tzlocal())
    assert ts.tz == dateutil.tz.tzlocal()
    assert "tz='tzlocal()')" in repr(ts)

    tz = timezones.maybe_get_tz('tzlocal()')
    assert tz == dateutil.tz.tzlocal()

    # get offset using normal datetime for test
    offset = dateutil.tz.tzlocal().utcoffset(datetime(2011, 1, 1))
    offset = offset.total_seconds() * 1000000000
    assert ts.value + offset == Timestamp('2011-01-01').value


@pytest.mark.parametrize('eastern, localize', [
    (pytz.timezone('US/Eastern'), lambda tz, x: tz.localize(x)),
    (dateutil.tz.gettz('US/Eastern'), lambda tz, x: x.replace(tzinfo=tz))])
def test_infer_tz(eastern, localize):
    utc = pytz.utc

    start_naive = datetime(2001, 1, 1)
    end_naive = datetime(2009, 1, 1)

    start = localize(eastern, start_naive)
    end = localize(eastern, end_naive)

    assert (timezones.infer_tzinfo(start, end) is
            tslib._localize_pydatetime(start_naive, eastern).tzinfo)
    assert (timezones.infer_tzinfo(start, None) is
            tslib._localize_pydatetime(start_naive, eastern).tzinfo)
    assert (timezones.infer_tzinfo(None, end) is
            tslib._localize_pydatetime(end_naive, eastern).tzinfo)

    start = utc.localize(start_naive)
    end = utc.localize(end_naive)
    assert timezones.infer_tzinfo(start, end) is utc

    end = tslib._localize_pydatetime(end_naive, eastern)
    with pytest.raises(Exception):
        timezones.infer_tzinfo(start, end)
    with pytest.raises(Exception):
        timezones.infer_tzinfo(end, start)
