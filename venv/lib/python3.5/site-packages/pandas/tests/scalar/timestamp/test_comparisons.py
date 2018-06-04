# -*- coding: utf-8 -*-
from datetime import datetime
import operator

import pytest
import numpy as np

from dateutil.tz import tzutc
from pytz import utc

from pandas.compat import long, PY2
from pandas import Timestamp


class TestTimestampComparison(object):
    def test_comparison_object_array(self):
        # GH#15183
        ts = Timestamp('2011-01-03 00:00:00-0500', tz='US/Eastern')
        other = Timestamp('2011-01-01 00:00:00-0500', tz='US/Eastern')
        naive = Timestamp('2011-01-01 00:00:00')

        arr = np.array([other, ts], dtype=object)
        res = arr == ts
        expected = np.array([False, True], dtype=bool)
        assert (res == expected).all()

        # 2D case
        arr = np.array([[other, ts],
                        [ts, other]],
                       dtype=object)
        res = arr != ts
        expected = np.array([[True, False], [False, True]], dtype=bool)
        assert res.shape == expected.shape
        assert (res == expected).all()

        # tzaware mismatch
        arr = np.array([naive], dtype=object)
        with pytest.raises(TypeError):
            arr < ts

    def test_comparison(self):
        # 5-18-2012 00:00:00.000
        stamp = long(1337299200000000000)

        val = Timestamp(stamp)

        assert val == val
        assert not val != val
        assert not val < val
        assert val <= val
        assert not val > val
        assert val >= val

        other = datetime(2012, 5, 18)
        assert val == other
        assert not val != other
        assert not val < other
        assert val <= other
        assert not val > other
        assert val >= other

        other = Timestamp(stamp + 100)

        assert val != other
        assert val != other
        assert val < other
        assert val <= other
        assert other > val
        assert other >= val

    def test_compare_invalid(self):
        # GH 8058
        val = Timestamp('20130101 12:01:02')
        assert not val == 'foo'
        assert not val == 10.0
        assert not val == 1
        assert not val == long(1)
        assert not val == []
        assert not val == {'foo': 1}
        assert not val == np.float64(1)
        assert not val == np.int64(1)

        assert val != 'foo'
        assert val != 10.0
        assert val != 1
        assert val != long(1)
        assert val != []
        assert val != {'foo': 1}
        assert val != np.float64(1)
        assert val != np.int64(1)

    def test_cant_compare_tz_naive_w_aware(self):
        # see gh-1404
        a = Timestamp('3/12/2012')
        b = Timestamp('3/12/2012', tz='utc')

        pytest.raises(Exception, a.__eq__, b)
        pytest.raises(Exception, a.__ne__, b)
        pytest.raises(Exception, a.__lt__, b)
        pytest.raises(Exception, a.__gt__, b)
        pytest.raises(Exception, b.__eq__, a)
        pytest.raises(Exception, b.__ne__, a)
        pytest.raises(Exception, b.__lt__, a)
        pytest.raises(Exception, b.__gt__, a)

        if PY2:
            pytest.raises(Exception, a.__eq__, b.to_pydatetime())
            pytest.raises(Exception, a.to_pydatetime().__eq__, b)
        else:
            assert not a == b.to_pydatetime()
            assert not a.to_pydatetime() == b

    def test_cant_compare_tz_naive_w_aware_explicit_pytz(self):
        # see gh-1404
        a = Timestamp('3/12/2012')
        b = Timestamp('3/12/2012', tz=utc)

        pytest.raises(Exception, a.__eq__, b)
        pytest.raises(Exception, a.__ne__, b)
        pytest.raises(Exception, a.__lt__, b)
        pytest.raises(Exception, a.__gt__, b)
        pytest.raises(Exception, b.__eq__, a)
        pytest.raises(Exception, b.__ne__, a)
        pytest.raises(Exception, b.__lt__, a)
        pytest.raises(Exception, b.__gt__, a)

        if PY2:
            pytest.raises(Exception, a.__eq__, b.to_pydatetime())
            pytest.raises(Exception, a.to_pydatetime().__eq__, b)
        else:
            assert not a == b.to_pydatetime()
            assert not a.to_pydatetime() == b

    def test_cant_compare_tz_naive_w_aware_dateutil(self):
        # see gh-1404
        a = Timestamp('3/12/2012')
        b = Timestamp('3/12/2012', tz=tzutc())

        pytest.raises(Exception, a.__eq__, b)
        pytest.raises(Exception, a.__ne__, b)
        pytest.raises(Exception, a.__lt__, b)
        pytest.raises(Exception, a.__gt__, b)
        pytest.raises(Exception, b.__eq__, a)
        pytest.raises(Exception, b.__ne__, a)
        pytest.raises(Exception, b.__lt__, a)
        pytest.raises(Exception, b.__gt__, a)

        if PY2:
            pytest.raises(Exception, a.__eq__, b.to_pydatetime())
            pytest.raises(Exception, a.to_pydatetime().__eq__, b)
        else:
            assert not a == b.to_pydatetime()
            assert not a.to_pydatetime() == b

    def test_timestamp_compare_scalars(self):
        # case where ndim == 0
        lhs = np.datetime64(datetime(2013, 12, 6))
        rhs = Timestamp('now')
        nat = Timestamp('nat')

        ops = {'gt': 'lt',
               'lt': 'gt',
               'ge': 'le',
               'le': 'ge',
               'eq': 'eq',
               'ne': 'ne'}

        for left, right in ops.items():
            left_f = getattr(operator, left)
            right_f = getattr(operator, right)
            expected = left_f(lhs, rhs)

            result = right_f(rhs, lhs)
            assert result == expected

            expected = left_f(rhs, nat)
            result = right_f(nat, rhs)
            assert result == expected

    def test_timestamp_compare_with_early_datetime(self):
        # e.g. datetime.min
        stamp = Timestamp('2012-01-01')

        assert not stamp == datetime.min
        assert not stamp == datetime(1600, 1, 1)
        assert not stamp == datetime(2700, 1, 1)
        assert stamp != datetime.min
        assert stamp != datetime(1600, 1, 1)
        assert stamp != datetime(2700, 1, 1)
        assert stamp > datetime(1600, 1, 1)
        assert stamp >= datetime(1600, 1, 1)
        assert stamp < datetime(2700, 1, 1)
        assert stamp <= datetime(2700, 1, 1)
