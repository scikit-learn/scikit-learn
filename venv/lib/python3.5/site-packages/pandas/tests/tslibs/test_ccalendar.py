# -*- coding: utf-8 -*-
from datetime import datetime

import numpy as np

from pandas._libs.tslibs import ccalendar


def test_get_day_of_year():
    assert ccalendar.get_day_of_year(2001, 3, 1) == 60
    assert ccalendar.get_day_of_year(2004, 3, 1) == 61
    assert ccalendar.get_day_of_year(1907, 12, 31) == 365
    assert ccalendar.get_day_of_year(2004, 12, 31) == 366

    dt = datetime.fromordinal(1 + np.random.randint(365 * 4000))
    result = ccalendar.get_day_of_year(dt.year, dt.month, dt.day)
    expected = (dt - dt.replace(month=1, day=1)).days + 1
    assert result == expected
