# -*- coding: utf-8 -*-
from datetime import timedelta
from pandas.compat import zip
from pandas import compat
import re

import numpy as np

from pandas.core.dtypes.generic import ABCSeries
from pandas.core.dtypes.common import (
    is_period_arraylike,
    is_timedelta64_dtype,
    is_datetime64_dtype)

from pandas.tseries.offsets import DateOffset

from pandas._libs.tslib import Timedelta

import pandas._libs.tslibs.frequencies as libfreqs
from pandas._libs.tslibs.frequencies import (  # noqa, semi-public API
    get_freq, get_base_alias, get_to_timestamp_base, get_freq_code,
    FreqGroup,
    is_subperiod, is_superperiod)

from pandas._libs.tslibs.resolution import (Resolution,
                                            _FrequencyInferer,
                                            _TimedeltaFrequencyInferer)

from pytz import AmbiguousTimeError


RESO_NS = 0
RESO_US = 1
RESO_MS = 2
RESO_SEC = 3
RESO_MIN = 4
RESO_HR = 5
RESO_DAY = 6

# ---------------------------------------------------------------------
# Offset names ("time rules") and related functions

from pandas._libs.tslibs.offsets import _offset_to_period_map  # noqa:E402
from pandas.tseries.offsets import (Nano, Micro, Milli, Second,  # noqa
                                    Minute, Hour,
                                    Day, BDay, CDay, Week, MonthBegin,
                                    MonthEnd, BMonthBegin, BMonthEnd,
                                    QuarterBegin, QuarterEnd, BQuarterBegin,
                                    BQuarterEnd, YearBegin, YearEnd,
                                    BYearBegin, BYearEnd, prefix_mapping)
try:
    cday = CDay()
except NotImplementedError:
    cday = None

#: cache of previously seen offsets
_offset_map = {}


def get_period_alias(offset_str):
    """ alias to closest period strings BQ->Q etc"""
    return _offset_to_period_map.get(offset_str, None)


_name_to_offset_map = {'days': Day(1),
                       'hours': Hour(1),
                       'minutes': Minute(1),
                       'seconds': Second(1),
                       'milliseconds': Milli(1),
                       'microseconds': Micro(1),
                       'nanoseconds': Nano(1)}


def to_offset(freq):
    """
    Return DateOffset object from string or tuple representation
    or datetime.timedelta object

    Parameters
    ----------
    freq : str, tuple, datetime.timedelta, DateOffset or None

    Returns
    -------
    delta : DateOffset
        None if freq is None

    Raises
    ------
    ValueError
        If freq is an invalid frequency

    See Also
    --------
    pandas.DateOffset

    Examples
    --------
    >>> to_offset('5min')
    <5 * Minutes>

    >>> to_offset('1D1H')
    <25 * Hours>

    >>> to_offset(('W', 2))
    <2 * Weeks: weekday=6>

    >>> to_offset((2, 'B'))
    <2 * BusinessDays>

    >>> to_offset(datetime.timedelta(days=1))
    <Day>

    >>> to_offset(Hour())
    <Hour>
    """
    if freq is None:
        return None

    if isinstance(freq, DateOffset):
        return freq

    if isinstance(freq, tuple):
        name = freq[0]
        stride = freq[1]
        if isinstance(stride, compat.string_types):
            name, stride = stride, name
        name, _ = libfreqs._base_and_stride(name)
        delta = get_offset(name) * stride

    elif isinstance(freq, timedelta):
        delta = None
        freq = Timedelta(freq)
        try:
            for name in freq.components._fields:
                offset = _name_to_offset_map[name]
                stride = getattr(freq.components, name)
                if stride != 0:
                    offset = stride * offset
                    if delta is None:
                        delta = offset
                    else:
                        delta = delta + offset
        except Exception:
            raise ValueError(libfreqs._INVALID_FREQ_ERROR.format(freq))

    else:
        delta = None
        stride_sign = None
        try:
            splitted = re.split(libfreqs.opattern, freq)
            if splitted[-1] != '' and not splitted[-1].isspace():
                # the last element must be blank
                raise ValueError('last element must be blank')
            for sep, stride, name in zip(splitted[0::4], splitted[1::4],
                                         splitted[2::4]):
                if sep != '' and not sep.isspace():
                    raise ValueError('separator must be spaces')
                prefix = libfreqs._lite_rule_alias.get(name) or name
                if stride_sign is None:
                    stride_sign = -1 if stride.startswith('-') else 1
                if not stride:
                    stride = 1
                if prefix in Resolution._reso_str_bump_map.keys():
                    stride, name = Resolution.get_stride_from_decimal(
                        float(stride), prefix
                    )
                stride = int(stride)
                offset = get_offset(name)
                offset = offset * int(np.fabs(stride) * stride_sign)
                if delta is None:
                    delta = offset
                else:
                    delta = delta + offset
        except Exception:
            raise ValueError(libfreqs._INVALID_FREQ_ERROR.format(freq))

    if delta is None:
        raise ValueError(libfreqs._INVALID_FREQ_ERROR.format(freq))

    return delta


def get_offset(name):
    """
    Return DateOffset object associated with rule name

    Examples
    --------
    get_offset('EOM') --> BMonthEnd(1)
    """
    if name not in libfreqs._dont_uppercase:
        name = name.upper()
        name = libfreqs._lite_rule_alias.get(name, name)
        name = libfreqs._lite_rule_alias.get(name.lower(), name)
    else:
        name = libfreqs._lite_rule_alias.get(name, name)

    if name not in _offset_map:
        try:
            split = name.split('-')
            klass = prefix_mapping[split[0]]
            # handles case where there's no suffix (and will TypeError if too
            # many '-')
            offset = klass._from_name(*split[1:])
        except (ValueError, TypeError, KeyError):
            # bad prefix or suffix
            raise ValueError(libfreqs._INVALID_FREQ_ERROR.format(name))
        # cache
        _offset_map[name] = offset
    # do not return cache because it's mutable
    return _offset_map[name].copy()


getOffset = get_offset

# ---------------------------------------------------------------------
# Period codes


def infer_freq(index, warn=True):
    """
    Infer the most likely frequency given the input index. If the frequency is
    uncertain, a warning will be printed.

    Parameters
    ----------
    index : DatetimeIndex or TimedeltaIndex
      if passed a Series will use the values of the series (NOT THE INDEX)
    warn : boolean, default True

    Returns
    -------
    freq : string or None
        None if no discernible frequency
        TypeError if the index is not datetime-like
        ValueError if there are less than three values.
    """
    import pandas as pd

    if isinstance(index, ABCSeries):
        values = index._values
        if not (is_datetime64_dtype(values) or
                is_timedelta64_dtype(values) or
                values.dtype == object):
            raise TypeError("cannot infer freq from a non-convertible dtype "
                            "on a Series of {dtype}".format(dtype=index.dtype))
        index = values

    if is_period_arraylike(index):
        raise TypeError("PeriodIndex given. Check the `freq` attribute "
                        "instead of using infer_freq.")
    elif isinstance(index, pd.TimedeltaIndex):
        inferer = _TimedeltaFrequencyInferer(index, warn=warn)
        return inferer.get_freq()

    if isinstance(index, pd.Index) and not isinstance(index, pd.DatetimeIndex):
        if isinstance(index, (pd.Int64Index, pd.Float64Index)):
            raise TypeError("cannot infer freq from a non-convertible index "
                            "type {type}".format(type=type(index)))
        index = index.values

    if not isinstance(index, pd.DatetimeIndex):
        try:
            index = pd.DatetimeIndex(index)
        except AmbiguousTimeError:
            index = pd.DatetimeIndex(index.asi8)

    inferer = _FrequencyInferer(index, warn=warn)
    return inferer.get_freq()
