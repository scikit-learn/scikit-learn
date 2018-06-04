from datetime import datetime, timedelta, time
import numpy as np
from collections import MutableMapping

from pandas._libs import tslib
from pandas._libs.tslibs.strptime import array_strptime
from pandas._libs.tslibs import parsing, conversion
from pandas._libs.tslibs.parsing import (  # noqa
    parse_time_string,
    DateParseError,
    _format_is_iso,
    _guess_datetime_format)

from pandas.core.dtypes.common import (
    _ensure_object,
    is_datetime64_ns_dtype,
    is_datetime64_dtype,
    is_datetime64tz_dtype,
    is_integer_dtype,
    is_integer,
    is_float,
    is_list_like,
    is_scalar,
    is_numeric_dtype)
from pandas.core.dtypes.generic import (
    ABCIndexClass, ABCSeries,
    ABCDataFrame)
from pandas.core.dtypes.missing import notna
from pandas.core import algorithms


def _guess_datetime_format_for_array(arr, **kwargs):
    # Try to guess the format based on the first non-NaN element
    non_nan_elements = notna(arr).nonzero()[0]
    if len(non_nan_elements):
        return _guess_datetime_format(arr[non_nan_elements[0]], **kwargs)


def _maybe_cache(arg, format, cache, tz, convert_listlike):
    """
    Create a cache of unique dates from an array of dates

    Parameters
    ----------
    arg : integer, float, string, datetime, list, tuple, 1-d array, Series
    format : string
        Strftime format to parse time
    cache : boolean
        True attempts to create a cache of converted values
    tz : string
        Timezone of the dates
    convert_listlike : function
        Conversion function to apply on dates

    Returns
    -------
    cache_array : Series
        Cache of converted, unique dates. Can be empty
    """
    from pandas import Series
    cache_array = Series()
    if cache:
        # Perform a quicker unique check
        from pandas import Index
        if not Index(arg).is_unique:
            unique_dates = algorithms.unique(arg)
            cache_dates = convert_listlike(unique_dates, True, format, tz=tz)
            cache_array = Series(cache_dates, index=unique_dates)
    return cache_array


def _convert_and_box_cache(arg, cache_array, box, errors, name=None):
    """
    Convert array of dates with a cache and box the result

    Parameters
    ----------
    arg : integer, float, string, datetime, list, tuple, 1-d array, Series
    cache_array : Series
        Cache of converted, unique dates
    box : boolean
        True boxes result as an Index-like, False returns an ndarray
    errors : string
        'ignore' plus box=True will convert result to Index
    name : string, default None
        Name for a DatetimeIndex

    Returns
    -------
    result : datetime of converted dates
        Returns:

        - Index-like if box=True
        - ndarray if box=False
    """
    from pandas import Series, DatetimeIndex, Index
    result = Series(arg).map(cache_array)
    if box:
        if errors == 'ignore':
            return Index(result)
        else:
            return DatetimeIndex(result, name=name)
    return result.values


def to_datetime(arg, errors='raise', dayfirst=False, yearfirst=False,
                utc=None, box=True, format=None, exact=True,
                unit=None, infer_datetime_format=False, origin='unix',
                cache=False):
    """
    Convert argument to datetime.

    Parameters
    ----------
    arg : integer, float, string, datetime, list, tuple, 1-d array, Series

        .. versionadded:: 0.18.1

           or DataFrame/dict-like

    errors : {'ignore', 'raise', 'coerce'}, default 'raise'

        - If 'raise', then invalid parsing will raise an exception
        - If 'coerce', then invalid parsing will be set as NaT
        - If 'ignore', then invalid parsing will return the input
    dayfirst : boolean, default False
        Specify a date parse order if `arg` is str or its list-likes.
        If True, parses dates with the day first, eg 10/11/12 is parsed as
        2012-11-10.
        Warning: dayfirst=True is not strict, but will prefer to parse
        with day first (this is a known bug, based on dateutil behavior).
    yearfirst : boolean, default False
        Specify a date parse order if `arg` is str or its list-likes.

        - If True parses dates with the year first, eg 10/11/12 is parsed as
          2010-11-12.
        - If both dayfirst and yearfirst are True, yearfirst is preceded (same
          as dateutil).

        Warning: yearfirst=True is not strict, but will prefer to parse
        with year first (this is a known bug, based on dateutil beahavior).

        .. versionadded:: 0.16.1

    utc : boolean, default None
        Return UTC DatetimeIndex if True (converting any tz-aware
        datetime.datetime objects as well).
    box : boolean, default True

        - If True returns a DatetimeIndex
        - If False returns ndarray of values.
    format : string, default None
        strftime to parse time, eg "%d/%m/%Y", note that "%f" will parse
        all the way up to nanoseconds.
    exact : boolean, True by default

        - If True, require an exact format match.
        - If False, allow the format to match anywhere in the target string.

    unit : string, default 'ns'
        unit of the arg (D,s,ms,us,ns) denote the unit, which is an
        integer or float number. This will be based off the origin.
        Example, with unit='ms' and origin='unix' (the default), this
        would calculate the number of milliseconds to the unix epoch start.
    infer_datetime_format : boolean, default False
        If True and no `format` is given, attempt to infer the format of the
        datetime strings, and if it can be inferred, switch to a faster
        method of parsing them. In some cases this can increase the parsing
        speed by ~5-10x.
    origin : scalar, default is 'unix'
        Define the reference date. The numeric values would be parsed as number
        of units (defined by `unit`) since this reference date.

        - If 'unix' (or POSIX) time; origin is set to 1970-01-01.
        - If 'julian', unit must be 'D', and origin is set to beginning of
          Julian Calendar. Julian day number 0 is assigned to the day starting
          at noon on January 1, 4713 BC.
        - If Timestamp convertible, origin is set to Timestamp identified by
          origin.

        .. versionadded:: 0.20.0
    cache : boolean, default False
        If True, use a cache of unique, converted dates to apply the datetime
        conversion. May produce sigificant speed-up when parsing duplicate date
        strings, especially ones with timezone offsets.

        .. versionadded:: 0.23.0

    Returns
    -------
    ret : datetime if parsing succeeded.
        Return type depends on input:

        - list-like: DatetimeIndex
        - Series: Series of datetime64 dtype
        - scalar: Timestamp

        In case when it is not possible to return designated types (e.g. when
        any element of input is before Timestamp.min or after Timestamp.max)
        return will have datetime.datetime type (or corresponding
        array/Series).

    Examples
    --------
    Assembling a datetime from multiple columns of a DataFrame. The keys can be
    common abbreviations like ['year', 'month', 'day', 'minute', 'second',
    'ms', 'us', 'ns']) or plurals of the same

    >>> df = pd.DataFrame({'year': [2015, 2016],
                           'month': [2, 3],
                           'day': [4, 5]})
    >>> pd.to_datetime(df)
    0   2015-02-04
    1   2016-03-05
    dtype: datetime64[ns]

    If a date does not meet the `timestamp limitations
    <http://pandas.pydata.org/pandas-docs/stable/timeseries.html
    #timeseries-timestamp-limits>`_, passing errors='ignore'
    will return the original input instead of raising any exception.

    Passing errors='coerce' will force an out-of-bounds date to NaT,
    in addition to forcing non-dates (or non-parseable dates) to NaT.

    >>> pd.to_datetime('13000101', format='%Y%m%d', errors='ignore')
    datetime.datetime(1300, 1, 1, 0, 0)
    >>> pd.to_datetime('13000101', format='%Y%m%d', errors='coerce')
    NaT

    Passing infer_datetime_format=True can often-times speedup a parsing
    if its not an ISO8601 format exactly, but in a regular format.

    >>> s = pd.Series(['3/11/2000', '3/12/2000', '3/13/2000']*1000)

    >>> s.head()
    0    3/11/2000
    1    3/12/2000
    2    3/13/2000
    3    3/11/2000
    4    3/12/2000
    dtype: object

    >>> %timeit pd.to_datetime(s,infer_datetime_format=True)
    100 loops, best of 3: 10.4 ms per loop

    >>> %timeit pd.to_datetime(s,infer_datetime_format=False)
    1 loop, best of 3: 471 ms per loop

    Using a unix epoch time

    >>> pd.to_datetime(1490195805, unit='s')
    Timestamp('2017-03-22 15:16:45')
    >>> pd.to_datetime(1490195805433502912, unit='ns')
    Timestamp('2017-03-22 15:16:45.433502912')

    .. warning:: For float arg, precision rounding might happen. To prevent
        unexpected behavior use a fixed-width exact type.

    Using a non-unix epoch origin

    >>> pd.to_datetime([1, 2, 3], unit='D',
                       origin=pd.Timestamp('1960-01-01'))
    0    1960-01-02
    1    1960-01-03
    2    1960-01-04

    See also
    --------
    pandas.DataFrame.astype : Cast argument to a specified dtype.
    pandas.to_timedelta : Convert argument to timedelta.
    """
    from pandas.core.indexes.datetimes import DatetimeIndex

    tz = 'utc' if utc else None

    def _convert_listlike(arg, box, format, name=None, tz=tz):

        if isinstance(arg, (list, tuple)):
            arg = np.array(arg, dtype='O')

        # these are shortcutable
        if is_datetime64tz_dtype(arg):
            if not isinstance(arg, DatetimeIndex):
                return DatetimeIndex(arg, tz=tz, name=name)
            if utc:
                arg = arg.tz_convert(None).tz_localize('UTC')
            return arg

        elif is_datetime64_ns_dtype(arg):
            if box and not isinstance(arg, DatetimeIndex):
                try:
                    return DatetimeIndex(arg, tz=tz, name=name)
                except ValueError:
                    pass

            return arg

        elif unit is not None:
            if format is not None:
                raise ValueError("cannot specify both format and unit")
            arg = getattr(arg, 'values', arg)
            result = tslib.array_with_unit_to_datetime(arg, unit,
                                                       errors=errors)
            if box:
                if errors == 'ignore':
                    from pandas import Index
                    return Index(result)

                return DatetimeIndex(result, tz=tz, name=name)
            return result
        elif getattr(arg, 'ndim', 1) > 1:
            raise TypeError('arg must be a string, datetime, list, tuple, '
                            '1-d array, or Series')

        arg = _ensure_object(arg)
        require_iso8601 = False

        if infer_datetime_format and format is None:
            format = _guess_datetime_format_for_array(arg, dayfirst=dayfirst)

        if format is not None:
            # There is a special fast-path for iso8601 formatted
            # datetime strings, so in those cases don't use the inferred
            # format because this path makes process slower in this
            # special case
            format_is_iso8601 = _format_is_iso(format)
            if format_is_iso8601:
                require_iso8601 = not infer_datetime_format
                format = None

        try:
            result = None

            if format is not None:
                # shortcut formatting here
                if format == '%Y%m%d':
                    try:
                        result = _attempt_YYYYMMDD(arg, errors=errors)
                    except:
                        raise ValueError("cannot convert the input to "
                                         "'%Y%m%d' date format")

                # fallback
                if result is None:
                    try:
                        result = array_strptime(arg, format, exact=exact,
                                                errors=errors)
                    except tslib.OutOfBoundsDatetime:
                        if errors == 'raise':
                            raise
                        result = arg
                    except ValueError:
                        # if format was inferred, try falling back
                        # to array_to_datetime - terminate here
                        # for specified formats
                        if not infer_datetime_format:
                            if errors == 'raise':
                                raise
                            result = arg

            if result is None and (format is None or infer_datetime_format):
                result = tslib.array_to_datetime(
                    arg,
                    errors=errors,
                    utc=utc,
                    dayfirst=dayfirst,
                    yearfirst=yearfirst,
                    require_iso8601=require_iso8601
                )

            if is_datetime64_dtype(result) and box:
                result = DatetimeIndex(result, tz=tz, name=name)
            return result

        except ValueError as e:
            try:
                values, tz = conversion.datetime_to_datetime64(arg)
                return DatetimeIndex._simple_new(values, name=name, tz=tz)
            except (ValueError, TypeError):
                raise e

    if arg is None:
        return None

    # handle origin
    if origin == 'julian':

        original = arg
        j0 = tslib.Timestamp(0).to_julian_date()
        if unit != 'D':
            raise ValueError("unit must be 'D' for origin='julian'")
        try:
            arg = arg - j0
        except:
            raise ValueError("incompatible 'arg' type for given "
                             "'origin'='julian'")

        # premptively check this for a nice range
        j_max = tslib.Timestamp.max.to_julian_date() - j0
        j_min = tslib.Timestamp.min.to_julian_date() - j0
        if np.any(arg > j_max) or np.any(arg < j_min):
            raise tslib.OutOfBoundsDatetime(
                "{original} is Out of Bounds for "
                "origin='julian'".format(original=original))

    elif origin not in ['unix', 'julian']:

        # arg must be a numeric
        original = arg
        if not ((is_scalar(arg) and (is_integer(arg) or is_float(arg))) or
                is_numeric_dtype(np.asarray(arg))):
            raise ValueError(
                "'{arg}' is not compatible with origin='{origin}'; "
                "it must be numeric with a unit specified ".format(
                    arg=arg,
                    origin=origin))

        # we are going to offset back to unix / epoch time
        try:
            offset = tslib.Timestamp(origin)
        except tslib.OutOfBoundsDatetime:
            raise tslib.OutOfBoundsDatetime(
                "origin {origin} is Out of Bounds".format(origin=origin))
        except ValueError:
            raise ValueError("origin {origin} cannot be converted "
                             "to a Timestamp".format(origin=origin))

        if offset.tz is not None:
            raise ValueError(
                "origin offset {} must be tz-naive".format(offset))
        offset -= tslib.Timestamp(0)

        # convert the offset to the unit of the arg
        # this should be lossless in terms of precision
        offset = offset // tslib.Timedelta(1, unit=unit)

        # scalars & ndarray-like can handle the addition
        if is_list_like(arg) and not isinstance(
                arg, (ABCSeries, ABCIndexClass, np.ndarray)):
            arg = np.asarray(arg)
        arg = arg + offset

    if isinstance(arg, tslib.Timestamp):
        result = arg
    elif isinstance(arg, ABCSeries):
        cache_array = _maybe_cache(arg, format, cache, tz, _convert_listlike)
        if not cache_array.empty:
            result = arg.map(cache_array)
        else:
            from pandas import Series
            values = _convert_listlike(arg._values, True, format)
            result = Series(values, index=arg.index, name=arg.name)
    elif isinstance(arg, (ABCDataFrame, MutableMapping)):
        result = _assemble_from_unit_mappings(arg, errors=errors)
    elif isinstance(arg, ABCIndexClass):
        cache_array = _maybe_cache(arg, format, cache, tz, _convert_listlike)
        if not cache_array.empty:
            result = _convert_and_box_cache(arg, cache_array, box, errors,
                                            name=arg.name)
        else:
            result = _convert_listlike(arg, box, format, name=arg.name)
    elif is_list_like(arg):
        cache_array = _maybe_cache(arg, format, cache, tz, _convert_listlike)
        if not cache_array.empty:
            result = _convert_and_box_cache(arg, cache_array, box, errors)
        else:
            result = _convert_listlike(arg, box, format)
    else:
        result = _convert_listlike(np.array([arg]), box, format)[0]

    return result


# mappings for assembling units
_unit_map = {'year': 'year',
             'years': 'year',
             'month': 'month',
             'months': 'month',
             'day': 'day',
             'days': 'day',
             'hour': 'h',
             'hours': 'h',
             'minute': 'm',
             'minutes': 'm',
             'second': 's',
             'seconds': 's',
             'ms': 'ms',
             'millisecond': 'ms',
             'milliseconds': 'ms',
             'us': 'us',
             'microsecond': 'us',
             'microseconds': 'us',
             'ns': 'ns',
             'nanosecond': 'ns',
             'nanoseconds': 'ns'
             }


def _assemble_from_unit_mappings(arg, errors):
    """
    assemble the unit specified fields from the arg (DataFrame)
    Return a Series for actual parsing

    Parameters
    ----------
    arg : DataFrame
    errors : {'ignore', 'raise', 'coerce'}, default 'raise'

        - If 'raise', then invalid parsing will raise an exception
        - If 'coerce', then invalid parsing will be set as NaT
        - If 'ignore', then invalid parsing will return the input

    Returns
    -------
    Series
    """
    from pandas import to_timedelta, to_numeric, DataFrame
    arg = DataFrame(arg)
    if not arg.columns.is_unique:
        raise ValueError("cannot assemble with duplicate keys")

    # replace passed unit with _unit_map
    def f(value):
        if value in _unit_map:
            return _unit_map[value]

        # m is case significant
        if value.lower() in _unit_map:
            return _unit_map[value.lower()]

        return value

    unit = {k: f(k) for k in arg.keys()}
    unit_rev = {v: k for k, v in unit.items()}

    # we require at least Ymd
    required = ['year', 'month', 'day']
    req = sorted(list(set(required) - set(unit_rev.keys())))
    if len(req):
        raise ValueError("to assemble mappings requires at least that "
                         "[year, month, day] be specified: [{required}] "
                         "is missing".format(required=','.join(req)))

    # keys we don't recognize
    excess = sorted(list(set(unit_rev.keys()) - set(_unit_map.values())))
    if len(excess):
        raise ValueError("extra keys have been passed "
                         "to the datetime assemblage: "
                         "[{excess}]".format(excess=','.join(excess)))

    def coerce(values):
        # we allow coercion to if errors allows
        values = to_numeric(values, errors=errors)

        # prevent overflow in case of int8 or int16
        if is_integer_dtype(values):
            values = values.astype('int64', copy=False)
        return values

    values = (coerce(arg[unit_rev['year']]) * 10000 +
              coerce(arg[unit_rev['month']]) * 100 +
              coerce(arg[unit_rev['day']]))
    try:
        values = to_datetime(values, format='%Y%m%d', errors=errors)
    except (TypeError, ValueError) as e:
        raise ValueError("cannot assemble the "
                         "datetimes: {error}".format(error=e))

    for u in ['h', 'm', 's', 'ms', 'us', 'ns']:
        value = unit_rev.get(u)
        if value is not None and value in arg:
            try:
                values += to_timedelta(coerce(arg[value]),
                                       unit=u,
                                       errors=errors)
            except (TypeError, ValueError) as e:
                raise ValueError("cannot assemble the datetimes [{value}]: "
                                 "{error}".format(value=value, error=e))

    return values


def _attempt_YYYYMMDD(arg, errors):
    """ try to parse the YYYYMMDD/%Y%m%d format, try to deal with NaT-like,
        arg is a passed in as an object dtype, but could really be ints/strings
        with nan-like/or floats (e.g. with nan)

    Parameters
    ----------
    arg : passed value
    errors : 'raise','ignore','coerce'
    """

    def calc(carg):
        # calculate the actual result
        carg = carg.astype(object)
        parsed = parsing.try_parse_year_month_day(carg / 10000,
                                                  carg / 100 % 100,
                                                  carg % 100)
        return tslib.array_to_datetime(parsed, errors=errors)

    def calc_with_mask(carg, mask):
        result = np.empty(carg.shape, dtype='M8[ns]')
        iresult = result.view('i8')
        iresult[~mask] = tslib.iNaT
        result[mask] = calc(carg[mask].astype(np.float64).astype(np.int64)). \
            astype('M8[ns]')
        return result

    # try intlike / strings that are ints
    try:
        return calc(arg.astype(np.int64))
    except:
        pass

    # a float with actual np.nan
    try:
        carg = arg.astype(np.float64)
        return calc_with_mask(carg, notna(carg))
    except:
        pass

    # string with NaN-like
    try:
        mask = ~algorithms.isin(arg, list(tslib.nat_strings))
        return calc_with_mask(arg, mask)
    except:
        pass

    return None


# Fixed time formats for time parsing
_time_formats = ["%H:%M", "%H%M", "%I:%M%p", "%I%M%p",
                 "%H:%M:%S", "%H%M%S", "%I:%M:%S%p", "%I%M%S%p"]


def _guess_time_format_for_array(arr):
    # Try to guess the format based on the first non-NaN element
    non_nan_elements = notna(arr).nonzero()[0]
    if len(non_nan_elements):
        element = arr[non_nan_elements[0]]
        for time_format in _time_formats:
            try:
                datetime.strptime(element, time_format)
                return time_format
            except ValueError:
                pass

    return None


def to_time(arg, format=None, infer_time_format=False, errors='raise'):
    """
    Parse time strings to time objects using fixed strptime formats ("%H:%M",
    "%H%M", "%I:%M%p", "%I%M%p", "%H:%M:%S", "%H%M%S", "%I:%M:%S%p",
    "%I%M%S%p")

    Use infer_time_format if all the strings are in the same format to speed
    up conversion.

    Parameters
    ----------
    arg : string in time format, datetime.time, list, tuple, 1-d array,  Series
    format : str, default None
        Format used to convert arg into a time object.  If None, fixed formats
        are used.
    infer_time_format: bool, default False
        Infer the time format based on the first non-NaN element.  If all
        strings are in the same format, this will speed up conversion.
    errors : {'ignore', 'raise', 'coerce'}, default 'raise'
        - If 'raise', then invalid parsing will raise an exception
        - If 'coerce', then invalid parsing will be set as None
        - If 'ignore', then invalid parsing will return the input

    Returns
    -------
    datetime.time
    """
    from pandas.core.series import Series

    def _convert_listlike(arg, format):

        if isinstance(arg, (list, tuple)):
            arg = np.array(arg, dtype='O')

        elif getattr(arg, 'ndim', 1) > 1:
            raise TypeError('arg must be a string, datetime, list, tuple, '
                            '1-d array, or Series')

        arg = _ensure_object(arg)

        if infer_time_format and format is None:
            format = _guess_time_format_for_array(arg)

        times = []
        if format is not None:
            for element in arg:
                try:
                    times.append(datetime.strptime(element, format).time())
                except (ValueError, TypeError):
                    if errors == 'raise':
                        msg = ("Cannot convert {element} to a time with given "
                               "format {format}").format(element=element,
                                                         format=format)
                        raise ValueError(msg)
                    elif errors == 'ignore':
                        return arg
                    else:
                        times.append(None)
        else:
            formats = _time_formats[:]
            format_found = False
            for element in arg:
                time_object = None
                for time_format in formats:
                    try:
                        time_object = datetime.strptime(element,
                                                        time_format).time()
                        if not format_found:
                            # Put the found format in front
                            fmt = formats.pop(formats.index(time_format))
                            formats.insert(0, fmt)
                            format_found = True
                        break
                    except (ValueError, TypeError):
                        continue

                if time_object is not None:
                    times.append(time_object)
                elif errors == 'raise':
                    raise ValueError("Cannot convert arg {arg} to "
                                     "a time".format(arg=arg))
                elif errors == 'ignore':
                    return arg
                else:
                    times.append(None)

        return times

    if arg is None:
        return arg
    elif isinstance(arg, time):
        return arg
    elif isinstance(arg, Series):
        values = _convert_listlike(arg._values, format)
        return Series(values, index=arg.index, name=arg.name)
    elif isinstance(arg, ABCIndexClass):
        return _convert_listlike(arg, format)
    elif is_list_like(arg):
        return _convert_listlike(arg, format)

    return _convert_listlike(np.array([arg]), format)[0]


def format(dt):
    """Returns date in YYYYMMDD format."""
    return dt.strftime('%Y%m%d')


OLE_TIME_ZERO = datetime(1899, 12, 30, 0, 0, 0)


def ole2datetime(oledt):
    """function for converting excel date to normal date format"""
    val = float(oledt)

    # Excel has a bug where it thinks the date 2/29/1900 exists
    # we just reject any date before 3/1/1900.
    if val < 61:
        msg = "Value is outside of acceptable range: {value}".format(value=val)
        raise ValueError(msg)

    return OLE_TIME_ZERO + timedelta(days=val)
