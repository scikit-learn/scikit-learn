#!/usr/env/bin python

"""
self-contained to write legacy storage (pickle/msgpack) files

To use this script. Create an environment where you want
generate pickles, say its for 0.18.1, with your pandas clone
in ~/pandas

. activate pandas_0.18.1
cd ~/

$ python pandas/pandas/tests/io/generate_legacy_storage_files.py \
    pandas/pandas/tests/io/data/legacy_pickle/0.18.1/ pickle

This script generates a storage file for the current arch, system,
and python version
  pandas version: 0.18.1
  output dir    : pandas/pandas/tests/io/data/legacy_pickle/0.18.1/
  storage format: pickle
created pickle file: 0.18.1_x86_64_darwin_3.5.2.pickle

The idea here is you are using the *current* version of the
generate_legacy_storage_files with an *older* version of pandas to
generate a pickle file. We will then check this file into a current
branch, and test using test_pickle.py. This will load the *older*
pickles and test versus the current data that is generated
(with master). These are then compared.

If we have cases where we changed the signature (e.g. we renamed
offset -> freq in Timestamp). Then we have to conditionally execute
in the generate_legacy_storage_files.py to make it
run under the older AND the newer version.

"""

from __future__ import print_function
from warnings import catch_warnings
from distutils.version import LooseVersion
from pandas import (Series, DataFrame, Panel,
                    SparseSeries, SparseDataFrame,
                    Index, MultiIndex, bdate_range, to_msgpack,
                    date_range, period_range, timedelta_range,
                    Timestamp, NaT, Categorical, Period)
from pandas.tseries.offsets import (
    DateOffset, Hour, Minute, Day,
    MonthBegin, MonthEnd, YearBegin,
    YearEnd, Week, WeekOfMonth, LastWeekOfMonth,
    BusinessDay, BusinessHour, CustomBusinessDay, FY5253,
    Easter,
    SemiMonthEnd, SemiMonthBegin,
    QuarterBegin, QuarterEnd)
from pandas.compat import u
import os
import sys
import numpy as np
import pandas
import platform as pl
from datetime import timedelta

_loose_version = LooseVersion(pandas.__version__)


def _create_sp_series():
    nan = np.nan

    # nan-based
    arr = np.arange(15, dtype=np.float64)
    arr[7:12] = nan
    arr[-1:] = nan

    bseries = SparseSeries(arr, kind='block')
    bseries.name = u'bseries'
    return bseries


def _create_sp_tsseries():
    nan = np.nan

    # nan-based
    arr = np.arange(15, dtype=np.float64)
    arr[7:12] = nan
    arr[-1:] = nan

    date_index = bdate_range('1/1/2011', periods=len(arr))
    bseries = SparseSeries(arr, index=date_index, kind='block')
    bseries.name = u'btsseries'
    return bseries


def _create_sp_frame():
    nan = np.nan

    data = {u'A': [nan, nan, nan, 0, 1, 2, 3, 4, 5, 6],
            u'B': [0, 1, 2, nan, nan, nan, 3, 4, 5, 6],
            u'C': np.arange(10).astype(np.int64),
            u'D': [0, 1, 2, 3, 4, 5, nan, nan, nan, nan]}

    dates = bdate_range('1/1/2011', periods=10)
    return SparseDataFrame(data, index=dates)


def create_data():
    """ create the pickle/msgpack data """

    data = {
        u'A': [0., 1., 2., 3., np.nan],
        u'B': [0, 1, 0, 1, 0],
        u'C': [u'foo1', u'foo2', u'foo3', u'foo4', u'foo5'],
        u'D': date_range('1/1/2009', periods=5),
        u'E': [0., 1, Timestamp('20100101'), u'foo', 2.]
    }

    scalars = dict(timestamp=Timestamp('20130101'),
                   period=Period('2012', 'M'))

    index = dict(int=Index(np.arange(10)),
                 date=date_range('20130101', periods=10),
                 period=period_range('2013-01-01', freq='M', periods=10),
                 float=Index(np.arange(10, dtype=np.float64)),
                 uint=Index(np.arange(10, dtype=np.uint64)),
                 timedelta=timedelta_range('00:00:00', freq='30T', periods=10))

    if _loose_version >= LooseVersion('0.18'):
        from pandas import RangeIndex
        index['range'] = RangeIndex(10)

    if _loose_version >= LooseVersion('0.21'):
        from pandas import interval_range
        index['interval'] = interval_range(0, periods=10)

    mi = dict(reg2=MultiIndex.from_tuples(
        tuple(zip(*[[u'bar', u'bar', u'baz', u'baz', u'foo',
                     u'foo', u'qux', u'qux'],
                    [u'one', u'two', u'one', u'two', u'one',
                     u'two', u'one', u'two']])),
        names=[u'first', u'second']))

    series = dict(float=Series(data[u'A']),
                  int=Series(data[u'B']),
                  mixed=Series(data[u'E']),
                  ts=Series(np.arange(10).astype(np.int64),
                            index=date_range('20130101', periods=10)),
                  mi=Series(np.arange(5).astype(np.float64),
                            index=MultiIndex.from_tuples(
                                tuple(zip(*[[1, 1, 2, 2, 2],
                                            [3, 4, 3, 4, 5]])),
                                names=[u'one', u'two'])),
                  dup=Series(np.arange(5).astype(np.float64),
                             index=[u'A', u'B', u'C', u'D', u'A']),
                  cat=Series(Categorical([u'foo', u'bar', u'baz'])),
                  dt=Series(date_range('20130101', periods=5)),
                  dt_tz=Series(date_range('20130101', periods=5,
                                          tz='US/Eastern')),
                  period=Series([Period('2000Q1')] * 5))

    mixed_dup_df = DataFrame(data)
    mixed_dup_df.columns = list(u"ABCDA")
    frame = dict(float=DataFrame({u'A': series[u'float'],
                                  u'B': series[u'float'] + 1}),
                 int=DataFrame({u'A': series[u'int'],
                                u'B': series[u'int'] + 1}),
                 mixed=DataFrame({k: data[k]
                                  for k in [u'A', u'B', u'C', u'D']}),
                 mi=DataFrame({u'A': np.arange(5).astype(np.float64),
                               u'B': np.arange(5).astype(np.int64)},
                              index=MultiIndex.from_tuples(
                                  tuple(zip(*[[u'bar', u'bar', u'baz',
                                               u'baz', u'baz'],
                                              [u'one', u'two', u'one',
                                               u'two', u'three']])),
                                  names=[u'first', u'second'])),
                 dup=DataFrame(np.arange(15).reshape(5, 3).astype(np.float64),
                               columns=[u'A', u'B', u'A']),
                 cat_onecol=DataFrame({u'A': Categorical([u'foo', u'bar'])}),
                 cat_and_float=DataFrame({
                     u'A': Categorical([u'foo', u'bar', u'baz']),
                     u'B': np.arange(3).astype(np.int64)}),
                 mixed_dup=mixed_dup_df,
                 dt_mixed_tzs=DataFrame({
                     u'A': Timestamp('20130102', tz='US/Eastern'),
                     u'B': Timestamp('20130603', tz='CET')}, index=range(5)),
                 dt_mixed2_tzs=DataFrame({
                     u'A': Timestamp('20130102', tz='US/Eastern'),
                     u'B': Timestamp('20130603', tz='CET'),
                     u'C': Timestamp('20130603', tz='UTC')}, index=range(5))
                 )

    with catch_warnings(record=True):
        mixed_dup_panel = Panel({u'ItemA': frame[u'float'],
                                 u'ItemB': frame[u'int']})
        mixed_dup_panel.items = [u'ItemA', u'ItemA']
        panel = dict(float=Panel({u'ItemA': frame[u'float'],
                                  u'ItemB': frame[u'float'] + 1}),
                     dup=Panel(
                         np.arange(30).reshape(3, 5, 2).astype(np.float64),
                         items=[u'A', u'B', u'A']),
                     mixed_dup=mixed_dup_panel)

    cat = dict(int8=Categorical(list('abcdefg')),
               int16=Categorical(np.arange(1000)),
               int32=Categorical(np.arange(10000)))

    timestamp = dict(normal=Timestamp('2011-01-01'),
                     nat=NaT,
                     tz=Timestamp('2011-01-01', tz='US/Eastern'))

    if _loose_version < LooseVersion('0.19.2'):
        timestamp['freq'] = Timestamp('2011-01-01', offset='D')
        timestamp['both'] = Timestamp('2011-01-01', tz='Asia/Tokyo',
                                      offset='M')
    else:
        timestamp['freq'] = Timestamp('2011-01-01', freq='D')
        timestamp['both'] = Timestamp('2011-01-01', tz='Asia/Tokyo',
                                      freq='M')

    off = {'DateOffset': DateOffset(years=1),
           'DateOffset_h_ns': DateOffset(hour=6, nanoseconds=5824),
           'BusinessDay': BusinessDay(offset=timedelta(seconds=9)),
           'BusinessHour': BusinessHour(normalize=True, n=6, end='15:14'),
           'CustomBusinessDay': CustomBusinessDay(weekmask='Mon Fri'),
           'SemiMonthBegin': SemiMonthBegin(day_of_month=9),
           'SemiMonthEnd': SemiMonthEnd(day_of_month=24),
           'MonthBegin': MonthBegin(1),
           'MonthEnd': MonthEnd(1),
           'QuarterBegin': QuarterBegin(1),
           'QuarterEnd': QuarterEnd(1),
           'Day': Day(1),
           'YearBegin': YearBegin(1),
           'YearEnd': YearEnd(1),
           'Week': Week(1),
           'Week_Tues': Week(2, normalize=False, weekday=1),
           'WeekOfMonth': WeekOfMonth(week=3, weekday=4),
           'LastWeekOfMonth': LastWeekOfMonth(n=1, weekday=3),
           'FY5253': FY5253(n=2, weekday=6, startingMonth=7, variation="last"),
           'Easter': Easter(),
           'Hour': Hour(1),
           'Minute': Minute(1)}

    return dict(series=series,
                frame=frame,
                panel=panel,
                index=index,
                scalars=scalars,
                mi=mi,
                sp_series=dict(float=_create_sp_series(),
                               ts=_create_sp_tsseries()),
                sp_frame=dict(float=_create_sp_frame()),
                cat=cat,
                timestamp=timestamp,
                offsets=off)


def create_pickle_data():
    data = create_data()

    # Pre-0.14.1 versions generated non-unpicklable mixed-type frames and
    # panels if their columns/items were non-unique.
    if _loose_version < LooseVersion('0.14.1'):
        del data['frame']['mixed_dup']
        del data['panel']['mixed_dup']
    if _loose_version < LooseVersion('0.17.0'):
        del data['series']['period']
        del data['scalars']['period']
    return data


def _u(x):
    return {u(k): _u(x[k]) for k in x} if isinstance(x, dict) else x


def create_msgpack_data():
    data = create_data()
    if _loose_version < LooseVersion('0.17.0'):
        del data['frame']['mixed_dup']
        del data['panel']['mixed_dup']
        del data['frame']['dup']
        del data['panel']['dup']
    if _loose_version < LooseVersion('0.18.0'):
        del data['series']['dt_tz']
        del data['frame']['dt_mixed_tzs']
    # Not supported
    del data['sp_series']
    del data['sp_frame']
    del data['series']['cat']
    del data['series']['period']
    del data['frame']['cat_onecol']
    del data['frame']['cat_and_float']
    del data['scalars']['period']
    if _loose_version < LooseVersion('0.23.0'):
        del data['index']['interval']
    del data['offsets']
    return _u(data)


def platform_name():
    return '_'.join([str(pandas.__version__), str(pl.machine()),
                     str(pl.system().lower()), str(pl.python_version())])


def write_legacy_pickles(output_dir):

    # make sure we are < 0.13 compat (in py3)
    try:
        from pandas.compat import zip, cPickle as pickle  # noqa
    except:
        import pickle

    version = pandas.__version__

    print("This script generates a storage file for the current arch, system, "
          "and python version")
    print("  pandas version: {0}".format(version))
    print("  output dir    : {0}".format(output_dir))
    print("  storage format: pickle")

    pth = '{0}.pickle'.format(platform_name())

    fh = open(os.path.join(output_dir, pth), 'wb')
    pickle.dump(create_pickle_data(), fh, pickle.HIGHEST_PROTOCOL)
    fh.close()

    print("created pickle file: %s" % pth)


def write_legacy_msgpack(output_dir, compress):

    version = pandas.__version__

    print("This script generates a storage file for the current arch, "
          "system, and python version")
    print("  pandas version: {0}".format(version))
    print("  output dir    : {0}".format(output_dir))
    print("  storage format: msgpack")
    pth = '{0}.msgpack'.format(platform_name())
    to_msgpack(os.path.join(output_dir, pth), create_msgpack_data(),
               compress=compress)

    print("created msgpack file: %s" % pth)


def write_legacy_file():
    # force our cwd to be the first searched
    sys.path.insert(0, '.')

    if not (3 <= len(sys.argv) <= 4):
        exit("Specify output directory and storage type: generate_legacy_"
             "storage_files.py <output_dir> <storage_type> "
             "<msgpack_compress_type>")

    output_dir = str(sys.argv[1])
    storage_type = str(sys.argv[2])
    try:
        compress_type = str(sys.argv[3])
    except IndexError:
        compress_type = None

    if storage_type == 'pickle':
        write_legacy_pickles(output_dir=output_dir)
    elif storage_type == 'msgpack':
        write_legacy_msgpack(output_dir=output_dir, compress=compress_type)
    else:
        exit("storage_type must be one of {'pickle', 'msgpack'}")


if __name__ == '__main__':
    write_legacy_file()
