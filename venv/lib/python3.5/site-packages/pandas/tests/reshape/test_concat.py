from warnings import catch_warnings
from itertools import combinations, product

import datetime as dt
import dateutil
import numpy as np
from numpy.random import randn

from datetime import datetime
from pandas.compat import StringIO, iteritems, PY2
import pandas as pd
from pandas import (DataFrame, concat,
                    read_csv, isna, Series, date_range,
                    Index, Panel, MultiIndex, Timestamp,
                    DatetimeIndex, Categorical)
from pandas.core.dtypes.dtypes import CategoricalDtype
from pandas.util import testing as tm
from pandas.util.testing import (assert_frame_equal,
                                 makeCustomDataframe as mkdf)

import pytest


@pytest.fixture(params=[True, False])
def sort(request):
    """Boolean sort keyword for concat and DataFrame.append."""
    return request.param


@pytest.fixture(params=[True, False, None])
def sort_with_none(request):
    """Boolean sort keyword for concat and DataFrame.append.

    Includes the default of None
    """
    # TODO: Replace with sort once keyword changes.
    return request.param


class ConcatenateBase(object):

    def setup_method(self, method):
        self.frame = DataFrame(tm.getSeriesData())
        self.mixed_frame = self.frame.copy()
        self.mixed_frame['foo'] = 'bar'


class TestConcatAppendCommon(ConcatenateBase):

    """
    Test common dtype coercion rules between concat and append.
    """

    def setup_method(self, method):

        dt_data = [pd.Timestamp('2011-01-01'),
                   pd.Timestamp('2011-01-02'),
                   pd.Timestamp('2011-01-03')]
        tz_data = [pd.Timestamp('2011-01-01', tz='US/Eastern'),
                   pd.Timestamp('2011-01-02', tz='US/Eastern'),
                   pd.Timestamp('2011-01-03', tz='US/Eastern')]

        td_data = [pd.Timedelta('1 days'),
                   pd.Timedelta('2 days'),
                   pd.Timedelta('3 days')]

        period_data = [pd.Period('2011-01', freq='M'),
                       pd.Period('2011-02', freq='M'),
                       pd.Period('2011-03', freq='M')]

        self.data = {'bool': [True, False, True],
                     'int64': [1, 2, 3],
                     'float64': [1.1, np.nan, 3.3],
                     'category': pd.Categorical(['X', 'Y', 'Z']),
                     'object': ['a', 'b', 'c'],
                     'datetime64[ns]': dt_data,
                     'datetime64[ns, US/Eastern]': tz_data,
                     'timedelta64[ns]': td_data,
                     'period[M]': period_data}

    def _check_expected_dtype(self, obj, label):
        """
        Check whether obj has expected dtype depending on label
        considering not-supported dtypes
        """
        if isinstance(obj, pd.Index):
            if label == 'bool':
                assert obj.dtype == 'object'
            else:
                assert obj.dtype == label
        elif isinstance(obj, pd.Series):
            if label.startswith('period'):
                assert obj.dtype == 'object'
            else:
                assert obj.dtype == label
        else:
            raise ValueError

    def test_dtypes(self):
        # to confirm test case covers intended dtypes
        for typ, vals in iteritems(self.data):
            self._check_expected_dtype(pd.Index(vals), typ)
            self._check_expected_dtype(pd.Series(vals), typ)

    def test_concatlike_same_dtypes(self):
        # GH 13660
        for typ1, vals1 in iteritems(self.data):

            vals2 = vals1
            vals3 = vals1

            if typ1 == 'category':
                exp_data = pd.Categorical(list(vals1) + list(vals2))
                exp_data3 = pd.Categorical(list(vals1) + list(vals2) +
                                           list(vals3))
            else:
                exp_data = vals1 + vals2
                exp_data3 = vals1 + vals2 + vals3

            # ----- Index ----- #

            # index.append
            res = pd.Index(vals1).append(pd.Index(vals2))
            exp = pd.Index(exp_data)
            tm.assert_index_equal(res, exp)

            # 3 elements
            res = pd.Index(vals1).append([pd.Index(vals2), pd.Index(vals3)])
            exp = pd.Index(exp_data3)
            tm.assert_index_equal(res, exp)

            # index.append name mismatch
            i1 = pd.Index(vals1, name='x')
            i2 = pd.Index(vals2, name='y')
            res = i1.append(i2)
            exp = pd.Index(exp_data)
            tm.assert_index_equal(res, exp)

            # index.append name match
            i1 = pd.Index(vals1, name='x')
            i2 = pd.Index(vals2, name='x')
            res = i1.append(i2)
            exp = pd.Index(exp_data, name='x')
            tm.assert_index_equal(res, exp)

            # cannot append non-index
            with tm.assert_raises_regex(TypeError,
                                        'all inputs must be Index'):
                pd.Index(vals1).append(vals2)

            with tm.assert_raises_regex(TypeError,
                                        'all inputs must be Index'):
                pd.Index(vals1).append([pd.Index(vals2), vals3])

            # ----- Series ----- #

            # series.append
            res = pd.Series(vals1).append(pd.Series(vals2),
                                          ignore_index=True)
            exp = pd.Series(exp_data)
            tm.assert_series_equal(res, exp, check_index_type=True)

            # concat
            res = pd.concat([pd.Series(vals1), pd.Series(vals2)],
                            ignore_index=True)
            tm.assert_series_equal(res, exp, check_index_type=True)

            # 3 elements
            res = pd.Series(vals1).append([pd.Series(vals2), pd.Series(vals3)],
                                          ignore_index=True)
            exp = pd.Series(exp_data3)
            tm.assert_series_equal(res, exp)

            res = pd.concat([pd.Series(vals1), pd.Series(vals2),
                             pd.Series(vals3)], ignore_index=True)
            tm.assert_series_equal(res, exp)

            # name mismatch
            s1 = pd.Series(vals1, name='x')
            s2 = pd.Series(vals2, name='y')
            res = s1.append(s2, ignore_index=True)
            exp = pd.Series(exp_data)
            tm.assert_series_equal(res, exp, check_index_type=True)

            res = pd.concat([s1, s2], ignore_index=True)
            tm.assert_series_equal(res, exp, check_index_type=True)

            # name match
            s1 = pd.Series(vals1, name='x')
            s2 = pd.Series(vals2, name='x')
            res = s1.append(s2, ignore_index=True)
            exp = pd.Series(exp_data, name='x')
            tm.assert_series_equal(res, exp, check_index_type=True)

            res = pd.concat([s1, s2], ignore_index=True)
            tm.assert_series_equal(res, exp, check_index_type=True)

            # cannot append non-index
            msg = (r'cannot concatenate object of type \"(.+?)\";'
                   ' only pd.Series, pd.DataFrame, and pd.Panel'
                   r' \(deprecated\) objs are valid')
            with tm.assert_raises_regex(TypeError, msg):
                pd.Series(vals1).append(vals2)

            with tm.assert_raises_regex(TypeError, msg):
                pd.Series(vals1).append([pd.Series(vals2), vals3])

            with tm.assert_raises_regex(TypeError, msg):
                pd.concat([pd.Series(vals1), vals2])

            with tm.assert_raises_regex(TypeError, msg):
                pd.concat([pd.Series(vals1), pd.Series(vals2), vals3])

    def test_concatlike_dtypes_coercion(self):
        # GH 13660
        for typ1, vals1 in iteritems(self.data):
            for typ2, vals2 in iteritems(self.data):

                vals3 = vals2

                # basically infer
                exp_index_dtype = None
                exp_series_dtype = None

                if typ1 == typ2:
                    # same dtype is tested in test_concatlike_same_dtypes
                    continue
                elif typ1 == 'category' or typ2 == 'category':
                    # ToDo: suspicious
                    continue

                # specify expected dtype
                if typ1 == 'bool' and typ2 in ('int64', 'float64'):
                    # series coerces to numeric based on numpy rule
                    # index doesn't because bool is object dtype
                    exp_series_dtype = typ2
                elif typ2 == 'bool' and typ1 in ('int64', 'float64'):
                    exp_series_dtype = typ1
                elif (typ1 == 'datetime64[ns, US/Eastern]' or
                      typ2 == 'datetime64[ns, US/Eastern]' or
                      typ1 == 'timedelta64[ns]' or
                      typ2 == 'timedelta64[ns]'):
                    exp_index_dtype = object
                    exp_series_dtype = object

                exp_data = vals1 + vals2
                exp_data3 = vals1 + vals2 + vals3

                # ----- Index ----- #

                # index.append
                res = pd.Index(vals1).append(pd.Index(vals2))
                exp = pd.Index(exp_data, dtype=exp_index_dtype)
                tm.assert_index_equal(res, exp)

                # 3 elements
                res = pd.Index(vals1).append([pd.Index(vals2),
                                              pd.Index(vals3)])
                exp = pd.Index(exp_data3, dtype=exp_index_dtype)
                tm.assert_index_equal(res, exp)

                # ----- Series ----- #

                # series.append
                res = pd.Series(vals1).append(pd.Series(vals2),
                                              ignore_index=True)
                exp = pd.Series(exp_data, dtype=exp_series_dtype)
                tm.assert_series_equal(res, exp, check_index_type=True)

                # concat
                res = pd.concat([pd.Series(vals1), pd.Series(vals2)],
                                ignore_index=True)
                tm.assert_series_equal(res, exp, check_index_type=True)

                # 3 elements
                res = pd.Series(vals1).append([pd.Series(vals2),
                                               pd.Series(vals3)],
                                              ignore_index=True)
                exp = pd.Series(exp_data3, dtype=exp_series_dtype)
                tm.assert_series_equal(res, exp)

                res = pd.concat([pd.Series(vals1), pd.Series(vals2),
                                 pd.Series(vals3)], ignore_index=True)
                tm.assert_series_equal(res, exp)

    def test_concatlike_common_coerce_to_pandas_object(self):
        # GH 13626
        # result must be Timestamp/Timedelta, not datetime.datetime/timedelta
        dti = pd.DatetimeIndex(['2011-01-01', '2011-01-02'])
        tdi = pd.TimedeltaIndex(['1 days', '2 days'])

        exp = pd.Index([pd.Timestamp('2011-01-01'),
                        pd.Timestamp('2011-01-02'),
                        pd.Timedelta('1 days'),
                        pd.Timedelta('2 days')])

        res = dti.append(tdi)
        tm.assert_index_equal(res, exp)
        assert isinstance(res[0], pd.Timestamp)
        assert isinstance(res[-1], pd.Timedelta)

        dts = pd.Series(dti)
        tds = pd.Series(tdi)
        res = dts.append(tds)
        tm.assert_series_equal(res, pd.Series(exp, index=[0, 1, 0, 1]))
        assert isinstance(res.iloc[0], pd.Timestamp)
        assert isinstance(res.iloc[-1], pd.Timedelta)

        res = pd.concat([dts, tds])
        tm.assert_series_equal(res, pd.Series(exp, index=[0, 1, 0, 1]))
        assert isinstance(res.iloc[0], pd.Timestamp)
        assert isinstance(res.iloc[-1], pd.Timedelta)

    def test_concatlike_datetimetz(self, tz_aware_fixture):
        tz = tz_aware_fixture
        # GH 7795
        dti1 = pd.DatetimeIndex(['2011-01-01', '2011-01-02'], tz=tz)
        dti2 = pd.DatetimeIndex(['2012-01-01', '2012-01-02'], tz=tz)

        exp = pd.DatetimeIndex(['2011-01-01', '2011-01-02',
                                '2012-01-01', '2012-01-02'], tz=tz)

        res = dti1.append(dti2)
        tm.assert_index_equal(res, exp)

        dts1 = pd.Series(dti1)
        dts2 = pd.Series(dti2)
        res = dts1.append(dts2)
        tm.assert_series_equal(res, pd.Series(exp, index=[0, 1, 0, 1]))

        res = pd.concat([dts1, dts2])
        tm.assert_series_equal(res, pd.Series(exp, index=[0, 1, 0, 1]))

    @pytest.mark.parametrize('tz',
                             ['UTC', 'US/Eastern', 'Asia/Tokyo', 'EST5EDT'])
    def test_concatlike_datetimetz_short(self, tz):
        # GH 7795
        ix1 = pd.DatetimeIndex(start='2014-07-15', end='2014-07-17',
                               freq='D', tz=tz)
        ix2 = pd.DatetimeIndex(['2014-07-11', '2014-07-21'], tz=tz)
        df1 = pd.DataFrame(0, index=ix1, columns=['A', 'B'])
        df2 = pd.DataFrame(0, index=ix2, columns=['A', 'B'])

        exp_idx = pd.DatetimeIndex(['2014-07-15', '2014-07-16',
                                    '2014-07-17', '2014-07-11',
                                    '2014-07-21'], tz=tz)
        exp = pd.DataFrame(0, index=exp_idx, columns=['A', 'B'])

        tm.assert_frame_equal(df1.append(df2), exp)
        tm.assert_frame_equal(pd.concat([df1, df2]), exp)

    def test_concatlike_datetimetz_to_object(self, tz_aware_fixture):
        tz = tz_aware_fixture
        # GH 13660

        # different tz coerces to object
        dti1 = pd.DatetimeIndex(['2011-01-01', '2011-01-02'], tz=tz)
        dti2 = pd.DatetimeIndex(['2012-01-01', '2012-01-02'])

        exp = pd.Index([pd.Timestamp('2011-01-01', tz=tz),
                        pd.Timestamp('2011-01-02', tz=tz),
                        pd.Timestamp('2012-01-01'),
                        pd.Timestamp('2012-01-02')], dtype=object)

        res = dti1.append(dti2)
        tm.assert_index_equal(res, exp)

        dts1 = pd.Series(dti1)
        dts2 = pd.Series(dti2)
        res = dts1.append(dts2)
        tm.assert_series_equal(res, pd.Series(exp, index=[0, 1, 0, 1]))

        res = pd.concat([dts1, dts2])
        tm.assert_series_equal(res, pd.Series(exp, index=[0, 1, 0, 1]))

        # different tz
        dti3 = pd.DatetimeIndex(['2012-01-01', '2012-01-02'],
                                tz='US/Pacific')

        exp = pd.Index([pd.Timestamp('2011-01-01', tz=tz),
                        pd.Timestamp('2011-01-02', tz=tz),
                        pd.Timestamp('2012-01-01', tz='US/Pacific'),
                        pd.Timestamp('2012-01-02', tz='US/Pacific')],
                       dtype=object)

        res = dti1.append(dti3)
        # tm.assert_index_equal(res, exp)

        dts1 = pd.Series(dti1)
        dts3 = pd.Series(dti3)
        res = dts1.append(dts3)
        tm.assert_series_equal(res, pd.Series(exp, index=[0, 1, 0, 1]))

        res = pd.concat([dts1, dts3])
        tm.assert_series_equal(res, pd.Series(exp, index=[0, 1, 0, 1]))

    def test_concatlike_common_period(self):
        # GH 13660
        pi1 = pd.PeriodIndex(['2011-01', '2011-02'], freq='M')
        pi2 = pd.PeriodIndex(['2012-01', '2012-02'], freq='M')

        exp = pd.PeriodIndex(['2011-01', '2011-02', '2012-01',
                              '2012-02'], freq='M')

        res = pi1.append(pi2)
        tm.assert_index_equal(res, exp)

        ps1 = pd.Series(pi1)
        ps2 = pd.Series(pi2)
        res = ps1.append(ps2)
        tm.assert_series_equal(res, pd.Series(exp, index=[0, 1, 0, 1]))

        res = pd.concat([ps1, ps2])
        tm.assert_series_equal(res, pd.Series(exp, index=[0, 1, 0, 1]))

    def test_concatlike_common_period_diff_freq_to_object(self):
        # GH 13221
        pi1 = pd.PeriodIndex(['2011-01', '2011-02'], freq='M')
        pi2 = pd.PeriodIndex(['2012-01-01', '2012-02-01'], freq='D')

        exp = pd.Index([pd.Period('2011-01', freq='M'),
                        pd.Period('2011-02', freq='M'),
                        pd.Period('2012-01-01', freq='D'),
                        pd.Period('2012-02-01', freq='D')], dtype=object)

        res = pi1.append(pi2)
        tm.assert_index_equal(res, exp)

        ps1 = pd.Series(pi1)
        ps2 = pd.Series(pi2)
        res = ps1.append(ps2)
        tm.assert_series_equal(res, pd.Series(exp, index=[0, 1, 0, 1]))

        res = pd.concat([ps1, ps2])
        tm.assert_series_equal(res, pd.Series(exp, index=[0, 1, 0, 1]))

    def test_concatlike_common_period_mixed_dt_to_object(self):
        # GH 13221
        # different datetimelike
        pi1 = pd.PeriodIndex(['2011-01', '2011-02'], freq='M')
        tdi = pd.TimedeltaIndex(['1 days', '2 days'])
        exp = pd.Index([pd.Period('2011-01', freq='M'),
                        pd.Period('2011-02', freq='M'),
                        pd.Timedelta('1 days'),
                        pd.Timedelta('2 days')], dtype=object)

        res = pi1.append(tdi)
        tm.assert_index_equal(res, exp)

        ps1 = pd.Series(pi1)
        tds = pd.Series(tdi)
        res = ps1.append(tds)
        tm.assert_series_equal(res, pd.Series(exp, index=[0, 1, 0, 1]))

        res = pd.concat([ps1, tds])
        tm.assert_series_equal(res, pd.Series(exp, index=[0, 1, 0, 1]))

        # inverse
        exp = pd.Index([pd.Timedelta('1 days'),
                        pd.Timedelta('2 days'),
                        pd.Period('2011-01', freq='M'),
                        pd.Period('2011-02', freq='M')], dtype=object)

        res = tdi.append(pi1)
        tm.assert_index_equal(res, exp)

        ps1 = pd.Series(pi1)
        tds = pd.Series(tdi)
        res = tds.append(ps1)
        tm.assert_series_equal(res, pd.Series(exp, index=[0, 1, 0, 1]))

        res = pd.concat([tds, ps1])
        tm.assert_series_equal(res, pd.Series(exp, index=[0, 1, 0, 1]))

    def test_concat_categorical(self):
        # GH 13524

        # same categories -> category
        s1 = pd.Series([1, 2, np.nan], dtype='category')
        s2 = pd.Series([2, 1, 2], dtype='category')

        exp = pd.Series([1, 2, np.nan, 2, 1, 2], dtype='category')
        tm.assert_series_equal(pd.concat([s1, s2], ignore_index=True), exp)
        tm.assert_series_equal(s1.append(s2, ignore_index=True), exp)

        # partially different categories => not-category
        s1 = pd.Series([3, 2], dtype='category')
        s2 = pd.Series([2, 1], dtype='category')

        exp = pd.Series([3, 2, 2, 1])
        tm.assert_series_equal(pd.concat([s1, s2], ignore_index=True), exp)
        tm.assert_series_equal(s1.append(s2, ignore_index=True), exp)

        # completely different categories (same dtype) => not-category
        s1 = pd.Series([10, 11, np.nan], dtype='category')
        s2 = pd.Series([np.nan, 1, 3, 2], dtype='category')

        exp = pd.Series([10, 11, np.nan, np.nan, 1, 3, 2])
        tm.assert_series_equal(pd.concat([s1, s2], ignore_index=True), exp)
        tm.assert_series_equal(s1.append(s2, ignore_index=True), exp)

    def test_union_categorical_same_categories_different_order(self):
        # https://github.com/pandas-dev/pandas/issues/19096
        a = pd.Series(Categorical(['a', 'b', 'c'], categories=['a', 'b', 'c']))
        b = pd.Series(Categorical(['a', 'b', 'c'], categories=['b', 'a', 'c']))
        result = pd.concat([a, b], ignore_index=True)
        expected = pd.Series(Categorical(['a', 'b', 'c', 'a', 'b', 'c'],
                                         categories=['a', 'b', 'c']))
        tm.assert_series_equal(result, expected)

    def test_concat_categorical_coercion(self):
        # GH 13524

        # category + not-category => not-category
        s1 = pd.Series([1, 2, np.nan], dtype='category')
        s2 = pd.Series([2, 1, 2])

        exp = pd.Series([1, 2, np.nan, 2, 1, 2])
        tm.assert_series_equal(pd.concat([s1, s2], ignore_index=True), exp)
        tm.assert_series_equal(s1.append(s2, ignore_index=True), exp)

        # result shouldn't be affected by 1st elem dtype
        exp = pd.Series([2, 1, 2, 1, 2, np.nan])
        tm.assert_series_equal(pd.concat([s2, s1], ignore_index=True), exp)
        tm.assert_series_equal(s2.append(s1, ignore_index=True), exp)

        # all values are not in category => not-category
        s1 = pd.Series([3, 2], dtype='category')
        s2 = pd.Series([2, 1])

        exp = pd.Series([3, 2, 2, 1])
        tm.assert_series_equal(pd.concat([s1, s2], ignore_index=True), exp)
        tm.assert_series_equal(s1.append(s2, ignore_index=True), exp)

        exp = pd.Series([2, 1, 3, 2])
        tm.assert_series_equal(pd.concat([s2, s1], ignore_index=True), exp)
        tm.assert_series_equal(s2.append(s1, ignore_index=True), exp)

        # completely different categories => not-category
        s1 = pd.Series([10, 11, np.nan], dtype='category')
        s2 = pd.Series([1, 3, 2])

        exp = pd.Series([10, 11, np.nan, 1, 3, 2])
        tm.assert_series_equal(pd.concat([s1, s2], ignore_index=True), exp)
        tm.assert_series_equal(s1.append(s2, ignore_index=True), exp)

        exp = pd.Series([1, 3, 2, 10, 11, np.nan])
        tm.assert_series_equal(pd.concat([s2, s1], ignore_index=True), exp)
        tm.assert_series_equal(s2.append(s1, ignore_index=True), exp)

        # different dtype => not-category
        s1 = pd.Series([10, 11, np.nan], dtype='category')
        s2 = pd.Series(['a', 'b', 'c'])

        exp = pd.Series([10, 11, np.nan, 'a', 'b', 'c'])
        tm.assert_series_equal(pd.concat([s1, s2], ignore_index=True), exp)
        tm.assert_series_equal(s1.append(s2, ignore_index=True), exp)

        exp = pd.Series(['a', 'b', 'c', 10, 11, np.nan])
        tm.assert_series_equal(pd.concat([s2, s1], ignore_index=True), exp)
        tm.assert_series_equal(s2.append(s1, ignore_index=True), exp)

        # if normal series only contains NaN-likes => not-category
        s1 = pd.Series([10, 11], dtype='category')
        s2 = pd.Series([np.nan, np.nan, np.nan])

        exp = pd.Series([10, 11, np.nan, np.nan, np.nan])
        tm.assert_series_equal(pd.concat([s1, s2], ignore_index=True), exp)
        tm.assert_series_equal(s1.append(s2, ignore_index=True), exp)

        exp = pd.Series([np.nan, np.nan, np.nan, 10, 11])
        tm.assert_series_equal(pd.concat([s2, s1], ignore_index=True), exp)
        tm.assert_series_equal(s2.append(s1, ignore_index=True), exp)

    def test_concat_categorical_3elem_coercion(self):
        # GH 13524

        # mixed dtypes => not-category
        s1 = pd.Series([1, 2, np.nan], dtype='category')
        s2 = pd.Series([2, 1, 2], dtype='category')
        s3 = pd.Series([1, 2, 1, 2, np.nan])

        exp = pd.Series([1, 2, np.nan, 2, 1, 2, 1, 2, 1, 2, np.nan])
        tm.assert_series_equal(pd.concat([s1, s2, s3], ignore_index=True), exp)
        tm.assert_series_equal(s1.append([s2, s3], ignore_index=True), exp)

        exp = pd.Series([1, 2, 1, 2, np.nan, 1, 2, np.nan, 2, 1, 2])
        tm.assert_series_equal(pd.concat([s3, s1, s2], ignore_index=True), exp)
        tm.assert_series_equal(s3.append([s1, s2], ignore_index=True), exp)

        # values are all in either category => not-category
        s1 = pd.Series([4, 5, 6], dtype='category')
        s2 = pd.Series([1, 2, 3], dtype='category')
        s3 = pd.Series([1, 3, 4])

        exp = pd.Series([4, 5, 6, 1, 2, 3, 1, 3, 4])
        tm.assert_series_equal(pd.concat([s1, s2, s3], ignore_index=True), exp)
        tm.assert_series_equal(s1.append([s2, s3], ignore_index=True), exp)

        exp = pd.Series([1, 3, 4, 4, 5, 6, 1, 2, 3])
        tm.assert_series_equal(pd.concat([s3, s1, s2], ignore_index=True), exp)
        tm.assert_series_equal(s3.append([s1, s2], ignore_index=True), exp)

        # values are all in either category => not-category
        s1 = pd.Series([4, 5, 6], dtype='category')
        s2 = pd.Series([1, 2, 3], dtype='category')
        s3 = pd.Series([10, 11, 12])

        exp = pd.Series([4, 5, 6, 1, 2, 3, 10, 11, 12])
        tm.assert_series_equal(pd.concat([s1, s2, s3], ignore_index=True), exp)
        tm.assert_series_equal(s1.append([s2, s3], ignore_index=True), exp)

        exp = pd.Series([10, 11, 12, 4, 5, 6, 1, 2, 3])
        tm.assert_series_equal(pd.concat([s3, s1, s2], ignore_index=True), exp)
        tm.assert_series_equal(s3.append([s1, s2], ignore_index=True), exp)

    def test_concat_categorical_multi_coercion(self):
        # GH 13524

        s1 = pd.Series([1, 3], dtype='category')
        s2 = pd.Series([3, 4], dtype='category')
        s3 = pd.Series([2, 3])
        s4 = pd.Series([2, 2], dtype='category')
        s5 = pd.Series([1, np.nan])
        s6 = pd.Series([1, 3, 2], dtype='category')

        # mixed dtype, values are all in categories => not-category
        exp = pd.Series([1, 3, 3, 4, 2, 3, 2, 2, 1, np.nan, 1, 3, 2])
        res = pd.concat([s1, s2, s3, s4, s5, s6], ignore_index=True)
        tm.assert_series_equal(res, exp)
        res = s1.append([s2, s3, s4, s5, s6], ignore_index=True)
        tm.assert_series_equal(res, exp)

        exp = pd.Series([1, 3, 2, 1, np.nan, 2, 2, 2, 3, 3, 4, 1, 3])
        res = pd.concat([s6, s5, s4, s3, s2, s1], ignore_index=True)
        tm.assert_series_equal(res, exp)
        res = s6.append([s5, s4, s3, s2, s1], ignore_index=True)
        tm.assert_series_equal(res, exp)

    def test_concat_categorical_ordered(self):
        # GH 13524

        s1 = pd.Series(pd.Categorical([1, 2, np.nan], ordered=True))
        s2 = pd.Series(pd.Categorical([2, 1, 2], ordered=True))

        exp = pd.Series(pd.Categorical([1, 2, np.nan, 2, 1, 2], ordered=True))
        tm.assert_series_equal(pd.concat([s1, s2], ignore_index=True), exp)
        tm.assert_series_equal(s1.append(s2, ignore_index=True), exp)

        exp = pd.Series(pd.Categorical([1, 2, np.nan, 2, 1, 2, 1, 2, np.nan],
                                       ordered=True))
        tm.assert_series_equal(pd.concat([s1, s2, s1], ignore_index=True), exp)
        tm.assert_series_equal(s1.append([s2, s1], ignore_index=True), exp)

    def test_concat_categorical_coercion_nan(self):
        # GH 13524

        # some edge cases
        # category + not-category => not category
        s1 = pd.Series(np.array([np.nan, np.nan], dtype=np.float64),
                       dtype='category')
        s2 = pd.Series([np.nan, 1])

        exp = pd.Series([np.nan, np.nan, np.nan, 1])
        tm.assert_series_equal(pd.concat([s1, s2], ignore_index=True), exp)
        tm.assert_series_equal(s1.append(s2, ignore_index=True), exp)

        s1 = pd.Series([1, np.nan], dtype='category')
        s2 = pd.Series([np.nan, np.nan])

        exp = pd.Series([1, np.nan, np.nan, np.nan])
        tm.assert_series_equal(pd.concat([s1, s2], ignore_index=True), exp)
        tm.assert_series_equal(s1.append(s2, ignore_index=True), exp)

        # mixed dtype, all nan-likes => not-category
        s1 = pd.Series([np.nan, np.nan], dtype='category')
        s2 = pd.Series([np.nan, np.nan])

        exp = pd.Series([np.nan, np.nan, np.nan, np.nan])
        tm.assert_series_equal(pd.concat([s1, s2], ignore_index=True), exp)
        tm.assert_series_equal(s1.append(s2, ignore_index=True), exp)
        tm.assert_series_equal(pd.concat([s2, s1], ignore_index=True), exp)
        tm.assert_series_equal(s2.append(s1, ignore_index=True), exp)

        # all category nan-likes => category
        s1 = pd.Series([np.nan, np.nan], dtype='category')
        s2 = pd.Series([np.nan, np.nan], dtype='category')

        exp = pd.Series([np.nan, np.nan, np.nan, np.nan], dtype='category')

        tm.assert_series_equal(pd.concat([s1, s2], ignore_index=True), exp)
        tm.assert_series_equal(s1.append(s2, ignore_index=True), exp)

    def test_concat_categorical_empty(self):
        # GH 13524

        s1 = pd.Series([], dtype='category')
        s2 = pd.Series([1, 2], dtype='category')

        tm.assert_series_equal(pd.concat([s1, s2], ignore_index=True), s2)
        tm.assert_series_equal(s1.append(s2, ignore_index=True), s2)

        tm.assert_series_equal(pd.concat([s2, s1], ignore_index=True), s2)
        tm.assert_series_equal(s2.append(s1, ignore_index=True), s2)

        s1 = pd.Series([], dtype='category')
        s2 = pd.Series([], dtype='category')

        tm.assert_series_equal(pd.concat([s1, s2], ignore_index=True), s2)
        tm.assert_series_equal(s1.append(s2, ignore_index=True), s2)

        s1 = pd.Series([], dtype='category')
        s2 = pd.Series([], dtype='object')

        # different dtype => not-category
        tm.assert_series_equal(pd.concat([s1, s2], ignore_index=True), s2)
        tm.assert_series_equal(s1.append(s2, ignore_index=True), s2)
        tm.assert_series_equal(pd.concat([s2, s1], ignore_index=True), s2)
        tm.assert_series_equal(s2.append(s1, ignore_index=True), s2)

        s1 = pd.Series([], dtype='category')
        s2 = pd.Series([np.nan, np.nan])

        # empty Series is ignored
        exp = pd.Series([np.nan, np.nan])
        tm.assert_series_equal(pd.concat([s1, s2], ignore_index=True), exp)
        tm.assert_series_equal(s1.append(s2, ignore_index=True), exp)

        tm.assert_series_equal(pd.concat([s2, s1], ignore_index=True), exp)
        tm.assert_series_equal(s2.append(s1, ignore_index=True), exp)


class TestAppend(ConcatenateBase):

    def test_append(self, sort):
        begin_index = self.frame.index[:5]
        end_index = self.frame.index[5:]

        begin_frame = self.frame.reindex(begin_index)
        end_frame = self.frame.reindex(end_index)

        appended = begin_frame.append(end_frame)
        tm.assert_almost_equal(appended['A'], self.frame['A'])

        del end_frame['A']
        partial_appended = begin_frame.append(end_frame, sort=sort)
        assert 'A' in partial_appended

        partial_appended = end_frame.append(begin_frame, sort=sort)
        assert 'A' in partial_appended

        # mixed type handling
        appended = self.mixed_frame[:5].append(self.mixed_frame[5:])
        tm.assert_frame_equal(appended, self.mixed_frame)

        # what to test here
        mixed_appended = self.mixed_frame[:5].append(self.frame[5:], sort=sort)
        mixed_appended2 = self.frame[:5].append(self.mixed_frame[5:],
                                                sort=sort)

        # all equal except 'foo' column
        tm.assert_frame_equal(
            mixed_appended.reindex(columns=['A', 'B', 'C', 'D']),
            mixed_appended2.reindex(columns=['A', 'B', 'C', 'D']))

        # append empty
        empty = DataFrame({})

        appended = self.frame.append(empty)
        tm.assert_frame_equal(self.frame, appended)
        assert appended is not self.frame

        appended = empty.append(self.frame)
        tm.assert_frame_equal(self.frame, appended)
        assert appended is not self.frame

        # Overlap
        with pytest.raises(ValueError):
            self.frame.append(self.frame, verify_integrity=True)

        # see gh-6129: new columns
        df = DataFrame({'a': {'x': 1, 'y': 2}, 'b': {'x': 3, 'y': 4}})
        row = Series([5, 6, 7], index=['a', 'b', 'c'], name='z')
        expected = DataFrame({'a': {'x': 1, 'y': 2, 'z': 5}, 'b': {
                             'x': 3, 'y': 4, 'z': 6}, 'c': {'z': 7}})
        result = df.append(row)
        tm.assert_frame_equal(result, expected)

    def test_append_length0_frame(self, sort):
        df = DataFrame(columns=['A', 'B', 'C'])
        df3 = DataFrame(index=[0, 1], columns=['A', 'B'])
        df5 = df.append(df3, sort=sort)

        expected = DataFrame(index=[0, 1], columns=['A', 'B', 'C'])
        assert_frame_equal(df5, expected)

    def test_append_records(self):
        arr1 = np.zeros((2,), dtype=('i4,f4,a10'))
        arr1[:] = [(1, 2., 'Hello'), (2, 3., "World")]

        arr2 = np.zeros((3,), dtype=('i4,f4,a10'))
        arr2[:] = [(3, 4., 'foo'),
                   (5, 6., "bar"),
                   (7., 8., 'baz')]

        df1 = DataFrame(arr1)
        df2 = DataFrame(arr2)

        result = df1.append(df2, ignore_index=True)
        expected = DataFrame(np.concatenate((arr1, arr2)))
        assert_frame_equal(result, expected)

    # rewrite sort fixture, since we also want to test default of None
    def test_append_sorts(self, sort_with_none):
        df1 = pd.DataFrame({"a": [1, 2], "b": [1, 2]}, columns=['b', 'a'])
        df2 = pd.DataFrame({"a": [1, 2], 'c': [3, 4]}, index=[2, 3])

        if sort_with_none is None:
            # only warn if not explicitly specified
            # don't check stacklevel since its set for concat, and append
            # has an extra stack.
            ctx = tm.assert_produces_warning(FutureWarning,
                                             check_stacklevel=False)
        else:
            ctx = tm.assert_produces_warning(None)

        with ctx:
            result = df1.append(df2, sort=sort_with_none)

        # for None / True
        expected = pd.DataFrame({"b": [1, 2, None, None],
                                 "a": [1, 2, 1, 2],
                                 "c": [None, None, 3, 4]},
                                columns=['a', 'b', 'c'])
        if sort_with_none is False:
            expected = expected[['b', 'a', 'c']]
        tm.assert_frame_equal(result, expected)

    def test_append_different_columns(self, sort):
        df = DataFrame({'bools': np.random.randn(10) > 0,
                        'ints': np.random.randint(0, 10, 10),
                        'floats': np.random.randn(10),
                        'strings': ['foo', 'bar'] * 5})

        a = df[:5].loc[:, ['bools', 'ints', 'floats']]
        b = df[5:].loc[:, ['strings', 'ints', 'floats']]

        appended = a.append(b, sort=sort)
        assert isna(appended['strings'][0:4]).all()
        assert isna(appended['bools'][5:]).all()

    def test_append_many(self, sort):
        chunks = [self.frame[:5], self.frame[5:10],
                  self.frame[10:15], self.frame[15:]]

        result = chunks[0].append(chunks[1:])
        tm.assert_frame_equal(result, self.frame)

        chunks[-1] = chunks[-1].copy()
        chunks[-1]['foo'] = 'bar'
        result = chunks[0].append(chunks[1:], sort=sort)
        tm.assert_frame_equal(result.loc[:, self.frame.columns], self.frame)
        assert (result['foo'][15:] == 'bar').all()
        assert result['foo'][:15].isna().all()

    def test_append_preserve_index_name(self):
        # #980
        df1 = DataFrame(data=None, columns=['A', 'B', 'C'])
        df1 = df1.set_index(['A'])
        df2 = DataFrame(data=[[1, 4, 7], [2, 5, 8], [3, 6, 9]],
                        columns=['A', 'B', 'C'])
        df2 = df2.set_index(['A'])

        result = df1.append(df2)
        assert result.index.name == 'A'

    indexes_can_append = [
        pd.RangeIndex(3),
        pd.Index([4, 5, 6]),
        pd.Index([4.5, 5.5, 6.5]),
        pd.Index(list('abc')),
        pd.CategoricalIndex('A B C'.split()),
        pd.CategoricalIndex('D E F'.split(), ordered=True),
        pd.DatetimeIndex([dt.datetime(2013, 1, 3, 0, 0),
                          dt.datetime(2013, 1, 3, 6, 10),
                          dt.datetime(2013, 1, 3, 7, 12)]),
    ]

    indexes_cannot_append_with_other = [
        pd.IntervalIndex.from_breaks([0, 1, 2, 3]),
        pd.MultiIndex.from_arrays(['A B C'.split(), 'D E F'.split()]),
    ]

    all_indexes = indexes_can_append + indexes_cannot_append_with_other

    @pytest.mark.parametrize("index",
                             all_indexes,
                             ids=lambda x: x.__class__.__name__)
    def test_append_same_columns_type(self, index):
        # GH18359

        # df wider than ser
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], columns=index)
        ser_index = index[:2]
        ser = pd.Series([7, 8], index=ser_index, name=2)
        result = df.append(ser)
        expected = pd.DataFrame([[1., 2., 3.], [4, 5, 6], [7, 8, np.nan]],
                                index=[0, 1, 2],
                                columns=index)
        assert_frame_equal(result, expected)

        # ser wider than df
        ser_index = index
        index = index[:2]
        df = pd.DataFrame([[1, 2], [4, 5]], columns=index)
        ser = pd.Series([7, 8, 9], index=ser_index, name=2)
        result = df.append(ser)
        expected = pd.DataFrame([[1, 2, np.nan], [4, 5, np.nan], [7, 8, 9]],
                                index=[0, 1, 2],
                                columns=ser_index)
        assert_frame_equal(result, expected)

    @pytest.mark.parametrize("df_columns, series_index",
                             combinations(indexes_can_append, r=2),
                             ids=lambda x: x.__class__.__name__)
    def test_append_different_columns_types(self, df_columns, series_index):
        # GH18359
        # See also test 'test_append_different_columns_types_raises' below
        # for errors raised when appending

        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], columns=df_columns)
        ser = pd.Series([7, 8, 9], index=series_index, name=2)

        result = df.append(ser)
        idx_diff = ser.index.difference(df_columns)
        combined_columns = Index(df_columns.tolist()).append(idx_diff)
        expected = pd.DataFrame([[1., 2., 3., np.nan, np.nan, np.nan],
                                 [4, 5, 6, np.nan, np.nan, np.nan],
                                 [np.nan, np.nan, np.nan, 7, 8, 9]],
                                index=[0, 1, 2],
                                columns=combined_columns)
        assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "index_can_append, index_cannot_append_with_other",
        product(indexes_can_append, indexes_cannot_append_with_other),
        ids=lambda x: x.__class__.__name__)
    def test_append_different_columns_types_raises(
            self, index_can_append, index_cannot_append_with_other):
        # GH18359
        # Dataframe.append will raise if IntervalIndex/MultiIndex appends
        # or is appended to a different index type
        #
        # See also test 'test_append_different_columns_types' above for
        # appending without raising.

        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], columns=index_can_append)
        ser = pd.Series([7, 8, 9], index=index_cannot_append_with_other,
                        name=2)
        with pytest.raises(TypeError):
            df.append(ser)

        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]],
                          columns=index_cannot_append_with_other)
        ser = pd.Series([7, 8, 9], index=index_can_append, name=2)
        with pytest.raises(TypeError):
            df.append(ser)

    def test_append_dtype_coerce(self, sort):

        # GH 4993
        # appending with datetime will incorrectly convert datetime64

        df1 = DataFrame(index=[1, 2], data=[dt.datetime(2013, 1, 1, 0, 0),
                                            dt.datetime(2013, 1, 2, 0, 0)],
                        columns=['start_time'])
        df2 = DataFrame(index=[4, 5], data=[[dt.datetime(2013, 1, 3, 0, 0),
                                             dt.datetime(2013, 1, 3, 6, 10)],
                                            [dt.datetime(2013, 1, 4, 0, 0),
                                             dt.datetime(2013, 1, 4, 7, 10)]],
                        columns=['start_time', 'end_time'])

        expected = concat([Series([pd.NaT,
                                   pd.NaT,
                                   dt.datetime(2013, 1, 3, 6, 10),
                                   dt.datetime(2013, 1, 4, 7, 10)],
                                  name='end_time'),
                           Series([dt.datetime(2013, 1, 1, 0, 0),
                                   dt.datetime(2013, 1, 2, 0, 0),
                                   dt.datetime(2013, 1, 3, 0, 0),
                                   dt.datetime(2013, 1, 4, 0, 0)],
                                  name='start_time')],
                          axis=1, sort=sort)
        result = df1.append(df2, ignore_index=True, sort=sort)
        if sort:
            expected = expected[['end_time', 'start_time']]
        else:
            expected = expected[['start_time', 'end_time']]

        assert_frame_equal(result, expected)

    def test_append_missing_column_proper_upcast(self, sort):
        df1 = DataFrame({'A': np.array([1, 2, 3, 4], dtype='i8')})
        df2 = DataFrame({'B': np.array([True, False, True, False],
                                       dtype=bool)})

        appended = df1.append(df2, ignore_index=True, sort=sort)
        assert appended['A'].dtype == 'f8'
        assert appended['B'].dtype == 'O'


class TestConcatenate(ConcatenateBase):

    def test_concat_copy(self):
        df = DataFrame(np.random.randn(4, 3))
        df2 = DataFrame(np.random.randint(0, 10, size=4).reshape(4, 1))
        df3 = DataFrame({5: 'foo'}, index=range(4))

        # These are actual copies.
        result = concat([df, df2, df3], axis=1, copy=True)

        for b in result._data.blocks:
            assert b.values.base is None

        # These are the same.
        result = concat([df, df2, df3], axis=1, copy=False)

        for b in result._data.blocks:
            if b.is_float:
                assert b.values.base is df._data.blocks[0].values.base
            elif b.is_integer:
                assert b.values.base is df2._data.blocks[0].values.base
            elif b.is_object:
                assert b.values.base is not None

        # Float block was consolidated.
        df4 = DataFrame(np.random.randn(4, 1))
        result = concat([df, df2, df3, df4], axis=1, copy=False)
        for b in result._data.blocks:
            if b.is_float:
                assert b.values.base is None
            elif b.is_integer:
                assert b.values.base is df2._data.blocks[0].values.base
            elif b.is_object:
                assert b.values.base is not None

    def test_concat_with_group_keys(self):
        df = DataFrame(np.random.randn(4, 3))
        df2 = DataFrame(np.random.randn(4, 4))

        # axis=0
        df = DataFrame(np.random.randn(3, 4))
        df2 = DataFrame(np.random.randn(4, 4))

        result = concat([df, df2], keys=[0, 1])
        exp_index = MultiIndex.from_arrays([[0, 0, 0, 1, 1, 1, 1],
                                            [0, 1, 2, 0, 1, 2, 3]])
        expected = DataFrame(np.r_[df.values, df2.values],
                             index=exp_index)
        tm.assert_frame_equal(result, expected)

        result = concat([df, df], keys=[0, 1])
        exp_index2 = MultiIndex.from_arrays([[0, 0, 0, 1, 1, 1],
                                             [0, 1, 2, 0, 1, 2]])
        expected = DataFrame(np.r_[df.values, df.values],
                             index=exp_index2)
        tm.assert_frame_equal(result, expected)

        # axis=1
        df = DataFrame(np.random.randn(4, 3))
        df2 = DataFrame(np.random.randn(4, 4))

        result = concat([df, df2], keys=[0, 1], axis=1)
        expected = DataFrame(np.c_[df.values, df2.values],
                             columns=exp_index)
        tm.assert_frame_equal(result, expected)

        result = concat([df, df], keys=[0, 1], axis=1)
        expected = DataFrame(np.c_[df.values, df.values],
                             columns=exp_index2)
        tm.assert_frame_equal(result, expected)

    def test_concat_keys_specific_levels(self):
        df = DataFrame(np.random.randn(10, 4))
        pieces = [df.iloc[:, [0, 1]], df.iloc[:, [2]], df.iloc[:, [3]]]
        level = ['three', 'two', 'one', 'zero']
        result = concat(pieces, axis=1, keys=['one', 'two', 'three'],
                        levels=[level],
                        names=['group_key'])

        tm.assert_index_equal(result.columns.levels[0],
                              Index(level, name='group_key'))
        assert result.columns.names[0] == 'group_key'

    def test_concat_dataframe_keys_bug(self, sort):
        t1 = DataFrame({
            'value': Series([1, 2, 3], index=Index(['a', 'b', 'c'],
                                                   name='id'))})
        t2 = DataFrame({
            'value': Series([7, 8], index=Index(['a', 'b'], name='id'))})

        # it works
        result = concat([t1, t2], axis=1, keys=['t1', 't2'], sort=sort)
        assert list(result.columns) == [('t1', 'value'), ('t2', 'value')]

    def test_concat_series_partial_columns_names(self):
        # GH10698
        foo = Series([1, 2], name='foo')
        bar = Series([1, 2])
        baz = Series([4, 5])

        result = concat([foo, bar, baz], axis=1)
        expected = DataFrame({'foo': [1, 2], 0: [1, 2], 1: [
                             4, 5]}, columns=['foo', 0, 1])
        tm.assert_frame_equal(result, expected)

        result = concat([foo, bar, baz], axis=1, keys=[
                        'red', 'blue', 'yellow'])
        expected = DataFrame({'red': [1, 2], 'blue': [1, 2], 'yellow': [
                             4, 5]}, columns=['red', 'blue', 'yellow'])
        tm.assert_frame_equal(result, expected)

        result = concat([foo, bar, baz], axis=1, ignore_index=True)
        expected = DataFrame({0: [1, 2], 1: [1, 2], 2: [4, 5]})
        tm.assert_frame_equal(result, expected)

    def test_concat_dict(self):
        frames = {'foo': DataFrame(np.random.randn(4, 3)),
                  'bar': DataFrame(np.random.randn(4, 3)),
                  'baz': DataFrame(np.random.randn(4, 3)),
                  'qux': DataFrame(np.random.randn(4, 3))}

        sorted_keys = sorted(frames)

        result = concat(frames)
        expected = concat([frames[k] for k in sorted_keys], keys=sorted_keys)
        tm.assert_frame_equal(result, expected)

        result = concat(frames, axis=1)
        expected = concat([frames[k] for k in sorted_keys], keys=sorted_keys,
                          axis=1)
        tm.assert_frame_equal(result, expected)

        keys = ['baz', 'foo', 'bar']
        result = concat(frames, keys=keys)
        expected = concat([frames[k] for k in keys], keys=keys)
        tm.assert_frame_equal(result, expected)

    def test_concat_ignore_index(self, sort):
        frame1 = DataFrame({"test1": ["a", "b", "c"],
                            "test2": [1, 2, 3],
                            "test3": [4.5, 3.2, 1.2]})
        frame2 = DataFrame({"test3": [5.2, 2.2, 4.3]})
        frame1.index = Index(["x", "y", "z"])
        frame2.index = Index(["x", "y", "q"])

        v1 = concat([frame1, frame2], axis=1,
                    ignore_index=True, sort=sort)

        nan = np.nan
        expected = DataFrame([[nan, nan, nan, 4.3],
                              ['a', 1, 4.5, 5.2],
                              ['b', 2, 3.2, 2.2],
                              ['c', 3, 1.2, nan]],
                             index=Index(["q", "x", "y", "z"]))
        if not sort:
            expected = expected.loc[['x', 'y', 'z', 'q']]

        tm.assert_frame_equal(v1, expected)

    def test_concat_multiindex_with_keys(self):
        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])
        frame = DataFrame(np.random.randn(10, 3), index=index,
                          columns=Index(['A', 'B', 'C'], name='exp'))
        result = concat([frame, frame], keys=[0, 1], names=['iteration'])

        assert result.index.names == ('iteration',) + index.names
        tm.assert_frame_equal(result.loc[0], frame)
        tm.assert_frame_equal(result.loc[1], frame)
        assert result.index.nlevels == 3

    def test_concat_multiindex_with_tz(self):
        # GH 6606
        df = DataFrame({'dt': [datetime(2014, 1, 1),
                               datetime(2014, 1, 2),
                               datetime(2014, 1, 3)],
                        'b': ['A', 'B', 'C'],
                        'c': [1, 2, 3], 'd': [4, 5, 6]})
        df['dt'] = df['dt'].apply(lambda d: Timestamp(d, tz='US/Pacific'))
        df = df.set_index(['dt', 'b'])

        exp_idx1 = DatetimeIndex(['2014-01-01', '2014-01-02',
                                  '2014-01-03'] * 2,
                                 tz='US/Pacific', name='dt')
        exp_idx2 = Index(['A', 'B', 'C'] * 2, name='b')
        exp_idx = MultiIndex.from_arrays([exp_idx1, exp_idx2])
        expected = DataFrame({'c': [1, 2, 3] * 2, 'd': [4, 5, 6] * 2},
                             index=exp_idx, columns=['c', 'd'])

        result = concat([df, df])
        tm.assert_frame_equal(result, expected)

    def test_concat_multiindex_with_none_in_index_names(self):
        # GH 15787
        index = pd.MultiIndex.from_product([[1], range(5)],
                                           names=['level1', None])
        df = pd.DataFrame({'col': range(5)}, index=index, dtype=np.int32)

        result = concat([df, df], keys=[1, 2], names=['level2'])
        index = pd.MultiIndex.from_product([[1, 2], [1], range(5)],
                                           names=['level2', 'level1', None])
        expected = pd.DataFrame({'col': list(range(5)) * 2},
                                index=index, dtype=np.int32)
        assert_frame_equal(result, expected)

        result = concat([df, df[:2]], keys=[1, 2], names=['level2'])
        level2 = [1] * 5 + [2] * 2
        level1 = [1] * 7
        no_name = list(range(5)) + list(range(2))
        tuples = list(zip(level2, level1, no_name))
        index = pd.MultiIndex.from_tuples(tuples,
                                          names=['level2', 'level1', None])
        expected = pd.DataFrame({'col': no_name}, index=index,
                                dtype=np.int32)
        assert_frame_equal(result, expected)

    def test_concat_keys_and_levels(self):
        df = DataFrame(np.random.randn(1, 3))
        df2 = DataFrame(np.random.randn(1, 4))

        levels = [['foo', 'baz'], ['one', 'two']]
        names = ['first', 'second']
        result = concat([df, df2, df, df2],
                        keys=[('foo', 'one'), ('foo', 'two'),
                              ('baz', 'one'), ('baz', 'two')],
                        levels=levels,
                        names=names)
        expected = concat([df, df2, df, df2])
        exp_index = MultiIndex(levels=levels + [[0]],
                               labels=[[0, 0, 1, 1], [0, 1, 0, 1],
                                       [0, 0, 0, 0]],
                               names=names + [None])
        expected.index = exp_index

        tm.assert_frame_equal(result, expected)

        # no names
        result = concat([df, df2, df, df2],
                        keys=[('foo', 'one'), ('foo', 'two'),
                              ('baz', 'one'), ('baz', 'two')],
                        levels=levels)
        assert result.index.names == (None,) * 3

        # no levels
        result = concat([df, df2, df, df2],
                        keys=[('foo', 'one'), ('foo', 'two'),
                              ('baz', 'one'), ('baz', 'two')],
                        names=['first', 'second'])
        assert result.index.names == ('first', 'second') + (None,)
        tm.assert_index_equal(result.index.levels[0],
                              Index(['baz', 'foo'], name='first'))

    def test_concat_keys_levels_no_overlap(self):
        # GH #1406
        df = DataFrame(np.random.randn(1, 3), index=['a'])
        df2 = DataFrame(np.random.randn(1, 4), index=['b'])

        pytest.raises(ValueError, concat, [df, df],
                      keys=['one', 'two'], levels=[['foo', 'bar', 'baz']])

        pytest.raises(ValueError, concat, [df, df2],
                      keys=['one', 'two'], levels=[['foo', 'bar', 'baz']])

    def test_concat_rename_index(self):
        a = DataFrame(np.random.rand(3, 3),
                      columns=list('ABC'),
                      index=Index(list('abc'), name='index_a'))
        b = DataFrame(np.random.rand(3, 3),
                      columns=list('ABC'),
                      index=Index(list('abc'), name='index_b'))

        result = concat([a, b], keys=['key0', 'key1'],
                        names=['lvl0', 'lvl1'])

        exp = concat([a, b], keys=['key0', 'key1'], names=['lvl0'])
        names = list(exp.index.names)
        names[1] = 'lvl1'
        exp.index.set_names(names, inplace=True)

        tm.assert_frame_equal(result, exp)
        assert result.index.names == exp.index.names

    def test_crossed_dtypes_weird_corner(self):
        columns = ['A', 'B', 'C', 'D']
        df1 = DataFrame({'A': np.array([1, 2, 3, 4], dtype='f8'),
                         'B': np.array([1, 2, 3, 4], dtype='i8'),
                         'C': np.array([1, 2, 3, 4], dtype='f8'),
                         'D': np.array([1, 2, 3, 4], dtype='i8')},
                        columns=columns)

        df2 = DataFrame({'A': np.array([1, 2, 3, 4], dtype='i8'),
                         'B': np.array([1, 2, 3, 4], dtype='f8'),
                         'C': np.array([1, 2, 3, 4], dtype='i8'),
                         'D': np.array([1, 2, 3, 4], dtype='f8')},
                        columns=columns)

        appended = df1.append(df2, ignore_index=True)
        expected = DataFrame(np.concatenate([df1.values, df2.values], axis=0),
                             columns=columns)
        tm.assert_frame_equal(appended, expected)

        df = DataFrame(np.random.randn(1, 3), index=['a'])
        df2 = DataFrame(np.random.randn(1, 4), index=['b'])
        result = concat(
            [df, df2], keys=['one', 'two'], names=['first', 'second'])
        assert result.index.names == ('first', 'second')

    def test_dups_index(self):
        # GH 4771

        # single dtypes
        df = DataFrame(np.random.randint(0, 10, size=40).reshape(
            10, 4), columns=['A', 'A', 'C', 'C'])

        result = concat([df, df], axis=1)
        assert_frame_equal(result.iloc[:, :4], df)
        assert_frame_equal(result.iloc[:, 4:], df)

        result = concat([df, df], axis=0)
        assert_frame_equal(result.iloc[:10], df)
        assert_frame_equal(result.iloc[10:], df)

        # multi dtypes
        df = concat([DataFrame(np.random.randn(10, 4),
                               columns=['A', 'A', 'B', 'B']),
                     DataFrame(np.random.randint(0, 10, size=20)
                               .reshape(10, 2),
                               columns=['A', 'C'])],
                    axis=1)

        result = concat([df, df], axis=1)
        assert_frame_equal(result.iloc[:, :6], df)
        assert_frame_equal(result.iloc[:, 6:], df)

        result = concat([df, df], axis=0)
        assert_frame_equal(result.iloc[:10], df)
        assert_frame_equal(result.iloc[10:], df)

        # append
        result = df.iloc[0:8, :].append(df.iloc[8:])
        assert_frame_equal(result, df)

        result = df.iloc[0:8, :].append(df.iloc[8:9]).append(df.iloc[9:10])
        assert_frame_equal(result, df)

        expected = concat([df, df], axis=0)
        result = df.append(df)
        assert_frame_equal(result, expected)

    def test_with_mixed_tuples(self, sort):
        # 10697
        # columns have mixed tuples, so handle properly
        df1 = DataFrame({u'A': 'foo', (u'B', 1): 'bar'}, index=range(2))
        df2 = DataFrame({u'B': 'foo', (u'B', 1): 'bar'}, index=range(2))

        # it works
        concat([df1, df2], sort=sort)

    def test_handle_empty_objects(self, sort):
        df = DataFrame(np.random.randn(10, 4), columns=list('abcd'))

        baz = df[:5].copy()
        baz['foo'] = 'bar'
        empty = df[5:5]

        frames = [baz, empty, empty, df[5:]]
        concatted = concat(frames, axis=0, sort=sort)

        expected = df.reindex(columns=['a', 'b', 'c', 'd', 'foo'])
        expected['foo'] = expected['foo'].astype('O')
        expected.loc[0:4, 'foo'] = 'bar'

        tm.assert_frame_equal(concatted, expected)

        # empty as first element with time series
        # GH3259
        df = DataFrame(dict(A=range(10000)), index=date_range(
            '20130101', periods=10000, freq='s'))
        empty = DataFrame()
        result = concat([df, empty], axis=1)
        assert_frame_equal(result, df)
        result = concat([empty, df], axis=1)
        assert_frame_equal(result, df)

        result = concat([df, empty])
        assert_frame_equal(result, df)
        result = concat([empty, df])
        assert_frame_equal(result, df)

    def test_concat_mixed_objs(self):

        # concat mixed series/frames
        # G2385

        # axis 1
        index = date_range('01-Jan-2013', periods=10, freq='H')
        arr = np.arange(10, dtype='int64')
        s1 = Series(arr, index=index)
        s2 = Series(arr, index=index)
        df = DataFrame(arr.reshape(-1, 1), index=index)

        expected = DataFrame(np.repeat(arr, 2).reshape(-1, 2),
                             index=index, columns=[0, 0])
        result = concat([df, df], axis=1)
        assert_frame_equal(result, expected)

        expected = DataFrame(np.repeat(arr, 2).reshape(-1, 2),
                             index=index, columns=[0, 1])
        result = concat([s1, s2], axis=1)
        assert_frame_equal(result, expected)

        expected = DataFrame(np.repeat(arr, 3).reshape(-1, 3),
                             index=index, columns=[0, 1, 2])
        result = concat([s1, s2, s1], axis=1)
        assert_frame_equal(result, expected)

        expected = DataFrame(np.repeat(arr, 5).reshape(-1, 5),
                             index=index, columns=[0, 0, 1, 2, 3])
        result = concat([s1, df, s2, s2, s1], axis=1)
        assert_frame_equal(result, expected)

        # with names
        s1.name = 'foo'
        expected = DataFrame(np.repeat(arr, 3).reshape(-1, 3),
                             index=index, columns=['foo', 0, 0])
        result = concat([s1, df, s2], axis=1)
        assert_frame_equal(result, expected)

        s2.name = 'bar'
        expected = DataFrame(np.repeat(arr, 3).reshape(-1, 3),
                             index=index, columns=['foo', 0, 'bar'])
        result = concat([s1, df, s2], axis=1)
        assert_frame_equal(result, expected)

        # ignore index
        expected = DataFrame(np.repeat(arr, 3).reshape(-1, 3),
                             index=index, columns=[0, 1, 2])
        result = concat([s1, df, s2], axis=1, ignore_index=True)
        assert_frame_equal(result, expected)

        # axis 0
        expected = DataFrame(np.tile(arr, 3).reshape(-1, 1),
                             index=index.tolist() * 3, columns=[0])
        result = concat([s1, df, s2])
        assert_frame_equal(result, expected)

        expected = DataFrame(np.tile(arr, 3).reshape(-1, 1), columns=[0])
        result = concat([s1, df, s2], ignore_index=True)
        assert_frame_equal(result, expected)

        # invalid concatente of mixed dims
        with catch_warnings(record=True):
            panel = tm.makePanel()
            pytest.raises(ValueError, lambda: concat([panel, s1], axis=1))

    def test_empty_dtype_coerce(self):

        # xref to #12411
        # xref to #12045
        # xref to #11594
        # see below

        # 10571
        df1 = DataFrame(data=[[1, None], [2, None]], columns=['a', 'b'])
        df2 = DataFrame(data=[[3, None], [4, None]], columns=['a', 'b'])
        result = concat([df1, df2])
        expected = df1.dtypes
        tm.assert_series_equal(result.dtypes, expected)

    def test_dtype_coerceion(self):

        # 12411
        df = DataFrame({'date': [pd.Timestamp('20130101').tz_localize('UTC'),
                                 pd.NaT]})

        result = concat([df.iloc[[0]], df.iloc[[1]]])
        tm.assert_series_equal(result.dtypes, df.dtypes)

        # 12045
        import datetime
        df = DataFrame({'date': [datetime.datetime(2012, 1, 1),
                                 datetime.datetime(1012, 1, 2)]})
        result = concat([df.iloc[[0]], df.iloc[[1]]])
        tm.assert_series_equal(result.dtypes, df.dtypes)

        # 11594
        df = DataFrame({'text': ['some words'] + [None] * 9})
        result = concat([df.iloc[[0]], df.iloc[[1]]])
        tm.assert_series_equal(result.dtypes, df.dtypes)

    def test_panel_concat_other_axes(self):
        with catch_warnings(record=True):
            panel = tm.makePanel()

            p1 = panel.iloc[:, :5, :]
            p2 = panel.iloc[:, 5:, :]

            result = concat([p1, p2], axis=1)
            tm.assert_panel_equal(result, panel)

            p1 = panel.iloc[:, :, :2]
            p2 = panel.iloc[:, :, 2:]

            result = concat([p1, p2], axis=2)
            tm.assert_panel_equal(result, panel)

            # if things are a bit misbehaved
            p1 = panel.iloc[:2, :, :2]
            p2 = panel.iloc[:, :, 2:]
            p1['ItemC'] = 'baz'

            result = concat([p1, p2], axis=2)

            expected = panel.copy()
            expected['ItemC'] = expected['ItemC'].astype('O')
            expected.loc['ItemC', :, :2] = 'baz'
            tm.assert_panel_equal(result, expected)

    def test_panel_concat_buglet(self, sort):
        with catch_warnings(record=True):
            # #2257
            def make_panel():
                index = 5
                cols = 3

                def df():
                    return DataFrame(np.random.randn(index, cols),
                                     index=["I%s" % i for i in range(index)],
                                     columns=["C%s" % i for i in range(cols)])
                return Panel(dict(("Item%s" % x, df())
                                  for x in ['A', 'B', 'C']))

            panel1 = make_panel()
            panel2 = make_panel()

            panel2 = panel2.rename_axis(dict((x, "%s_1" % x)
                                             for x in panel2.major_axis),
                                        axis=1)

            panel3 = panel2.rename_axis(lambda x: '%s_1' % x, axis=1)
            panel3 = panel3.rename_axis(lambda x: '%s_1' % x, axis=2)

            # it works!
            concat([panel1, panel3], axis=1, verify_integrity=True, sort=sort)

    def test_concat_series(self):

        ts = tm.makeTimeSeries()
        ts.name = 'foo'

        pieces = [ts[:5], ts[5:15], ts[15:]]

        result = concat(pieces)
        tm.assert_series_equal(result, ts)
        assert result.name == ts.name

        result = concat(pieces, keys=[0, 1, 2])
        expected = ts.copy()

        ts.index = DatetimeIndex(np.array(ts.index.values, dtype='M8[ns]'))

        exp_labels = [np.repeat([0, 1, 2], [len(x) for x in pieces]),
                      np.arange(len(ts))]
        exp_index = MultiIndex(levels=[[0, 1, 2], ts.index],
                               labels=exp_labels)
        expected.index = exp_index
        tm.assert_series_equal(result, expected)

    def test_concat_series_axis1(self, sort=sort):
        ts = tm.makeTimeSeries()

        pieces = [ts[:-2], ts[2:], ts[2:-2]]

        result = concat(pieces, axis=1)
        expected = DataFrame(pieces).T
        assert_frame_equal(result, expected)

        result = concat(pieces, keys=['A', 'B', 'C'], axis=1)
        expected = DataFrame(pieces, index=['A', 'B', 'C']).T
        assert_frame_equal(result, expected)

        # preserve series names, #2489
        s = Series(randn(5), name='A')
        s2 = Series(randn(5), name='B')

        result = concat([s, s2], axis=1)
        expected = DataFrame({'A': s, 'B': s2})
        assert_frame_equal(result, expected)

        s2.name = None
        result = concat([s, s2], axis=1)
        tm.assert_index_equal(result.columns,
                              Index(['A', 0], dtype='object'))

        # must reindex, #2603
        s = Series(randn(3), index=['c', 'a', 'b'], name='A')
        s2 = Series(randn(4), index=['d', 'a', 'b', 'c'], name='B')
        result = concat([s, s2], axis=1, sort=sort)
        expected = DataFrame({'A': s, 'B': s2})
        assert_frame_equal(result, expected)

    def test_concat_single_with_key(self):
        df = DataFrame(np.random.randn(10, 4))

        result = concat([df], keys=['foo'])
        expected = concat([df, df], keys=['foo', 'bar'])
        tm.assert_frame_equal(result, expected[:10])

    def test_concat_exclude_none(self):
        df = DataFrame(np.random.randn(10, 4))

        pieces = [df[:5], None, None, df[5:]]
        result = concat(pieces)
        tm.assert_frame_equal(result, df)
        pytest.raises(ValueError, concat, [None, None])

    def test_concat_datetime64_block(self):
        from pandas.core.indexes.datetimes import date_range

        rng = date_range('1/1/2000', periods=10)

        df = DataFrame({'time': rng})

        result = concat([df, df])
        assert (result.iloc[:10]['time'] == rng).all()
        assert (result.iloc[10:]['time'] == rng).all()

    def test_concat_timedelta64_block(self):
        from pandas import to_timedelta

        rng = to_timedelta(np.arange(10), unit='s')

        df = DataFrame({'time': rng})

        result = concat([df, df])
        assert (result.iloc[:10]['time'] == rng).all()
        assert (result.iloc[10:]['time'] == rng).all()

    def test_concat_keys_with_none(self):
        # #1649
        df0 = DataFrame([[10, 20, 30], [10, 20, 30], [10, 20, 30]])

        result = concat(dict(a=None, b=df0, c=df0[:2], d=df0[:1], e=df0))
        expected = concat(dict(b=df0, c=df0[:2], d=df0[:1], e=df0))
        tm.assert_frame_equal(result, expected)

        result = concat([None, df0, df0[:2], df0[:1], df0],
                        keys=['a', 'b', 'c', 'd', 'e'])
        expected = concat([df0, df0[:2], df0[:1], df0],
                          keys=['b', 'c', 'd', 'e'])
        tm.assert_frame_equal(result, expected)

    def test_concat_bug_1719(self):
        ts1 = tm.makeTimeSeries()
        ts2 = tm.makeTimeSeries()[::2]

        # to join with union
        # these two are of different length!
        left = concat([ts1, ts2], join='outer', axis=1)
        right = concat([ts2, ts1], join='outer', axis=1)

        assert len(left) == len(right)

    def test_concat_bug_2972(self):
        ts0 = Series(np.zeros(5))
        ts1 = Series(np.ones(5))
        ts0.name = ts1.name = 'same name'
        result = concat([ts0, ts1], axis=1)

        expected = DataFrame({0: ts0, 1: ts1})
        expected.columns = ['same name', 'same name']
        assert_frame_equal(result, expected)

    def test_concat_bug_3602(self):

        # GH 3602, duplicate columns
        df1 = DataFrame({'firmNo': [0, 0, 0, 0], 'prc': [6, 6, 6, 6],
                         'stringvar': ['rrr', 'rrr', 'rrr', 'rrr']})
        df2 = DataFrame({'C': [9, 10, 11, 12], 'misc': [1, 2, 3, 4],
                         'prc': [6, 6, 6, 6]})
        expected = DataFrame([[0, 6, 'rrr', 9, 1, 6],
                              [0, 6, 'rrr', 10, 2, 6],
                              [0, 6, 'rrr', 11, 3, 6],
                              [0, 6, 'rrr', 12, 4, 6]])
        expected.columns = ['firmNo', 'prc', 'stringvar', 'C', 'misc', 'prc']

        result = concat([df1, df2], axis=1)
        assert_frame_equal(result, expected)

    def test_concat_inner_join_empty(self):
        # GH 15328
        df_empty = pd.DataFrame()
        df_a = pd.DataFrame({'a': [1, 2]}, index=[0, 1], dtype='int64')
        df_expected = pd.DataFrame({'a': []}, index=[], dtype='int64')

        for how, expected in [('inner', df_expected), ('outer', df_a)]:
            result = pd.concat([df_a, df_empty], axis=1, join=how)
            assert_frame_equal(result, expected)

    def test_concat_series_axis1_same_names_ignore_index(self):
        dates = date_range('01-Jan-2013', '01-Jan-2014', freq='MS')[0:-1]
        s1 = Series(randn(len(dates)), index=dates, name='value')
        s2 = Series(randn(len(dates)), index=dates, name='value')

        result = concat([s1, s2], axis=1, ignore_index=True)
        expected = Index([0, 1])

        tm.assert_index_equal(result.columns, expected)

    def test_concat_iterables(self):
        from collections import deque, Iterable

        # GH8645 check concat works with tuples, list, generators, and weird
        # stuff like deque and custom iterables
        df1 = DataFrame([1, 2, 3])
        df2 = DataFrame([4, 5, 6])
        expected = DataFrame([1, 2, 3, 4, 5, 6])
        assert_frame_equal(concat((df1, df2), ignore_index=True), expected)
        assert_frame_equal(concat([df1, df2], ignore_index=True), expected)
        assert_frame_equal(concat((df for df in (df1, df2)),
                                  ignore_index=True), expected)
        assert_frame_equal(
            concat(deque((df1, df2)), ignore_index=True), expected)

        class CustomIterator1(object):

            def __len__(self):
                return 2

            def __getitem__(self, index):
                try:
                    return {0: df1, 1: df2}[index]
                except KeyError:
                    raise IndexError
        assert_frame_equal(pd.concat(CustomIterator1(),
                                     ignore_index=True), expected)

        class CustomIterator2(Iterable):

            def __iter__(self):
                yield df1
                yield df2
        assert_frame_equal(pd.concat(CustomIterator2(),
                                     ignore_index=True), expected)

    def test_concat_invalid(self):

        # trying to concat a ndframe with a non-ndframe
        df1 = mkdf(10, 2)
        for obj in [1, dict(), [1, 2], (1, 2)]:
            pytest.raises(TypeError, lambda x: concat([df1, obj]))

    def test_concat_invalid_first_argument(self):
        df1 = mkdf(10, 2)
        df2 = mkdf(10, 2)
        pytest.raises(TypeError, concat, df1, df2)

        # generator ok though
        concat(DataFrame(np.random.rand(5, 5)) for _ in range(3))

        # text reader ok
        # GH6583
        data = """index,A,B,C,D
foo,2,3,4,5
bar,7,8,9,10
baz,12,13,14,15
qux,12,13,14,15
foo2,12,13,14,15
bar2,12,13,14,15
"""

        reader = read_csv(StringIO(data), chunksize=1)
        result = concat(reader, ignore_index=True)
        expected = read_csv(StringIO(data))
        assert_frame_equal(result, expected)

    def test_concat_NaT_series(self):
        # GH 11693
        # test for merging NaT series with datetime series.
        x = Series(date_range('20151124 08:00', '20151124 09:00',
                              freq='1h', tz='US/Eastern'))
        y = Series(pd.NaT, index=[0, 1], dtype='datetime64[ns, US/Eastern]')
        expected = Series([x[0], x[1], pd.NaT, pd.NaT])

        result = concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)

        # all NaT with tz
        expected = Series(pd.NaT, index=range(4),
                          dtype='datetime64[ns, US/Eastern]')
        result = pd.concat([y, y], ignore_index=True)
        tm.assert_series_equal(result, expected)

        # without tz
        x = pd.Series(pd.date_range('20151124 08:00',
                                    '20151124 09:00', freq='1h'))
        y = pd.Series(pd.date_range('20151124 10:00',
                                    '20151124 11:00', freq='1h'))
        y[:] = pd.NaT
        expected = pd.Series([x[0], x[1], pd.NaT, pd.NaT])
        result = pd.concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)

        # all NaT without tz
        x[:] = pd.NaT
        expected = pd.Series(pd.NaT, index=range(4),
                             dtype='datetime64[ns]')
        result = pd.concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)

    def test_concat_tz_frame(self):
        df2 = DataFrame(dict(A=pd.Timestamp('20130102', tz='US/Eastern'),
                             B=pd.Timestamp('20130603', tz='CET')),
                        index=range(5))

        # concat
        df3 = pd.concat([df2.A.to_frame(), df2.B.to_frame()], axis=1)
        assert_frame_equal(df2, df3)

    def test_concat_tz_series(self):
        # gh-11755: tz and no tz
        x = Series(date_range('20151124 08:00',
                              '20151124 09:00',
                              freq='1h', tz='UTC'))
        y = Series(date_range('2012-01-01', '2012-01-02'))
        expected = Series([x[0], x[1], y[0], y[1]],
                          dtype='object')
        result = concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)

        # gh-11887: concat tz and object
        x = Series(date_range('20151124 08:00',
                              '20151124 09:00',
                              freq='1h', tz='UTC'))
        y = Series(['a', 'b'])
        expected = Series([x[0], x[1], y[0], y[1]],
                          dtype='object')
        result = concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)

        # see gh-12217 and gh-12306
        # Concatenating two UTC times
        first = pd.DataFrame([[datetime(2016, 1, 1)]])
        first[0] = first[0].dt.tz_localize('UTC')

        second = pd.DataFrame([[datetime(2016, 1, 2)]])
        second[0] = second[0].dt.tz_localize('UTC')

        result = pd.concat([first, second])
        assert result[0].dtype == 'datetime64[ns, UTC]'

        # Concatenating two London times
        first = pd.DataFrame([[datetime(2016, 1, 1)]])
        first[0] = first[0].dt.tz_localize('Europe/London')

        second = pd.DataFrame([[datetime(2016, 1, 2)]])
        second[0] = second[0].dt.tz_localize('Europe/London')

        result = pd.concat([first, second])
        assert result[0].dtype == 'datetime64[ns, Europe/London]'

        # Concatenating 2+1 London times
        first = pd.DataFrame([[datetime(2016, 1, 1)], [datetime(2016, 1, 2)]])
        first[0] = first[0].dt.tz_localize('Europe/London')

        second = pd.DataFrame([[datetime(2016, 1, 3)]])
        second[0] = second[0].dt.tz_localize('Europe/London')

        result = pd.concat([first, second])
        assert result[0].dtype == 'datetime64[ns, Europe/London]'

        # Concat'ing 1+2 London times
        first = pd.DataFrame([[datetime(2016, 1, 1)]])
        first[0] = first[0].dt.tz_localize('Europe/London')

        second = pd.DataFrame([[datetime(2016, 1, 2)], [datetime(2016, 1, 3)]])
        second[0] = second[0].dt.tz_localize('Europe/London')

        result = pd.concat([first, second])
        assert result[0].dtype == 'datetime64[ns, Europe/London]'

    def test_concat_tz_series_with_datetimelike(self):
        # see gh-12620: tz and timedelta
        x = [pd.Timestamp('2011-01-01', tz='US/Eastern'),
             pd.Timestamp('2011-02-01', tz='US/Eastern')]
        y = [pd.Timedelta('1 day'), pd.Timedelta('2 day')]
        result = concat([pd.Series(x), pd.Series(y)], ignore_index=True)
        tm.assert_series_equal(result, pd.Series(x + y, dtype='object'))

        # tz and period
        y = [pd.Period('2011-03', freq='M'), pd.Period('2011-04', freq='M')]
        result = concat([pd.Series(x), pd.Series(y)], ignore_index=True)
        tm.assert_series_equal(result, pd.Series(x + y, dtype='object'))

    def test_concat_tz_series_tzlocal(self):
        # see gh-13583
        x = [pd.Timestamp('2011-01-01', tz=dateutil.tz.tzlocal()),
             pd.Timestamp('2011-02-01', tz=dateutil.tz.tzlocal())]
        y = [pd.Timestamp('2012-01-01', tz=dateutil.tz.tzlocal()),
             pd.Timestamp('2012-02-01', tz=dateutil.tz.tzlocal())]

        result = concat([pd.Series(x), pd.Series(y)], ignore_index=True)
        tm.assert_series_equal(result, pd.Series(x + y))
        assert result.dtype == 'datetime64[ns, tzlocal()]'

    @pytest.mark.parametrize('tz1', [None, 'UTC'])
    @pytest.mark.parametrize('tz2', [None, 'UTC'])
    @pytest.mark.parametrize('s', [pd.NaT, pd.Timestamp('20150101')])
    def test_concat_NaT_dataframes_all_NaT_axis_0(self, tz1, tz2, s):
        # GH 12396

        # tz-naive
        first = pd.DataFrame([[pd.NaT], [pd.NaT]]).apply(
            lambda x: x.dt.tz_localize(tz1))
        second = pd.DataFrame([s]).apply(lambda x: x.dt.tz_localize(tz2))

        result = pd.concat([first, second], axis=0)
        expected = pd.DataFrame(pd.Series(
            [pd.NaT, pd.NaT, s], index=[0, 1, 0]))
        expected = expected.apply(lambda x: x.dt.tz_localize(tz2))
        if tz1 != tz2:
            expected = expected.astype(object)

        assert_frame_equal(result, expected)

    @pytest.mark.parametrize('tz1', [None, 'UTC'])
    @pytest.mark.parametrize('tz2', [None, 'UTC'])
    def test_concat_NaT_dataframes_all_NaT_axis_1(self, tz1, tz2):
        # GH 12396

        first = pd.DataFrame(pd.Series([pd.NaT, pd.NaT]).dt.tz_localize(tz1))
        second = pd.DataFrame(pd.Series(
            [pd.NaT]).dt.tz_localize(tz2), columns=[1])
        expected = pd.DataFrame(
            {0: pd.Series([pd.NaT, pd.NaT]).dt.tz_localize(tz1),
             1: pd.Series([pd.NaT, pd.NaT]).dt.tz_localize(tz2)}
        )
        result = pd.concat([first, second], axis=1)
        assert_frame_equal(result, expected)

    @pytest.mark.parametrize('tz1', [None, 'UTC'])
    @pytest.mark.parametrize('tz2', [None, 'UTC'])
    def test_concat_NaT_series_dataframe_all_NaT(self, tz1, tz2):
        # GH 12396

        # tz-naive
        first = pd.Series([pd.NaT, pd.NaT]).dt.tz_localize(tz1)
        second = pd.DataFrame([[pd.Timestamp('2015/01/01', tz=tz2)],
                               [pd.Timestamp('2016/01/01', tz=tz2)]],
                              index=[2, 3])

        expected = pd.DataFrame([pd.NaT, pd.NaT,
                                 pd.Timestamp('2015/01/01', tz=tz2),
                                 pd.Timestamp('2016/01/01', tz=tz2)])
        if tz1 != tz2:
            expected = expected.astype(object)

        result = pd.concat([first, second])
        assert_frame_equal(result, expected)

    @pytest.mark.parametrize('tz', [None, 'UTC'])
    def test_concat_NaT_dataframes(self, tz):
        # GH 12396

        first = pd.DataFrame([[pd.NaT], [pd.NaT]])
        first = first.apply(lambda x: x.dt.tz_localize(tz))
        second = pd.DataFrame([[pd.Timestamp('2015/01/01', tz=tz)],
                               [pd.Timestamp('2016/01/01', tz=tz)]],
                              index=[2, 3])
        expected = pd.DataFrame([pd.NaT, pd.NaT,
                                 pd.Timestamp('2015/01/01', tz=tz),
                                 pd.Timestamp('2016/01/01', tz=tz)])

        result = pd.concat([first, second], axis=0)
        assert_frame_equal(result, expected)

    def test_concat_period_series(self):
        x = Series(pd.PeriodIndex(['2015-11-01', '2015-12-01'], freq='D'))
        y = Series(pd.PeriodIndex(['2015-10-01', '2016-01-01'], freq='D'))
        expected = Series([x[0], x[1], y[0], y[1]], dtype='object')
        result = concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)
        assert result.dtype == 'object'

        # different freq
        x = Series(pd.PeriodIndex(['2015-11-01', '2015-12-01'], freq='D'))
        y = Series(pd.PeriodIndex(['2015-10-01', '2016-01-01'], freq='M'))
        expected = Series([x[0], x[1], y[0], y[1]], dtype='object')
        result = concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)
        assert result.dtype == 'object'

        x = Series(pd.PeriodIndex(['2015-11-01', '2015-12-01'], freq='D'))
        y = Series(pd.PeriodIndex(['2015-11-01', '2015-12-01'], freq='M'))
        expected = Series([x[0], x[1], y[0], y[1]], dtype='object')
        result = concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)
        assert result.dtype == 'object'

        # non-period
        x = Series(pd.PeriodIndex(['2015-11-01', '2015-12-01'], freq='D'))
        y = Series(pd.DatetimeIndex(['2015-11-01', '2015-12-01']))
        expected = Series([x[0], x[1], y[0], y[1]], dtype='object')
        result = concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)
        assert result.dtype == 'object'

        x = Series(pd.PeriodIndex(['2015-11-01', '2015-12-01'], freq='D'))
        y = Series(['A', 'B'])
        expected = Series([x[0], x[1], y[0], y[1]], dtype='object')
        result = concat([x, y], ignore_index=True)
        tm.assert_series_equal(result, expected)
        assert result.dtype == 'object'

    def test_concat_empty_series(self):
        # GH 11082
        s1 = pd.Series([1, 2, 3], name='x')
        s2 = pd.Series(name='y')
        res = pd.concat([s1, s2], axis=1)
        exp = pd.DataFrame({'x': [1, 2, 3], 'y': [np.nan, np.nan, np.nan]})
        tm.assert_frame_equal(res, exp)

        s1 = pd.Series([1, 2, 3], name='x')
        s2 = pd.Series(name='y')
        res = pd.concat([s1, s2], axis=0)
        # name will be reset
        exp = pd.Series([1, 2, 3])
        tm.assert_series_equal(res, exp)

        # empty Series with no name
        s1 = pd.Series([1, 2, 3], name='x')
        s2 = pd.Series(name=None)
        res = pd.concat([s1, s2], axis=1)
        exp = pd.DataFrame({'x': [1, 2, 3], 0: [np.nan, np.nan, np.nan]},
                           columns=['x', 0])
        tm.assert_frame_equal(res, exp)

    @pytest.mark.parametrize('tz', [None, 'UTC'])
    @pytest.mark.parametrize('values', [[], [1, 2, 3]])
    def test_concat_empty_series_timelike(self, tz, values):
        # GH 18447

        first = Series([], dtype='M8[ns]').dt.tz_localize(tz)
        second = Series(values)
        expected = DataFrame(
            {0: pd.Series([pd.NaT] * len(values),
                          dtype='M8[ns]'
                          ).dt.tz_localize(tz),
             1: values})
        result = concat([first, second], axis=1)
        assert_frame_equal(result, expected)

    def test_default_index(self):
        # is_series and ignore_index
        s1 = pd.Series([1, 2, 3], name='x')
        s2 = pd.Series([4, 5, 6], name='y')
        res = pd.concat([s1, s2], axis=1, ignore_index=True)
        assert isinstance(res.columns, pd.RangeIndex)
        exp = pd.DataFrame([[1, 4], [2, 5], [3, 6]])
        # use check_index_type=True to check the result have
        # RangeIndex (default index)
        tm.assert_frame_equal(res, exp, check_index_type=True,
                              check_column_type=True)

        # is_series and all inputs have no names
        s1 = pd.Series([1, 2, 3])
        s2 = pd.Series([4, 5, 6])
        res = pd.concat([s1, s2], axis=1, ignore_index=False)
        assert isinstance(res.columns, pd.RangeIndex)
        exp = pd.DataFrame([[1, 4], [2, 5], [3, 6]])
        exp.columns = pd.RangeIndex(2)
        tm.assert_frame_equal(res, exp, check_index_type=True,
                              check_column_type=True)

        # is_dataframe and ignore_index
        df1 = pd.DataFrame({'A': [1, 2], 'B': [5, 6]})
        df2 = pd.DataFrame({'A': [3, 4], 'B': [7, 8]})

        res = pd.concat([df1, df2], axis=0, ignore_index=True)
        exp = pd.DataFrame([[1, 5], [2, 6], [3, 7], [4, 8]],
                           columns=['A', 'B'])
        tm.assert_frame_equal(res, exp, check_index_type=True,
                              check_column_type=True)

        res = pd.concat([df1, df2], axis=1, ignore_index=True)
        exp = pd.DataFrame([[1, 5, 3, 7], [2, 6, 4, 8]])
        tm.assert_frame_equal(res, exp, check_index_type=True,
                              check_column_type=True)

    def test_concat_multiindex_rangeindex(self):
        # GH13542
        # when multi-index levels are RangeIndex objects
        # there is a bug in concat with objects of len 1

        df = DataFrame(np.random.randn(9, 2))
        df.index = MultiIndex(levels=[pd.RangeIndex(3), pd.RangeIndex(3)],
                              labels=[np.repeat(np.arange(3), 3),
                                      np.tile(np.arange(3), 3)])

        res = concat([df.iloc[[2, 3, 4], :], df.iloc[[5], :]])
        exp = df.iloc[[2, 3, 4, 5], :]
        tm.assert_frame_equal(res, exp)

    def test_concat_multiindex_dfs_with_deepcopy(self):
        # GH 9967
        from copy import deepcopy
        example_multiindex1 = pd.MultiIndex.from_product([['a'], ['b']])
        example_dataframe1 = pd.DataFrame([0], index=example_multiindex1)

        example_multiindex2 = pd.MultiIndex.from_product([['a'], ['c']])
        example_dataframe2 = pd.DataFrame([1], index=example_multiindex2)

        example_dict = {'s1': example_dataframe1, 's2': example_dataframe2}
        expected_index = pd.MultiIndex(levels=[['s1', 's2'],
                                               ['a'],
                                               ['b', 'c']],
                                       labels=[[0, 1], [0, 0], [0, 1]],
                                       names=['testname', None, None])
        expected = pd.DataFrame([[0], [1]], index=expected_index)
        result_copy = pd.concat(deepcopy(example_dict), names=['testname'])
        tm.assert_frame_equal(result_copy, expected)
        result_no_copy = pd.concat(example_dict, names=['testname'])
        tm.assert_frame_equal(result_no_copy, expected)

    def test_categorical_concat_append(self):
        cat = Categorical(["a", "b"], categories=["a", "b"])
        vals = [1, 2]
        df = DataFrame({"cats": cat, "vals": vals})
        cat2 = Categorical(["a", "b", "a", "b"], categories=["a", "b"])
        vals2 = [1, 2, 1, 2]
        exp = DataFrame({"cats": cat2, "vals": vals2},
                        index=Index([0, 1, 0, 1]))

        tm.assert_frame_equal(pd.concat([df, df]), exp)
        tm.assert_frame_equal(df.append(df), exp)

        # GH 13524 can concat different categories
        cat3 = Categorical(["a", "b"], categories=["a", "b", "c"])
        vals3 = [1, 2]
        df_different_categories = DataFrame({"cats": cat3, "vals": vals3})

        res = pd.concat([df, df_different_categories], ignore_index=True)
        exp = DataFrame({"cats": list('abab'), "vals": [1, 2, 1, 2]})
        tm.assert_frame_equal(res, exp)

        res = df.append(df_different_categories, ignore_index=True)
        tm.assert_frame_equal(res, exp)

    def test_categorical_concat_dtypes(self):

        # GH8143
        index = ['cat', 'obj', 'num']
        cat = Categorical(['a', 'b', 'c'])
        obj = Series(['a', 'b', 'c'])
        num = Series([1, 2, 3])
        df = pd.concat([Series(cat), obj, num], axis=1, keys=index)

        result = df.dtypes == 'object'
        expected = Series([False, True, False], index=index)
        tm.assert_series_equal(result, expected)

        result = df.dtypes == 'int64'
        expected = Series([False, False, True], index=index)
        tm.assert_series_equal(result, expected)

        result = df.dtypes == 'category'
        expected = Series([True, False, False], index=index)
        tm.assert_series_equal(result, expected)

    def test_categorical_concat(self, sort):
        # See GH 10177
        df1 = DataFrame(np.arange(18, dtype='int64').reshape(6, 3),
                        columns=["a", "b", "c"])

        df2 = DataFrame(np.arange(14, dtype='int64').reshape(7, 2),
                        columns=["a", "c"])

        cat_values = ["one", "one", "two", "one", "two", "two", "one"]
        df2['h'] = Series(Categorical(cat_values))

        res = pd.concat((df1, df2), axis=0, ignore_index=True, sort=sort)
        exp = DataFrame({'a': [0, 3, 6, 9, 12, 15, 0, 2, 4, 6, 8, 10, 12],
                         'b': [1, 4, 7, 10, 13, 16, np.nan, np.nan, np.nan,
                               np.nan, np.nan, np.nan, np.nan],
                         'c': [2, 5, 8, 11, 14, 17, 1, 3, 5, 7, 9, 11, 13],
                         'h': [None] * 6 + cat_values})
        tm.assert_frame_equal(res, exp)

    def test_categorical_concat_gh7864(self):
        # GH 7864
        # make sure ordering is preserverd
        df = DataFrame({"id": [1, 2, 3, 4, 5, 6], "raw_grade": list('abbaae')})
        df["grade"] = Categorical(df["raw_grade"])
        df['grade'].cat.set_categories(['e', 'a', 'b'])

        df1 = df[0:3]
        df2 = df[3:]

        tm.assert_index_equal(df['grade'].cat.categories,
                              df1['grade'].cat.categories)
        tm.assert_index_equal(df['grade'].cat.categories,
                              df2['grade'].cat.categories)

        dfx = pd.concat([df1, df2])
        tm.assert_index_equal(df['grade'].cat.categories,
                              dfx['grade'].cat.categories)

        dfa = df1.append(df2)
        tm.assert_index_equal(df['grade'].cat.categories,
                              dfa['grade'].cat.categories)

    def test_categorical_concat_preserve(self):

        # GH 8641  series concat not preserving category dtype
        # GH 13524 can concat different categories
        s = Series(list('abc'), dtype='category')
        s2 = Series(list('abd'), dtype='category')

        exp = Series(list('abcabd'))
        res = pd.concat([s, s2], ignore_index=True)
        tm.assert_series_equal(res, exp)

        exp = Series(list('abcabc'), dtype='category')
        res = pd.concat([s, s], ignore_index=True)
        tm.assert_series_equal(res, exp)

        exp = Series(list('abcabc'), index=[0, 1, 2, 0, 1, 2],
                     dtype='category')
        res = pd.concat([s, s])
        tm.assert_series_equal(res, exp)

        a = Series(np.arange(6, dtype='int64'))
        b = Series(list('aabbca'))

        df2 = DataFrame({'A': a,
                         'B': b.astype(CategoricalDtype(list('cab')))})
        res = pd.concat([df2, df2])
        exp = DataFrame(
            {'A': pd.concat([a, a]),
             'B': pd.concat([b, b]).astype(CategoricalDtype(list('cab')))})
        tm.assert_frame_equal(res, exp)

    def test_categorical_index_preserver(self):

        a = Series(np.arange(6, dtype='int64'))
        b = Series(list('aabbca'))

        df2 = DataFrame({'A': a,
                         'B': b.astype(CategoricalDtype(list('cab')))
                         }).set_index('B')
        result = pd.concat([df2, df2])
        expected = DataFrame(
            {'A': pd.concat([a, a]),
             'B': pd.concat([b, b]).astype(CategoricalDtype(list('cab')))
             }).set_index('B')
        tm.assert_frame_equal(result, expected)

        # wrong catgories
        df3 = DataFrame({'A': a, 'B': Categorical(b, categories=list('abe'))
                         }).set_index('B')
        pytest.raises(TypeError, lambda: pd.concat([df2, df3]))

    def test_concat_categoricalindex(self):
        # GH 16111, categories that aren't lexsorted
        categories = [9, 0, 1, 2, 3]

        a = pd.Series(1, index=pd.CategoricalIndex([9, 0],
                                                   categories=categories))
        b = pd.Series(2, index=pd.CategoricalIndex([0, 1],
                                                   categories=categories))
        c = pd.Series(3, index=pd.CategoricalIndex([1, 2],
                                                   categories=categories))

        result = pd.concat([a, b, c], axis=1)

        exp_idx = pd.CategoricalIndex([0, 1, 2, 9])
        exp = pd.DataFrame({0: [1, np.nan, np.nan, 1],
                            1: [2, 2, np.nan, np.nan],
                            2: [np.nan, 3, 3, np.nan]},
                           columns=[0, 1, 2],
                           index=exp_idx)
        tm.assert_frame_equal(result, exp)

    def test_concat_order(self):
        # GH 17344
        dfs = [pd.DataFrame(index=range(3), columns=['a', 1, None])]
        dfs += [pd.DataFrame(index=range(3), columns=[None, 1, 'a'])
                for i in range(100)]

        result = pd.concat(dfs, sort=True).columns

        if PY2:
            # Different sort order between incomparable objects between
            # python 2 and python3 via Index.union.
            expected = dfs[1].columns
        else:
            expected = dfs[0].columns
        tm.assert_index_equal(result, expected)

    def test_concat_datetime_timezone(self):
        # GH 18523
        idx1 = pd.date_range('2011-01-01', periods=3, freq='H',
                             tz='Europe/Paris')
        idx2 = pd.date_range(start=idx1[0], end=idx1[-1], freq='H')
        df1 = pd.DataFrame({'a': [1, 2, 3]}, index=idx1)
        df2 = pd.DataFrame({'b': [1, 2, 3]}, index=idx2)
        result = pd.concat([df1, df2], axis=1)

        exp_idx = DatetimeIndex(['2011-01-01 00:00:00+01:00',
                                 '2011-01-01 01:00:00+01:00',
                                 '2011-01-01 02:00:00+01:00'],
                                freq='H'
                                ).tz_localize('UTC').tz_convert('Europe/Paris')

        expected = pd.DataFrame([[1, 1], [2, 2], [3, 3]],
                                index=exp_idx, columns=['a', 'b'])

        tm.assert_frame_equal(result, expected)

        idx3 = pd.date_range('2011-01-01', periods=3,
                             freq='H', tz='Asia/Tokyo')
        df3 = pd.DataFrame({'b': [1, 2, 3]}, index=idx3)
        result = pd.concat([df1, df3], axis=1)

        exp_idx = DatetimeIndex(['2010-12-31 15:00:00+00:00',
                                 '2010-12-31 16:00:00+00:00',
                                 '2010-12-31 17:00:00+00:00',
                                 '2010-12-31 23:00:00+00:00',
                                 '2011-01-01 00:00:00+00:00',
                                 '2011-01-01 01:00:00+00:00']
                                ).tz_localize('UTC')

        expected = pd.DataFrame([[np.nan, 1], [np.nan, 2], [np.nan, 3],
                                 [1, np.nan], [2, np.nan], [3, np.nan]],
                                index=exp_idx, columns=['a', 'b'])

        tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize('pdt', [pd.Series, pd.DataFrame, pd.Panel])
@pytest.mark.parametrize('dt', np.sctypes['float'])
def test_concat_no_unnecessary_upcast(dt, pdt):
    with catch_warnings(record=True):
        # GH 13247
        dims = pdt().ndim
        dfs = [pdt(np.array([1], dtype=dt, ndmin=dims)),
               pdt(np.array([np.nan], dtype=dt, ndmin=dims)),
               pdt(np.array([5], dtype=dt, ndmin=dims))]
        x = pd.concat(dfs)
        assert x.values.dtype == dt


@pytest.mark.parametrize('pdt', [pd.Series, pd.DataFrame, pd.Panel])
@pytest.mark.parametrize('dt', np.sctypes['int'])
def test_concat_will_upcast(dt, pdt):
    with catch_warnings(record=True):
        dims = pdt().ndim
        dfs = [pdt(np.array([1], dtype=dt, ndmin=dims)),
               pdt(np.array([np.nan], ndmin=dims)),
               pdt(np.array([5], dtype=dt, ndmin=dims))]
        x = pd.concat(dfs)
        assert x.values.dtype == 'float64'


def test_concat_empty_and_non_empty_frame_regression():
    # GH 18178 regression test
    df1 = pd.DataFrame({'foo': [1]})
    df2 = pd.DataFrame({'foo': []})
    expected = pd.DataFrame({'foo': [1.0]})
    result = pd.concat([df1, df2])
    assert_frame_equal(result, expected)


def test_concat_empty_and_non_empty_series_regression():
    # GH 18187 regression test
    s1 = pd.Series([1])
    s2 = pd.Series([])
    expected = s1
    result = pd.concat([s1, s2])
    tm.assert_series_equal(result, expected)


def test_concat_sorts_columns(sort_with_none):
    # GH-4588
    df1 = pd.DataFrame({"a": [1, 2], "b": [1, 2]}, columns=['b', 'a'])
    df2 = pd.DataFrame({"a": [3, 4], "c": [5, 6]})

    # for sort=True/None
    expected = pd.DataFrame({"a": [1, 2, 3, 4],
                             "b": [1, 2, None, None],
                             "c": [None, None, 5, 6]},
                            columns=['a', 'b', 'c'])

    if sort_with_none is False:
        expected = expected[['b', 'a', 'c']]

    if sort_with_none is None:
        # only warn if not explicitly specified
        ctx = tm.assert_produces_warning(FutureWarning)
    else:
        ctx = tm.assert_produces_warning(None)

    # default
    with ctx:
        result = pd.concat([df1, df2], ignore_index=True, sort=sort_with_none)
    tm.assert_frame_equal(result, expected)


def test_concat_sorts_index(sort_with_none):
    df1 = pd.DataFrame({"a": [1, 2, 3]}, index=['c', 'a', 'b'])
    df2 = pd.DataFrame({"b": [1, 2]}, index=['a', 'b'])

    # For True/None
    expected = pd.DataFrame({"a": [2, 3, 1], "b": [1, 2, None]},
                            index=['a', 'b', 'c'],
                            columns=['a', 'b'])
    if sort_with_none is False:
        expected = expected.loc[['c', 'a', 'b']]

    if sort_with_none is None:
        # only warn if not explicitly specified
        ctx = tm.assert_produces_warning(FutureWarning)
    else:
        ctx = tm.assert_produces_warning(None)

    # Warn and sort by default
    with ctx:
        result = pd.concat([df1, df2], axis=1, sort=sort_with_none)
    tm.assert_frame_equal(result, expected)


def test_concat_inner_sort(sort_with_none):
    # https://github.com/pandas-dev/pandas/pull/20613
    df1 = pd.DataFrame({"a": [1, 2], "b": [1, 2], "c": [1, 2]},
                       columns=['b', 'a', 'c'])
    df2 = pd.DataFrame({"a": [1, 2], 'b': [3, 4]}, index=[3, 4])

    with tm.assert_produces_warning(None):
        # unset sort should *not* warn for inner join
        # since that never sorted
        result = pd.concat([df1, df2], sort=sort_with_none,
                           join='inner',
                           ignore_index=True)

    expected = pd.DataFrame({"b": [1, 2, 3, 4], "a": [1, 2, 1, 2]},
                            columns=['b', 'a'])
    if sort_with_none is True:
        expected = expected[['a', 'b']]
    tm.assert_frame_equal(result, expected)


def test_concat_aligned_sort():
    # GH-4588
    df = pd.DataFrame({"c": [1, 2], "b": [3, 4], 'a': [5, 6]},
                      columns=['c', 'b', 'a'])
    result = pd.concat([df, df], sort=True, ignore_index=True)
    expected = pd.DataFrame({'a': [5, 6, 5, 6], 'b': [3, 4, 3, 4],
                             'c': [1, 2, 1, 2]},
                            columns=['a', 'b', 'c'])
    tm.assert_frame_equal(result, expected)

    result = pd.concat([df, df[['c', 'b']]], join='inner', sort=True,
                       ignore_index=True)
    expected = expected[['b', 'c']]
    tm.assert_frame_equal(result, expected)


def test_concat_aligned_sort_does_not_raise():
    # GH-4588
    # We catch TypeErrors from sorting internally and do not re-raise.
    df = pd.DataFrame({1: [1, 2], "a": [3, 4]}, columns=[1, 'a'])
    expected = pd.DataFrame({1: [1, 2, 1, 2], 'a': [3, 4, 3, 4]},
                            columns=[1, 'a'])
    result = pd.concat([df, df], ignore_index=True, sort=True)
    tm.assert_frame_equal(result, expected)
