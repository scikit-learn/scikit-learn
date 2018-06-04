# -*- coding: utf-8 -*-
"""
Tests for Timestamp parsing, aimed at pandas/_libs/tslibs/parsing.pyx
"""
from datetime import datetime
import numpy as np
import pytest
from dateutil.parser import parse

import pandas.util._test_decorators as td
from pandas import compat
from pandas.util import testing as tm
from pandas._libs.tslibs import parsing
from pandas._libs.tslibs.parsing import parse_time_string


class TestParseQuarters(object):

    def test_parse_time_string(self):
        (date, parsed, reso) = parse_time_string('4Q1984')
        (date_lower, parsed_lower, reso_lower) = parse_time_string('4q1984')
        assert date == date_lower
        assert parsed == parsed_lower
        assert reso == reso_lower

    def test_parse_time_quarter_w_dash(self):
        # https://github.com/pandas-dev/pandas/issue/9688
        pairs = [('1988-Q2', '1988Q2'), ('2Q-1988', '2Q1988')]

        for dashed, normal in pairs:
            (date_dash, parsed_dash, reso_dash) = parse_time_string(dashed)
            (date, parsed, reso) = parse_time_string(normal)

            assert date_dash == date
            assert parsed_dash == parsed
            assert reso_dash == reso

        pytest.raises(parsing.DateParseError, parse_time_string, "-2Q1992")
        pytest.raises(parsing.DateParseError, parse_time_string, "2-Q1992")
        pytest.raises(parsing.DateParseError, parse_time_string, "4-4Q1992")


class TestDatetimeParsingWrappers(object):
    def test_does_not_convert_mixed_integer(self):
        bad_date_strings = ('-50000', '999', '123.1234', 'm', 'T')

        for bad_date_string in bad_date_strings:
            assert not parsing._does_string_look_like_datetime(bad_date_string)

        good_date_strings = ('2012-01-01',
                             '01/01/2012',
                             'Mon Sep 16, 2013',
                             '01012012',
                             '0101',
                             '1-1')

        for good_date_string in good_date_strings:
            assert parsing._does_string_look_like_datetime(good_date_string)

    def test_parsers_quarterly_with_freq(self):
        msg = ('Incorrect quarterly string is given, quarter '
               'must be between 1 and 4: 2013Q5')
        with tm.assert_raises_regex(parsing.DateParseError, msg):
            parsing.parse_time_string('2013Q5')

        # GH 5418
        msg = ('Unable to retrieve month information from given freq: '
               'INVLD-L-DEC-SAT')
        with tm.assert_raises_regex(parsing.DateParseError, msg):
            parsing.parse_time_string('2013Q1', freq='INVLD-L-DEC-SAT')

        cases = {('2013Q2', None): datetime(2013, 4, 1),
                 ('2013Q2', 'A-APR'): datetime(2012, 8, 1),
                 ('2013-Q2', 'A-DEC'): datetime(2013, 4, 1)}

        for (date_str, freq), exp in compat.iteritems(cases):
            result, _, _ = parsing.parse_time_string(date_str, freq=freq)
            assert result == exp

    def test_parsers_quarter_invalid(self):

        cases = ['2Q 2005', '2Q-200A', '2Q-200', '22Q2005', '6Q-20', '2Q200.']
        for case in cases:
            pytest.raises(ValueError, parsing.parse_time_string, case)

    def test_parsers_monthfreq(self):
        cases = {'201101': datetime(2011, 1, 1, 0, 0),
                 '200005': datetime(2000, 5, 1, 0, 0)}

        for date_str, expected in compat.iteritems(cases):
            result1, _, _ = parsing.parse_time_string(date_str, freq='M')
            assert result1 == expected


class TestGuessDatetimeFormat(object):

    @td.skip_if_not_us_locale
    @pytest.mark.parametrize(
        "string, format",
        [
            ('20111230', '%Y%m%d'),
            ('2011-12-30', '%Y-%m-%d'),
            ('30-12-2011', '%d-%m-%Y'),
            ('2011-12-30 00:00:00', '%Y-%m-%d %H:%M:%S'),
            ('2011-12-30T00:00:00', '%Y-%m-%dT%H:%M:%S'),
            ('2011-12-30 00:00:00.000000',
             '%Y-%m-%d %H:%M:%S.%f')])
    def test_guess_datetime_format_with_parseable_formats(
            self, string, format):
        result = parsing._guess_datetime_format(string)
        assert result == format

    @pytest.mark.parametrize(
        "dayfirst, expected",
        [
            (True, "%d/%m/%Y"),
            (False, "%m/%d/%Y")])
    def test_guess_datetime_format_with_dayfirst(self, dayfirst, expected):
        ambiguous_string = '01/01/2011'
        result = parsing._guess_datetime_format(
            ambiguous_string, dayfirst=dayfirst)
        assert result == expected

    @td.skip_if_has_locale
    @pytest.mark.parametrize(
        "string, format",
        [
            ('30/Dec/2011', '%d/%b/%Y'),
            ('30/December/2011', '%d/%B/%Y'),
            ('30/Dec/2011 00:00:00', '%d/%b/%Y %H:%M:%S')])
    def test_guess_datetime_format_with_locale_specific_formats(
            self, string, format):
        result = parsing._guess_datetime_format(string)
        assert result == format

    def test_guess_datetime_format_invalid_inputs(self):
        # A datetime string must include a year, month and a day for it
        # to be guessable, in addition to being a string that looks like
        # a datetime
        invalid_dts = [
            '2013',
            '01/2013',
            '12:00:00',
            '1/1/1/1',
            'this_is_not_a_datetime',
            '51a',
            9,
            datetime(2011, 1, 1),
        ]

        for invalid_dt in invalid_dts:
            assert parsing._guess_datetime_format(invalid_dt) is None

    @pytest.mark.parametrize(
        "string, format",
        [
            ('2011-1-1', '%Y-%m-%d'),
            ('30-1-2011', '%d-%m-%Y'),
            ('1/1/2011', '%m/%d/%Y'),
            ('2011-1-1 00:00:00', '%Y-%m-%d %H:%M:%S'),
            ('2011-1-1 0:0:0', '%Y-%m-%d %H:%M:%S'),
            ('2011-1-3T00:00:0', '%Y-%m-%dT%H:%M:%S')])
    def test_guess_datetime_format_nopadding(self, string, format):
        # GH 11142
        result = parsing._guess_datetime_format(string)
        assert result == format


class TestArrayToDatetime(object):
    def test_try_parse_dates(self):
        arr = np.array(['5/1/2000', '6/1/2000', '7/1/2000'], dtype=object)

        result = parsing.try_parse_dates(arr, dayfirst=True)
        expected = np.array([parse(d, dayfirst=True) for d in arr])
        tm.assert_numpy_array_equal(result, expected)
