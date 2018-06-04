from __future__ import division, print_function, absolute_import

import datetime
import os
import sys
from os.path import join as pjoin

if sys.version_info[0] >= 3:
    from io import StringIO
else:
    from cStringIO import StringIO

import numpy as np

from numpy.testing import (assert_array_almost_equal,
                           assert_array_equal, assert_equal, assert_)
import pytest
from pytest import raises as assert_raises

from scipy.io.arff.arffread import loadarff
from scipy.io.arff.arffread import read_header, parse_type, ParseArffError


data_path = pjoin(os.path.dirname(__file__), 'data')

test1 = pjoin(data_path, 'test1.arff')
test2 = pjoin(data_path, 'test2.arff')
test3 = pjoin(data_path, 'test3.arff')

test4 = pjoin(data_path, 'test4.arff')
test5 = pjoin(data_path, 'test5.arff')
test6 = pjoin(data_path, 'test6.arff')
test7 = pjoin(data_path, 'test7.arff')
test8 = pjoin(data_path, 'test8.arff')
expect4_data = [(0.1, 0.2, 0.3, 0.4, 'class1'),
                (-0.1, -0.2, -0.3, -0.4, 'class2'),
                (1, 2, 3, 4, 'class3')]
expected_types = ['numeric', 'numeric', 'numeric', 'numeric', 'nominal']

missing = pjoin(data_path, 'missing.arff')
expect_missing_raw = np.array([[1, 5], [2, 4], [np.nan, np.nan]])
expect_missing = np.empty(3, [('yop', float), ('yap', float)])
expect_missing['yop'] = expect_missing_raw[:, 0]
expect_missing['yap'] = expect_missing_raw[:, 1]


class TestData(object):
    def test1(self):
        # Parsing trivial file with nothing.
        self._test(test4)

    def test2(self):
        # Parsing trivial file with some comments in the data section.
        self._test(test5)

    def test3(self):
        # Parsing trivial file with nominal attribute of 1 character.
        self._test(test6)

    def _test(self, test_file):
        data, meta = loadarff(test_file)
        for i in range(len(data)):
            for j in range(4):
                assert_array_almost_equal(expect4_data[i][j], data[i][j])
        assert_equal(meta.types(), expected_types)

    def test_filelike(self):
        # Test reading from file-like object (StringIO)
        f1 = open(test1)
        data1, meta1 = loadarff(f1)
        f1.close()
        f2 = open(test1)
        data2, meta2 = loadarff(StringIO(f2.read()))
        f2.close()
        assert_(data1 == data2)
        assert_(repr(meta1) == repr(meta2))

    @pytest.mark.skipif(sys.version_info < (3, 6),
                        reason='Passing path-like objects to IO functions requires Python >= 3.6')
    def test_path(self):
        # Test reading from `pathlib.Path` object
        from pathlib import Path

        with open(test1) as f1:
            data1, meta1 = loadarff(f1)

        data2, meta2 = loadarff(Path(test1))

        assert_(data1 == data2)
        assert_(repr(meta1) == repr(meta2))

class TestMissingData(object):
    def test_missing(self):
        data, meta = loadarff(missing)
        for i in ['yop', 'yap']:
            assert_array_almost_equal(data[i], expect_missing[i])


class TestNoData(object):
    def test_nodata(self):
        # The file nodata.arff has no data in the @DATA section.
        # Reading it should result in an array with length 0.
        nodata_filename = os.path.join(data_path, 'nodata.arff')
        data, meta = loadarff(nodata_filename)
        expected_dtype = np.dtype([('sepallength', '<f8'),
                                   ('sepalwidth', '<f8'),
                                   ('petallength', '<f8'),
                                   ('petalwidth', '<f8'),
                                   ('class', 'S15')])
        assert_equal(data.dtype, expected_dtype)
        assert_equal(data.size, 0)


class TestHeader(object):
    def test_type_parsing(self):
        # Test parsing type of attribute from their value.
        ofile = open(test2)
        rel, attrs = read_header(ofile)
        ofile.close()

        expected = ['numeric', 'numeric', 'numeric', 'numeric', 'numeric',
                    'numeric', 'string', 'string', 'nominal', 'nominal']

        for i in range(len(attrs)):
            assert_(parse_type(attrs[i][1]) == expected[i])

    def test_badtype_parsing(self):
        # Test parsing wrong type of attribute from their value.
        ofile = open(test3)
        rel, attrs = read_header(ofile)
        ofile.close()

        for name, value in attrs:
            assert_raises(ParseArffError, parse_type, value)

    def test_fullheader1(self):
        # Parsing trivial header with nothing.
        ofile = open(test1)
        rel, attrs = read_header(ofile)
        ofile.close()

        # Test relation
        assert_(rel == 'test1')

        # Test numerical attributes
        assert_(len(attrs) == 5)
        for i in range(4):
            assert_(attrs[i][0] == 'attr%d' % i)
            assert_(attrs[i][1] == 'REAL')

        # Test nominal attribute
        assert_(attrs[4][0] == 'class')
        assert_(attrs[4][1] == '{class0, class1, class2, class3}')

    def test_dateheader(self):
        ofile = open(test7)
        rel, attrs = read_header(ofile)
        ofile.close()

        assert_(rel == 'test7')

        assert_(len(attrs) == 5)

        assert_(attrs[0][0] == 'attr_year')
        assert_(attrs[0][1] == 'DATE yyyy')

        assert_(attrs[1][0] == 'attr_month')
        assert_(attrs[1][1] == 'DATE yyyy-MM')

        assert_(attrs[2][0] == 'attr_date')
        assert_(attrs[2][1] == 'DATE yyyy-MM-dd')

        assert_(attrs[3][0] == 'attr_datetime_local')
        assert_(attrs[3][1] == 'DATE "yyyy-MM-dd HH:mm"')

        assert_(attrs[4][0] == 'attr_datetime_missing')
        assert_(attrs[4][1] == 'DATE "yyyy-MM-dd HH:mm"')

    def test_dateheader_unsupported(self):
        ofile = open(test8)
        rel, attrs = read_header(ofile)
        ofile.close()

        assert_(rel == 'test8')

        assert_(len(attrs) == 2)
        assert_(attrs[0][0] == 'attr_datetime_utc')
        assert_(attrs[0][1] == 'DATE "yyyy-MM-dd HH:mm Z"')

        assert_(attrs[1][0] == 'attr_datetime_full')
        assert_(attrs[1][1] == 'DATE "yy-MM-dd HH:mm:ss z"')


class TestDateAttribute(object):
    def setup_method(self):
        self.data, self.meta = loadarff(test7)

    def test_year_attribute(self):
        expected = np.array([
            '1999',
            '2004',
            '1817',
            '2100',
            '2013',
            '1631'
        ], dtype='datetime64[Y]')

        assert_array_equal(self.data["attr_year"], expected)

    def test_month_attribute(self):
        expected = np.array([
            '1999-01',
            '2004-12',
            '1817-04',
            '2100-09',
            '2013-11',
            '1631-10'
        ], dtype='datetime64[M]')

        assert_array_equal(self.data["attr_month"], expected)

    def test_date_attribute(self):
        expected = np.array([
            '1999-01-31',
            '2004-12-01',
            '1817-04-28',
            '2100-09-10',
            '2013-11-30',
            '1631-10-15'
        ], dtype='datetime64[D]')

        assert_array_equal(self.data["attr_date"], expected)

    def test_datetime_local_attribute(self):
        expected = np.array([
            datetime.datetime(year=1999, month=1, day=31, hour=0, minute=1),
            datetime.datetime(year=2004, month=12, day=1, hour=23, minute=59),
            datetime.datetime(year=1817, month=4, day=28, hour=13, minute=0),
            datetime.datetime(year=2100, month=9, day=10, hour=12, minute=0),
            datetime.datetime(year=2013, month=11, day=30, hour=4, minute=55),
            datetime.datetime(year=1631, month=10, day=15, hour=20, minute=4)
        ], dtype='datetime64[m]')

        assert_array_equal(self.data["attr_datetime_local"], expected)

    def test_datetime_missing(self):
        expected = np.array([
            'nat',
            '2004-12-01T23:59',
            'nat',
            'nat',
            '2013-11-30T04:55',
            '1631-10-15T20:04'
        ], dtype='datetime64[m]')

        assert_array_equal(self.data["attr_datetime_missing"], expected)

    def test_datetime_timezone(self):
        assert_raises(ValueError, loadarff, test8)
