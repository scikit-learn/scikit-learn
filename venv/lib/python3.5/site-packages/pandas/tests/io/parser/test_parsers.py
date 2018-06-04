# -*- coding: utf-8 -*-

import os
import pandas.util.testing as tm

from pandas import read_csv, read_table, DataFrame
import pandas.core.common as com
from pandas._libs.tslib import Timestamp
from pandas.compat import StringIO

from .common import ParserTests
from .header import HeaderTests
from .comment import CommentTests
from .dialect import DialectTests
from .quoting import QuotingTests
from .usecols import UsecolsTests
from .skiprows import SkipRowsTests
from .index_col import IndexColTests
from .na_values import NAvaluesTests
from .converters import ConverterTests
from .c_parser_only import CParserTests
from .parse_dates import ParseDatesTests
from .compression import CompressionTests
from .mangle_dupes import DupeColumnTests
from .multithread import MultithreadTests
from .python_parser_only import PythonParserTests
from .dtypes import DtypeTests


class BaseParser(CommentTests, CompressionTests,
                 ConverterTests, DialectTests,
                 DtypeTests, DupeColumnTests,
                 HeaderTests, IndexColTests,
                 MultithreadTests, NAvaluesTests,
                 ParseDatesTests, ParserTests,
                 SkipRowsTests, UsecolsTests,
                 QuotingTests):

    def read_csv(self, *args, **kwargs):
        raise NotImplementedError

    def read_table(self, *args, **kwargs):
        raise NotImplementedError

    def float_precision_choices(self):
        raise com.AbstractMethodError(self)

    def setup_method(self, method):
        self.dirpath = tm.get_data_path()
        self.csv1 = os.path.join(self.dirpath, 'test1.csv')
        self.csv2 = os.path.join(self.dirpath, 'test2.csv')
        self.xls1 = os.path.join(self.dirpath, 'test.xls')
        self.csv_shiftjs = os.path.join(self.dirpath, 'sauron.SHIFT_JIS.csv')


class TestCParserHighMemory(BaseParser, CParserTests):
    engine = 'c'
    low_memory = False
    float_precision_choices = [None, 'high', 'round_trip']

    def read_csv(self, *args, **kwds):
        kwds = kwds.copy()
        kwds['engine'] = self.engine
        kwds['low_memory'] = self.low_memory
        return read_csv(*args, **kwds)

    def read_table(self, *args, **kwds):
        kwds = kwds.copy()
        kwds['engine'] = self.engine
        kwds['low_memory'] = self.low_memory
        return read_table(*args, **kwds)


class TestCParserLowMemory(BaseParser, CParserTests):
    engine = 'c'
    low_memory = True
    float_precision_choices = [None, 'high', 'round_trip']

    def read_csv(self, *args, **kwds):
        kwds = kwds.copy()
        kwds['engine'] = self.engine
        kwds['low_memory'] = self.low_memory
        return read_csv(*args, **kwds)

    def read_table(self, *args, **kwds):
        kwds = kwds.copy()
        kwds['engine'] = self.engine
        kwds['low_memory'] = True
        return read_table(*args, **kwds)


class TestPythonParser(BaseParser, PythonParserTests):
    engine = 'python'
    float_precision_choices = [None]

    def read_csv(self, *args, **kwds):
        kwds = kwds.copy()
        kwds['engine'] = self.engine
        return read_csv(*args, **kwds)

    def read_table(self, *args, **kwds):
        kwds = kwds.copy()
        kwds['engine'] = self.engine
        return read_table(*args, **kwds)


class TestUnsortedUsecols(object):
    def test_override__set_noconvert_columns(self):
        # GH 17351 - usecols needs to be sorted in _setnoconvert_columns
        # based on the test_usecols_with_parse_dates test from usecols.py
        from pandas.io.parsers import CParserWrapper, TextFileReader

        s = """a,b,c,d,e
        0,1,20140101,0900,4
        0,1,20140102,1000,4"""

        parse_dates = [[1, 2]]
        cols = {
            'a': [0, 0],
            'c_d': [
                Timestamp('2014-01-01 09:00:00'),
                Timestamp('2014-01-02 10:00:00')
            ]
        }
        expected = DataFrame(cols, columns=['c_d', 'a'])

        class MyTextFileReader(TextFileReader):
            def __init__(self):
                self._currow = 0
                self.squeeze = False

        class MyCParserWrapper(CParserWrapper):
            def _set_noconvert_columns(self):
                if self.usecols_dtype == 'integer':
                    # self.usecols is a set, which is documented as unordered
                    # but in practice, a CPython set of integers is sorted.
                    # In other implementations this assumption does not hold.
                    # The following code simulates a different order, which
                    # before GH 17351 would cause the wrong columns to be
                    # converted via the parse_dates parameter
                    self.usecols = list(self.usecols)
                    self.usecols.reverse()
                return CParserWrapper._set_noconvert_columns(self)

        parser = MyTextFileReader()
        parser.options = {'usecols': [0, 2, 3],
                          'parse_dates': parse_dates,
                          'delimiter': ','}
        parser._engine = MyCParserWrapper(StringIO(s), **parser.options)
        df = parser.read()

        tm.assert_frame_equal(df, expected)
