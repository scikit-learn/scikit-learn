# pylint: disable=E1101
import os
import warnings
from datetime import datetime, date, time, timedelta
from distutils.version import LooseVersion
from functools import partial
from warnings import catch_warnings
from collections import OrderedDict

import numpy as np
import pytest
from numpy import nan

import pandas as pd
import pandas.util.testing as tm
import pandas.util._test_decorators as td
from pandas import DataFrame, Index, MultiIndex
from pandas.compat import u, range, map, BytesIO, iteritems, PY36
from pandas.core.config import set_option, get_option
from pandas.io.common import URLError
from pandas.io.excel import (
    ExcelFile, ExcelWriter, read_excel, _XlwtWriter, _OpenpyxlWriter,
    register_writer, _XlsxWriter
)
from pandas.io.formats.excel import ExcelFormatter
from pandas.io.parsers import read_csv
from pandas.util.testing import ensure_clean, makeCustomDataframe as mkdf


_seriesd = tm.getSeriesData()
_tsd = tm.getTimeSeriesData()
_frame = DataFrame(_seriesd)[:10]
_frame2 = DataFrame(_seriesd, columns=['D', 'C', 'B', 'A'])[:10]
_tsframe = tm.makeTimeDataFrame()[:5]
_mixed_frame = _frame.copy()
_mixed_frame['foo'] = 'bar'


@td.skip_if_no('xlrd', '0.9')
class SharedItems(object):

    def setup_method(self, method):
        self.dirpath = tm.get_data_path()
        self.frame = _frame.copy()
        self.frame2 = _frame2.copy()
        self.tsframe = _tsframe.copy()
        self.mixed_frame = _mixed_frame.copy()

    def get_csv_refdf(self, basename):
        """
        Obtain the reference data from read_csv with the Python engine.
        Test data path is defined by pandas.util.testing.get_data_path()

        Parameters
        ----------

        basename : str
            File base name, excluding file extension.

        Returns
        -------

        dfref : DataFrame
        """
        pref = os.path.join(self.dirpath, basename + '.csv')
        dfref = read_csv(pref, index_col=0, parse_dates=True, engine='python')
        return dfref

    def get_excelfile(self, basename, ext):
        """
        Return test data ExcelFile instance. Test data path is defined by
        pandas.util.testing.get_data_path()

        Parameters
        ----------

        basename : str
            File base name, excluding file extension.

        Returns
        -------

        excel : io.excel.ExcelFile
        """
        return ExcelFile(os.path.join(self.dirpath, basename + ext))

    def get_exceldf(self, basename, ext, *args, **kwds):
        """
        Return test data DataFrame. Test data path is defined by
        pandas.util.testing.get_data_path()

        Parameters
        ----------

        basename : str
            File base name, excluding file extension.

        Returns
        -------

        df : DataFrame
        """
        pth = os.path.join(self.dirpath, basename + ext)
        return read_excel(pth, *args, **kwds)


class ReadingTestsBase(SharedItems):
    # This is based on ExcelWriterBase

    def test_usecols_int(self, ext):

        dfref = self.get_csv_refdf('test1')
        dfref = dfref.reindex(columns=['A', 'B', 'C'])
        df1 = self.get_exceldf('test1', ext, 'Sheet1', index_col=0, usecols=3)
        df2 = self.get_exceldf('test1', ext, 'Sheet2', skiprows=[1],
                               index_col=0, usecols=3)

        with tm.assert_produces_warning(FutureWarning):
            df3 = self.get_exceldf('test1', ext, 'Sheet2', skiprows=[1],
                                   index_col=0, parse_cols=3)

        # TODO add index to xls file)
        tm.assert_frame_equal(df1, dfref, check_names=False)
        tm.assert_frame_equal(df2, dfref, check_names=False)
        tm.assert_frame_equal(df3, dfref, check_names=False)

    def test_usecols_list(self, ext):

        dfref = self.get_csv_refdf('test1')
        dfref = dfref.reindex(columns=['B', 'C'])
        df1 = self.get_exceldf('test1', ext, 'Sheet1', index_col=0,
                               usecols=[0, 2, 3])
        df2 = self.get_exceldf('test1', ext, 'Sheet2', skiprows=[1],
                               index_col=0, usecols=[0, 2, 3])

        with tm.assert_produces_warning(FutureWarning):
            df3 = self.get_exceldf('test1', ext, 'Sheet2', skiprows=[1],
                                   index_col=0, parse_cols=[0, 2, 3])

        # TODO add index to xls file)
        tm.assert_frame_equal(df1, dfref, check_names=False)
        tm.assert_frame_equal(df2, dfref, check_names=False)
        tm.assert_frame_equal(df3, dfref, check_names=False)

    def test_usecols_str(self, ext):

        dfref = self.get_csv_refdf('test1')

        df1 = dfref.reindex(columns=['A', 'B', 'C'])
        df2 = self.get_exceldf('test1', ext, 'Sheet1', index_col=0,
                               usecols='A:D')
        df3 = self.get_exceldf('test1', ext, 'Sheet2', skiprows=[1],
                               index_col=0, usecols='A:D')

        with tm.assert_produces_warning(FutureWarning):
            df4 = self.get_exceldf('test1', ext, 'Sheet2', skiprows=[1],
                                   index_col=0, parse_cols='A:D')

        # TODO add index to xls, read xls ignores index name ?
        tm.assert_frame_equal(df2, df1, check_names=False)
        tm.assert_frame_equal(df3, df1, check_names=False)
        tm.assert_frame_equal(df4, df1, check_names=False)

        df1 = dfref.reindex(columns=['B', 'C'])
        df2 = self.get_exceldf('test1', ext, 'Sheet1', index_col=0,
                               usecols='A,C,D')
        df3 = self.get_exceldf('test1', ext, 'Sheet2', skiprows=[1],
                               index_col=0, usecols='A,C,D')
        # TODO add index to xls file
        tm.assert_frame_equal(df2, df1, check_names=False)
        tm.assert_frame_equal(df3, df1, check_names=False)

        df1 = dfref.reindex(columns=['B', 'C'])
        df2 = self.get_exceldf('test1', ext, 'Sheet1', index_col=0,
                               usecols='A,C:D')
        df3 = self.get_exceldf('test1', ext, 'Sheet2', skiprows=[1],
                               index_col=0, usecols='A,C:D')
        tm.assert_frame_equal(df2, df1, check_names=False)
        tm.assert_frame_equal(df3, df1, check_names=False)

    def test_excel_stop_iterator(self, ext):

        parsed = self.get_exceldf('test2', ext, 'Sheet1')
        expected = DataFrame([['aaaa', 'bbbbb']], columns=['Test', 'Test1'])
        tm.assert_frame_equal(parsed, expected)

    def test_excel_cell_error_na(self, ext):

        parsed = self.get_exceldf('test3', ext, 'Sheet1')
        expected = DataFrame([[np.nan]], columns=['Test'])
        tm.assert_frame_equal(parsed, expected)

    def test_excel_passes_na(self, ext):

        excel = self.get_excelfile('test4', ext)

        parsed = read_excel(excel, 'Sheet1', keep_default_na=False,
                            na_values=['apple'])
        expected = DataFrame([['NA'], [1], ['NA'], [np.nan], ['rabbit']],
                             columns=['Test'])
        tm.assert_frame_equal(parsed, expected)

        parsed = read_excel(excel, 'Sheet1', keep_default_na=True,
                            na_values=['apple'])
        expected = DataFrame([[np.nan], [1], [np.nan], [np.nan], ['rabbit']],
                             columns=['Test'])
        tm.assert_frame_equal(parsed, expected)

        # 13967
        excel = self.get_excelfile('test5', ext)

        parsed = read_excel(excel, 'Sheet1', keep_default_na=False,
                            na_values=['apple'])
        expected = DataFrame([['1.#QNAN'], [1], ['nan'], [np.nan], ['rabbit']],
                             columns=['Test'])
        tm.assert_frame_equal(parsed, expected)

        parsed = read_excel(excel, 'Sheet1', keep_default_na=True,
                            na_values=['apple'])
        expected = DataFrame([[np.nan], [1], [np.nan], [np.nan], ['rabbit']],
                             columns=['Test'])
        tm.assert_frame_equal(parsed, expected)

    def test_excel_table_sheet_by_index(self, ext):

        excel = self.get_excelfile('test1', ext)
        dfref = self.get_csv_refdf('test1')

        df1 = read_excel(excel, 0, index_col=0)
        df2 = read_excel(excel, 1, skiprows=[1], index_col=0)
        tm.assert_frame_equal(df1, dfref, check_names=False)
        tm.assert_frame_equal(df2, dfref, check_names=False)

        df1 = excel.parse(0, index_col=0)
        df2 = excel.parse(1, skiprows=[1], index_col=0)
        tm.assert_frame_equal(df1, dfref, check_names=False)
        tm.assert_frame_equal(df2, dfref, check_names=False)

        df3 = read_excel(excel, 0, index_col=0, skipfooter=1)
        tm.assert_frame_equal(df3, df1.iloc[:-1])

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            df4 = read_excel(excel, 0, index_col=0, skip_footer=1)
            tm.assert_frame_equal(df3, df4)

        df3 = excel.parse(0, index_col=0, skipfooter=1)
        tm.assert_frame_equal(df3, df1.iloc[:-1])

        import xlrd
        with pytest.raises(xlrd.XLRDError):
            read_excel(excel, 'asdf')

    def test_excel_table(self, ext):

        dfref = self.get_csv_refdf('test1')

        df1 = self.get_exceldf('test1', ext, 'Sheet1', index_col=0)
        df2 = self.get_exceldf('test1', ext, 'Sheet2', skiprows=[1],
                               index_col=0)
        # TODO add index to file
        tm.assert_frame_equal(df1, dfref, check_names=False)
        tm.assert_frame_equal(df2, dfref, check_names=False)

        df3 = self.get_exceldf('test1', ext, 'Sheet1', index_col=0,
                               skipfooter=1)
        tm.assert_frame_equal(df3, df1.iloc[:-1])

    def test_reader_special_dtypes(self, ext):

        expected = DataFrame.from_dict(OrderedDict([
            ("IntCol", [1, 2, -3, 4, 0]),
            ("FloatCol", [1.25, 2.25, 1.83, 1.92, 0.0000000005]),
            ("BoolCol", [True, False, True, True, False]),
            ("StrCol", [1, 2, 3, 4, 5]),
            # GH5394 - this is why convert_float isn't vectorized
            ("Str2Col", ["a", 3, "c", "d", "e"]),
            ("DateCol", [datetime(2013, 10, 30), datetime(2013, 10, 31),
                         datetime(1905, 1, 1), datetime(2013, 12, 14),
                         datetime(2015, 3, 14)])
        ]))
        basename = 'test_types'

        # should read in correctly and infer types
        actual = self.get_exceldf(basename, ext, 'Sheet1')
        tm.assert_frame_equal(actual, expected)

        # if not coercing number, then int comes in as float
        float_expected = expected.copy()
        float_expected["IntCol"] = float_expected["IntCol"].astype(float)
        float_expected.loc[float_expected.index[1], "Str2Col"] = 3.0
        actual = self.get_exceldf(basename, ext, 'Sheet1', convert_float=False)
        tm.assert_frame_equal(actual, float_expected)

        # check setting Index (assuming xls and xlsx are the same here)
        for icol, name in enumerate(expected.columns):
            actual = self.get_exceldf(basename, ext, 'Sheet1', index_col=icol)
            exp = expected.set_index(name)
            tm.assert_frame_equal(actual, exp)

        # convert_float and converters should be different but both accepted
        expected["StrCol"] = expected["StrCol"].apply(str)
        actual = self.get_exceldf(
            basename, ext, 'Sheet1', converters={"StrCol": str})
        tm.assert_frame_equal(actual, expected)

        no_convert_float = float_expected.copy()
        no_convert_float["StrCol"] = no_convert_float["StrCol"].apply(str)
        actual = self.get_exceldf(basename, ext, 'Sheet1', convert_float=False,
                                  converters={"StrCol": str})
        tm.assert_frame_equal(actual, no_convert_float)

    # GH8212 - support for converters and missing values
    def test_reader_converters(self, ext):

        basename = 'test_converters'

        expected = DataFrame.from_dict(OrderedDict([
            ("IntCol", [1, 2, -3, -1000, 0]),
            ("FloatCol", [12.5, np.nan, 18.3, 19.2, 0.000000005]),
            ("BoolCol", ['Found', 'Found', 'Found', 'Not found', 'Found']),
            ("StrCol", ['1', np.nan, '3', '4', '5']),
        ]))

        converters = {'IntCol': lambda x: int(x) if x != '' else -1000,
                      'FloatCol': lambda x: 10 * x if x else np.nan,
                      2: lambda x: 'Found' if x != '' else 'Not found',
                      3: lambda x: str(x) if x else '',
                      }

        # should read in correctly and set types of single cells (not array
        # dtypes)
        actual = self.get_exceldf(basename, ext, 'Sheet1',
                                  converters=converters)
        tm.assert_frame_equal(actual, expected)

    def test_reader_dtype(self, ext):
        # GH 8212
        basename = 'testdtype'
        actual = self.get_exceldf(basename, ext)

        expected = DataFrame({
            'a': [1, 2, 3, 4],
            'b': [2.5, 3.5, 4.5, 5.5],
            'c': [1, 2, 3, 4],
            'd': [1.0, 2.0, np.nan, 4.0]}).reindex(
                columns=['a', 'b', 'c', 'd'])

        tm.assert_frame_equal(actual, expected)

        actual = self.get_exceldf(basename, ext,
                                  dtype={'a': 'float64',
                                         'b': 'float32',
                                         'c': str})

        expected['a'] = expected['a'].astype('float64')
        expected['b'] = expected['b'].astype('float32')
        expected['c'] = ['001', '002', '003', '004']
        tm.assert_frame_equal(actual, expected)

        with pytest.raises(ValueError):
            actual = self.get_exceldf(basename, ext, dtype={'d': 'int64'})

    def test_reading_all_sheets(self, ext):
        # Test reading all sheetnames by setting sheetname to None,
        # Ensure a dict is returned.
        # See PR #9450
        basename = 'test_multisheet'
        dfs = self.get_exceldf(basename, ext, sheet_name=None)
        # ensure this is not alphabetical to test order preservation
        expected_keys = ['Charlie', 'Alpha', 'Beta']
        tm.assert_contains_all(expected_keys, dfs.keys())
        # Issue 9930
        # Ensure sheet order is preserved
        assert expected_keys == list(dfs.keys())

    def test_reading_multiple_specific_sheets(self, ext):
        # Test reading specific sheetnames by specifying a mixed list
        # of integers and strings, and confirm that duplicated sheet
        # references (positions/names) are removed properly.
        # Ensure a dict is returned
        # See PR #9450
        basename = 'test_multisheet'
        # Explicitly request duplicates. Only the set should be returned.
        expected_keys = [2, 'Charlie', 'Charlie']
        dfs = self.get_exceldf(basename, ext, sheet_name=expected_keys)
        expected_keys = list(set(expected_keys))
        tm.assert_contains_all(expected_keys, dfs.keys())
        assert len(expected_keys) == len(dfs.keys())

    def test_reading_all_sheets_with_blank(self, ext):
        # Test reading all sheetnames by setting sheetname to None,
        # In the case where some sheets are blank.
        # Issue #11711
        basename = 'blank_with_header'
        dfs = self.get_exceldf(basename, ext, sheet_name=None)
        expected_keys = ['Sheet1', 'Sheet2', 'Sheet3']
        tm.assert_contains_all(expected_keys, dfs.keys())

    # GH6403
    def test_read_excel_blank(self, ext):
        actual = self.get_exceldf('blank', ext, 'Sheet1')
        tm.assert_frame_equal(actual, DataFrame())

    def test_read_excel_blank_with_header(self, ext):
        expected = DataFrame(columns=['col_1', 'col_2'])
        actual = self.get_exceldf('blank_with_header', ext, 'Sheet1')
        tm.assert_frame_equal(actual, expected)

    @td.skip_if_no('openpyxl')
    @td.skip_if_no('xlwt')
    # GH 12292 : error when read one empty column from excel file
    def test_read_one_empty_col_no_header(self, ext):
        df = pd.DataFrame(
            [["", 1, 100],
             ["", 2, 200],
             ["", 3, 300],
             ["", 4, 400]]
        )
        with ensure_clean(ext) as path:
            df.to_excel(path, 'no_header', index=False, header=False)
            actual_header_none = read_excel(
                path,
                'no_header',
                usecols=[0],
                header=None
            )

            actual_header_zero = read_excel(
                path,
                'no_header',
                usecols=[0],
                header=0
            )
        expected = DataFrame()
        tm.assert_frame_equal(actual_header_none, expected)
        tm.assert_frame_equal(actual_header_zero, expected)

    @td.skip_if_no('openpyxl')
    @td.skip_if_no('xlwt')
    def test_read_one_empty_col_with_header(self, ext):
        df = pd.DataFrame(
            [["", 1, 100],
             ["", 2, 200],
             ["", 3, 300],
             ["", 4, 400]]
        )
        with ensure_clean(ext) as path:
            df.to_excel(path, 'with_header', index=False, header=True)
            actual_header_none = read_excel(
                path,
                'with_header',
                usecols=[0],
                header=None
            )

            actual_header_zero = read_excel(
                path,
                'with_header',
                usecols=[0],
                header=0
            )
        expected_header_none = DataFrame(pd.Series([0], dtype='int64'))
        tm.assert_frame_equal(actual_header_none, expected_header_none)
        expected_header_zero = DataFrame(columns=[0])
        tm.assert_frame_equal(actual_header_zero, expected_header_zero)

    @td.skip_if_no('openpyxl')
    @td.skip_if_no('xlwt')
    def test_set_column_names_in_parameter(self, ext):
        # GH 12870 : pass down column names associated with
        # keyword argument names
        refdf = pd.DataFrame([[1, 'foo'], [2, 'bar'],
                              [3, 'baz']], columns=['a', 'b'])

        with ensure_clean(ext) as pth:
            with ExcelWriter(pth) as writer:
                refdf.to_excel(writer, 'Data_no_head',
                               header=False, index=False)
                refdf.to_excel(writer, 'Data_with_head', index=False)

            refdf.columns = ['A', 'B']

            with ExcelFile(pth) as reader:
                xlsdf_no_head = read_excel(reader, 'Data_no_head',
                                           header=None, names=['A', 'B'])
                xlsdf_with_head = read_excel(reader, 'Data_with_head',
                                             index_col=None, names=['A', 'B'])

            tm.assert_frame_equal(xlsdf_no_head, refdf)
            tm.assert_frame_equal(xlsdf_with_head, refdf)

    def test_date_conversion_overflow(self, ext):
        # GH 10001 : pandas.ExcelFile ignore parse_dates=False
        expected = pd.DataFrame([[pd.Timestamp('2016-03-12'), 'Marc Johnson'],
                                 [pd.Timestamp('2016-03-16'), 'Jack Black'],
                                 [1e+20, 'Timothy Brown']],
                                columns=['DateColWithBigInt', 'StringCol'])

        result = self.get_exceldf('testdateoverflow', ext)
        tm.assert_frame_equal(result, expected)

    def test_sheet_name_and_sheetname(self, ext):
        # GH10559: Minor improvement: Change "sheet_name" to "sheetname"
        # GH10969: DOC: Consistent var names (sheetname vs sheet_name)
        # GH12604: CLN GH10559 Rename sheetname variable to sheet_name
        # GH20920: ExcelFile.parse() and pd.read_xlsx() have different
        #          behavior for "sheetname" argument
        dfref = self.get_csv_refdf('test1')
        df1 = self.get_exceldf('test1', ext,
                               sheet_name='Sheet1')  # doc
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            df2 = self.get_exceldf('test1', ext,
                                   sheetname='Sheet1')  # bkwrd compat

        excel = self.get_excelfile('test1', ext)
        df1_parse = excel.parse(sheet_name='Sheet1')    # doc
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            df2_parse = excel.parse(sheetname='Sheet1')  # bkwrd compat

        tm.assert_frame_equal(df1, dfref, check_names=False)
        tm.assert_frame_equal(df2, dfref, check_names=False)
        tm.assert_frame_equal(df1_parse, dfref, check_names=False)
        tm.assert_frame_equal(df2_parse, dfref, check_names=False)

    def test_sheet_name_both_raises(self, ext):
        with tm.assert_raises_regex(TypeError, "Cannot specify both"):
            self.get_exceldf('test1', ext, sheetname='Sheet1',
                             sheet_name='Sheet1')

        excel = self.get_excelfile('test1', ext)
        with tm.assert_raises_regex(TypeError, "Cannot specify both"):
            excel.parse(sheetname='Sheet1',
                        sheet_name='Sheet1')


@pytest.mark.parametrize("ext", ['.xls', '.xlsx', '.xlsm'])
class TestXlrdReader(ReadingTestsBase):
    """
    This is the base class for the xlrd tests, and 3 different file formats
    are supported: xls, xlsx, xlsm
    """

    def test_excel_read_buffer(self, ext):

        pth = os.path.join(self.dirpath, 'test1' + ext)
        expected = read_excel(pth, 'Sheet1', index_col=0)
        with open(pth, 'rb') as f:
            actual = read_excel(f, 'Sheet1', index_col=0)
            tm.assert_frame_equal(expected, actual)

        with open(pth, 'rb') as f:
            xls = ExcelFile(f)
            actual = read_excel(xls, 'Sheet1', index_col=0)
            tm.assert_frame_equal(expected, actual)

    @td.skip_if_no('xlwt')
    def test_read_xlrd_Book(self, ext):
        import xlrd

        df = self.frame
        with ensure_clean('.xls') as pth:
            df.to_excel(pth, "SheetA")
            book = xlrd.open_workbook(pth)

            with ExcelFile(book, engine="xlrd") as xl:
                result = read_excel(xl, "SheetA")
                tm.assert_frame_equal(df, result)

            result = read_excel(book, sheet_name="SheetA", engine="xlrd")
            tm.assert_frame_equal(df, result)

    @tm.network
    def test_read_from_http_url(self, ext):
        url = ('https://raw.github.com/pandas-dev/pandas/master/'
               'pandas/tests/io/data/test1' + ext)
        url_table = read_excel(url)
        local_table = self.get_exceldf('test1', ext)
        tm.assert_frame_equal(url_table, local_table)

    @td.skip_if_no('s3fs')
    def test_read_from_s3_url(self, ext):
        boto3 = pytest.importorskip('boto3')
        moto = pytest.importorskip('moto')

        with moto.mock_s3():
            conn = boto3.resource("s3", region_name="us-east-1")
            conn.create_bucket(Bucket="pandas-test")
            file_name = os.path.join(self.dirpath, 'test1' + ext)
            with open(file_name, 'rb') as f:
                conn.Bucket("pandas-test").put_object(Key="test1" + ext,
                                                      Body=f)

            url = ('s3://pandas-test/test1' + ext)
            url_table = read_excel(url)
            local_table = self.get_exceldf('test1', ext)
            tm.assert_frame_equal(url_table, local_table)

    @pytest.mark.slow
    def test_read_from_file_url(self, ext):

        # FILE
        localtable = os.path.join(self.dirpath, 'test1' + ext)
        local_table = read_excel(localtable)

        try:
            url_table = read_excel('file://localhost/' + localtable)
        except URLError:
            # fails on some systems
            import platform
            pytest.skip("failing on %s" %
                        ' '.join(platform.uname()).strip())

        tm.assert_frame_equal(url_table, local_table)

    @td.skip_if_no('pathlib')
    def test_read_from_pathlib_path(self, ext):

        # GH12655
        from pathlib import Path

        str_path = os.path.join(self.dirpath, 'test1' + ext)
        expected = read_excel(str_path, 'Sheet1', index_col=0)

        path_obj = Path(self.dirpath, 'test1' + ext)
        actual = read_excel(path_obj, 'Sheet1', index_col=0)

        tm.assert_frame_equal(expected, actual)

    @td.skip_if_no('py.path')
    def test_read_from_py_localpath(self, ext):

        # GH12655
        from py.path import local as LocalPath

        str_path = os.path.join(self.dirpath, 'test1' + ext)
        expected = read_excel(str_path, 'Sheet1', index_col=0)

        abs_dir = os.path.abspath(self.dirpath)
        path_obj = LocalPath(abs_dir).join('test1' + ext)
        actual = read_excel(path_obj, 'Sheet1', index_col=0)

        tm.assert_frame_equal(expected, actual)

    def test_reader_closes_file(self, ext):

        pth = os.path.join(self.dirpath, 'test1' + ext)
        f = open(pth, 'rb')
        with ExcelFile(f) as xlsx:
            # parses okay
            read_excel(xlsx, 'Sheet1', index_col=0)

        assert f.closed

    @td.skip_if_no('openpyxl')
    @td.skip_if_no('xlwt')
    def test_creating_and_reading_multiple_sheets(self, ext):
        # Test reading multiple sheets, from a runtime created excel file
        # with multiple sheets.
        # See PR #9450
        def tdf(sheetname):
            d, i = [11, 22, 33], [1, 2, 3]
            return DataFrame(d, i, columns=[sheetname])

        sheets = ['AAA', 'BBB', 'CCC']

        dfs = [tdf(s) for s in sheets]
        dfs = dict(zip(sheets, dfs))

        with ensure_clean(ext) as pth:
            with ExcelWriter(pth) as ew:
                for sheetname, df in iteritems(dfs):
                    df.to_excel(ew, sheetname)
            dfs_returned = read_excel(pth, sheet_name=sheets)
            for s in sheets:
                tm.assert_frame_equal(dfs[s], dfs_returned[s])

    def test_reader_seconds(self, ext):
        import xlrd

        # Test reading times with and without milliseconds. GH5945.
        if LooseVersion(xlrd.__VERSION__) >= LooseVersion("0.9.3"):
            # Xlrd >= 0.9.3 can handle Excel milliseconds.
            expected = DataFrame.from_dict({"Time": [time(1, 2, 3),
                                            time(2, 45, 56, 100000),
                                            time(4, 29, 49, 200000),
                                            time(6, 13, 42, 300000),
                                            time(7, 57, 35, 400000),
                                            time(9, 41, 28, 500000),
                                            time(11, 25, 21, 600000),
                                            time(13, 9, 14, 700000),
                                            time(14, 53, 7, 800000),
                                            time(16, 37, 0, 900000),
                                            time(18, 20, 54)]})
        else:
            # Xlrd < 0.9.3 rounds Excel milliseconds.
            expected = DataFrame.from_dict({"Time": [time(1, 2, 3),
                                            time(2, 45, 56),
                                            time(4, 29, 49),
                                            time(6, 13, 42),
                                            time(7, 57, 35),
                                            time(9, 41, 29),
                                            time(11, 25, 22),
                                            time(13, 9, 15),
                                            time(14, 53, 8),
                                            time(16, 37, 1),
                                            time(18, 20, 54)]})

        actual = self.get_exceldf('times_1900', ext, 'Sheet1')
        tm.assert_frame_equal(actual, expected)

        actual = self.get_exceldf('times_1904', ext, 'Sheet1')
        tm.assert_frame_equal(actual, expected)

    def test_read_excel_multiindex(self, ext):
        # GH 4679
        mi = MultiIndex.from_product([['foo', 'bar'], ['a', 'b']])
        mi_file = os.path.join(self.dirpath, 'testmultiindex' + ext)

        expected = DataFrame([[1, 2.5, pd.Timestamp('2015-01-01'), True],
                              [2, 3.5, pd.Timestamp('2015-01-02'), False],
                              [3, 4.5, pd.Timestamp('2015-01-03'), False],
                              [4, 5.5, pd.Timestamp('2015-01-04'), True]],
                             columns=mi)

        actual = read_excel(mi_file, 'mi_column', header=[0, 1])
        tm.assert_frame_equal(actual, expected)
        actual = read_excel(mi_file, 'mi_column', header=[0, 1], index_col=0)
        tm.assert_frame_equal(actual, expected)

        expected.columns = ['a', 'b', 'c', 'd']
        expected.index = mi
        actual = read_excel(mi_file, 'mi_index', index_col=[0, 1])
        tm.assert_frame_equal(actual, expected, check_names=False)

        expected.columns = mi
        actual = read_excel(mi_file, 'both', index_col=[0, 1], header=[0, 1])
        tm.assert_frame_equal(actual, expected, check_names=False)

        expected.index = mi.set_names(['ilvl1', 'ilvl2'])
        expected.columns = ['a', 'b', 'c', 'd']
        actual = read_excel(mi_file, 'mi_index_name', index_col=[0, 1])
        tm.assert_frame_equal(actual, expected)

        expected.index = list(range(4))
        expected.columns = mi.set_names(['c1', 'c2'])
        actual = read_excel(mi_file, 'mi_column_name',
                            header=[0, 1], index_col=0)
        tm.assert_frame_equal(actual, expected)

        # Issue #11317
        expected.columns = mi.set_levels(
            [1, 2], level=1).set_names(['c1', 'c2'])
        actual = read_excel(mi_file, 'name_with_int',
                            index_col=0, header=[0, 1])
        tm.assert_frame_equal(actual, expected)

        expected.columns = mi.set_names(['c1', 'c2'])
        expected.index = mi.set_names(['ilvl1', 'ilvl2'])
        actual = read_excel(mi_file, 'both_name',
                            index_col=[0, 1], header=[0, 1])
        tm.assert_frame_equal(actual, expected)

        actual = read_excel(mi_file, 'both_name',
                            index_col=[0, 1], header=[0, 1])
        tm.assert_frame_equal(actual, expected)

        actual = read_excel(mi_file, 'both_name_skiprows', index_col=[0, 1],
                            header=[0, 1], skiprows=2)
        tm.assert_frame_equal(actual, expected)

    @td.skip_if_no('xlsxwriter')
    def test_read_excel_multiindex_empty_level(self, ext):
        # GH 12453
        with ensure_clean('.xlsx') as path:
            df = DataFrame({
                ('One', 'x'): {0: 1},
                ('Two', 'X'): {0: 3},
                ('Two', 'Y'): {0: 7},
                ('Zero', ''): {0: 0}
            })

            expected = DataFrame({
                ('One', u'x'): {0: 1},
                ('Two', u'X'): {0: 3},
                ('Two', u'Y'): {0: 7},
                ('Zero', 'Unnamed: 3_level_1'): {0: 0}
            })

            df.to_excel(path)
            actual = pd.read_excel(path, header=[0, 1])
            tm.assert_frame_equal(actual, expected)

            df = pd.DataFrame({
                ('Beg', ''): {0: 0},
                ('Middle', 'x'): {0: 1},
                ('Tail', 'X'): {0: 3},
                ('Tail', 'Y'): {0: 7}
            })

            expected = pd.DataFrame({
                ('Beg', 'Unnamed: 0_level_1'): {0: 0},
                ('Middle', u'x'): {0: 1},
                ('Tail', u'X'): {0: 3},
                ('Tail', u'Y'): {0: 7}
            })

            df.to_excel(path)
            actual = pd.read_excel(path, header=[0, 1])
            tm.assert_frame_equal(actual, expected)

    @td.skip_if_no('xlsxwriter')
    def test_excel_multindex_roundtrip(self, ext):
        # GH 4679
        with ensure_clean('.xlsx') as pth:
            for c_idx_names in [True, False]:
                for r_idx_names in [True, False]:
                    for c_idx_levels in [1, 3]:
                        for r_idx_levels in [1, 3]:
                            # column index name can't be serialized unless
                            # MultiIndex
                            if (c_idx_levels == 1 and c_idx_names):
                                continue

                            # empty name case current read in as unnamed
                            # levels, not Nones
                            check_names = True
                            if not r_idx_names and r_idx_levels > 1:
                                check_names = False

                            df = mkdf(5, 5, c_idx_names,
                                      r_idx_names, c_idx_levels,
                                      r_idx_levels)
                            df.to_excel(pth)
                            act = pd.read_excel(
                                pth, index_col=list(range(r_idx_levels)),
                                header=list(range(c_idx_levels)))
                            tm.assert_frame_equal(
                                df, act, check_names=check_names)

                            df.iloc[0, :] = np.nan
                            df.to_excel(pth)
                            act = pd.read_excel(
                                pth, index_col=list(range(r_idx_levels)),
                                header=list(range(c_idx_levels)))
                            tm.assert_frame_equal(
                                df, act, check_names=check_names)

                            df.iloc[-1, :] = np.nan
                            df.to_excel(pth)
                            act = pd.read_excel(
                                pth, index_col=list(range(r_idx_levels)),
                                header=list(range(c_idx_levels)))
                            tm.assert_frame_equal(
                                df, act, check_names=check_names)

    def test_excel_old_index_format(self, ext):
        # see gh-4679
        filename = 'test_index_name_pre17' + ext
        in_file = os.path.join(self.dirpath, filename)

        # We detect headers to determine if index names exist, so
        # that "index" name in the "names" version of the data will
        # now be interpreted as rows that include null data.
        data = np.array([[None, None, None, None, None],
                         ['R0C0', 'R0C1', 'R0C2', 'R0C3', 'R0C4'],
                         ['R1C0', 'R1C1', 'R1C2', 'R1C3', 'R1C4'],
                         ['R2C0', 'R2C1', 'R2C2', 'R2C3', 'R2C4'],
                         ['R3C0', 'R3C1', 'R3C2', 'R3C3', 'R3C4'],
                         ['R4C0', 'R4C1', 'R4C2', 'R4C3', 'R4C4']])
        columns = ['C_l0_g0', 'C_l0_g1', 'C_l0_g2', 'C_l0_g3', 'C_l0_g4']
        mi = MultiIndex(levels=[['R0', 'R_l0_g0', 'R_l0_g1',
                                 'R_l0_g2', 'R_l0_g3', 'R_l0_g4'],
                                ['R1', 'R_l1_g0', 'R_l1_g1',
                                 'R_l1_g2', 'R_l1_g3', 'R_l1_g4']],
                        labels=[[0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5]],
                        names=[None, None])
        si = Index(['R0', 'R_l0_g0', 'R_l0_g1', 'R_l0_g2',
                    'R_l0_g3', 'R_l0_g4'], name=None)

        expected = pd.DataFrame(data, index=si, columns=columns)

        actual = pd.read_excel(in_file, 'single_names')
        tm.assert_frame_equal(actual, expected)

        expected.index = mi

        actual = pd.read_excel(in_file, 'multi_names')
        tm.assert_frame_equal(actual, expected)

        # The analogous versions of the "names" version data
        # where there are explicitly no names for the indices.
        data = np.array([['R0C0', 'R0C1', 'R0C2', 'R0C3', 'R0C4'],
                         ['R1C0', 'R1C1', 'R1C2', 'R1C3', 'R1C4'],
                         ['R2C0', 'R2C1', 'R2C2', 'R2C3', 'R2C4'],
                         ['R3C0', 'R3C1', 'R3C2', 'R3C3', 'R3C4'],
                         ['R4C0', 'R4C1', 'R4C2', 'R4C3', 'R4C4']])
        columns = ['C_l0_g0', 'C_l0_g1', 'C_l0_g2', 'C_l0_g3', 'C_l0_g4']
        mi = MultiIndex(levels=[['R_l0_g0', 'R_l0_g1', 'R_l0_g2',
                                 'R_l0_g3', 'R_l0_g4'],
                                ['R_l1_g0', 'R_l1_g1', 'R_l1_g2',
                                 'R_l1_g3', 'R_l1_g4']],
                        labels=[[0, 1, 2, 3, 4], [0, 1, 2, 3, 4]],
                        names=[None, None])
        si = Index(['R_l0_g0', 'R_l0_g1', 'R_l0_g2',
                    'R_l0_g3', 'R_l0_g4'], name=None)

        expected = pd.DataFrame(data, index=si, columns=columns)

        actual = pd.read_excel(in_file, 'single_no_names')
        tm.assert_frame_equal(actual, expected)

        expected.index = mi

        actual = pd.read_excel(in_file, 'multi_no_names', index_col=[0, 1])
        tm.assert_frame_equal(actual, expected, check_names=False)

    def test_read_excel_bool_header_arg(self, ext):
        # GH 6114
        for arg in [True, False]:
            with pytest.raises(TypeError):
                pd.read_excel(os.path.join(self.dirpath, 'test1' + ext),
                              header=arg)

    def test_read_excel_chunksize(self, ext):
        # GH 8011
        with pytest.raises(NotImplementedError):
            pd.read_excel(os.path.join(self.dirpath, 'test1' + ext),
                          chunksize=100)

    @td.skip_if_no('openpyxl')
    @td.skip_if_no('xlwt')
    def test_read_excel_parse_dates(self, ext):
        # GH 11544, 12051
        df = DataFrame(
            {'col': [1, 2, 3],
             'date_strings': pd.date_range('2012-01-01', periods=3)})
        df2 = df.copy()
        df2['date_strings'] = df2['date_strings'].dt.strftime('%m/%d/%Y')

        with ensure_clean(ext) as pth:
            df2.to_excel(pth)

            res = read_excel(pth)
            tm.assert_frame_equal(df2, res)

            # no index_col specified when parse_dates is True
            with tm.assert_produces_warning():
                res = read_excel(pth, parse_dates=True)
                tm.assert_frame_equal(df2, res)

            res = read_excel(pth, parse_dates=['date_strings'], index_col=0)
            tm.assert_frame_equal(df, res)

            dateparser = lambda x: pd.datetime.strptime(x, '%m/%d/%Y')
            res = read_excel(pth, parse_dates=['date_strings'],
                             date_parser=dateparser, index_col=0)
            tm.assert_frame_equal(df, res)

    def test_read_excel_skiprows_list(self, ext):
        # GH 4903
        actual = pd.read_excel(os.path.join(self.dirpath,
                                            'testskiprows' + ext),
                               'skiprows_list', skiprows=[0, 2])
        expected = DataFrame([[1, 2.5, pd.Timestamp('2015-01-01'), True],
                              [2, 3.5, pd.Timestamp('2015-01-02'), False],
                              [3, 4.5, pd.Timestamp('2015-01-03'), False],
                              [4, 5.5, pd.Timestamp('2015-01-04'), True]],
                             columns=['a', 'b', 'c', 'd'])
        tm.assert_frame_equal(actual, expected)

        actual = pd.read_excel(os.path.join(self.dirpath,
                                            'testskiprows' + ext),
                               'skiprows_list', skiprows=np.array([0, 2]))
        tm.assert_frame_equal(actual, expected)

    def test_read_excel_nrows(self, ext):
        # GH 16645
        num_rows_to_pull = 5
        actual = pd.read_excel(os.path.join(self.dirpath, 'test1' + ext),
                               nrows=num_rows_to_pull)
        expected = pd.read_excel(os.path.join(self.dirpath,
                                              'test1' + ext))
        expected = expected[:num_rows_to_pull]
        tm.assert_frame_equal(actual, expected)

    def test_read_excel_nrows_greater_than_nrows_in_file(self, ext):
        # GH 16645
        expected = pd.read_excel(os.path.join(self.dirpath,
                                              'test1' + ext))
        num_records_in_file = len(expected)
        num_rows_to_pull = num_records_in_file + 10
        actual = pd.read_excel(os.path.join(self.dirpath, 'test1' + ext),
                               nrows=num_rows_to_pull)
        tm.assert_frame_equal(actual, expected)

    def test_read_excel_nrows_non_integer_parameter(self, ext):
        # GH 16645
        msg = "'nrows' must be an integer >=0"
        with tm.assert_raises_regex(ValueError, msg):
            pd.read_excel(os.path.join(self.dirpath, 'test1' + ext),
                          nrows='5')

    def test_read_excel_squeeze(self, ext):
        # GH 12157
        f = os.path.join(self.dirpath, 'test_squeeze' + ext)

        actual = pd.read_excel(f, 'two_columns', index_col=0, squeeze=True)
        expected = pd.Series([2, 3, 4], [4, 5, 6], name='b')
        expected.index.name = 'a'
        tm.assert_series_equal(actual, expected)

        actual = pd.read_excel(f, 'two_columns', squeeze=True)
        expected = pd.DataFrame({'a': [4, 5, 6],
                                 'b': [2, 3, 4]})
        tm.assert_frame_equal(actual, expected)

        actual = pd.read_excel(f, 'one_column', squeeze=True)
        expected = pd.Series([1, 2, 3], name='a')
        tm.assert_series_equal(actual, expected)


class _WriterBase(SharedItems):

    @pytest.fixture(autouse=True)
    def set_engine_and_path(self, request, merge_cells, engine, ext):
        """Fixture to set engine and open file for use in each test case

        Rather than requiring `engine=...` to be provided explicitly as an
        argument in each test, this fixture sets a global option to dictate
        which engine should be used to write Excel files. After executing
        the test it rolls back said change to the global option.

        It also uses a context manager to open a temporary excel file for
        the function to write to, accessible via `self.path`

        Notes
        -----
        This fixture will run as part of each test method defined in the
        class and any subclasses, on account of the `autouse=True`
        argument
        """
        option_name = 'io.excel.{ext}.writer'.format(ext=ext.strip('.'))
        prev_engine = get_option(option_name)
        set_option(option_name, engine)
        with ensure_clean(ext) as path:
            self.path = path
            yield
        set_option(option_name, prev_engine)  # Roll back option change


@pytest.mark.parametrize("merge_cells", [True, False])
@pytest.mark.parametrize("engine,ext", [
    pytest.param('openpyxl', '.xlsx', marks=pytest.mark.skipif(
        not td.safe_import('openpyxl'), reason='No openpyxl')),
    pytest.param('openpyxl', '.xlsm', marks=pytest.mark.skipif(
        not td.safe_import('openpyxl'), reason='No openpyxl')),
    pytest.param('xlwt', '.xls', marks=pytest.mark.skipif(
        not td.safe_import('xlwt'), reason='No xlwt')),
    pytest.param('xlsxwriter', '.xlsx', marks=pytest.mark.skipif(
        not td.safe_import('xlsxwriter'), reason='No xlsxwriter'))
])
class TestExcelWriter(_WriterBase):
    # Base class for test cases to run with different Excel writers.

    def test_excel_sheet_by_name_raise(self, merge_cells, engine, ext):
        import xlrd

        gt = DataFrame(np.random.randn(10, 2))
        gt.to_excel(self.path)
        xl = ExcelFile(self.path)
        df = read_excel(xl, 0)
        tm.assert_frame_equal(gt, df)

        with pytest.raises(xlrd.XLRDError):
            read_excel(xl, '0')

    def test_excelwriter_contextmanager(self, merge_cells, engine, ext):
        with ExcelWriter(self.path) as writer:
            self.frame.to_excel(writer, 'Data1')
            self.frame2.to_excel(writer, 'Data2')

        with ExcelFile(self.path) as reader:
            found_df = read_excel(reader, 'Data1')
            found_df2 = read_excel(reader, 'Data2')
            tm.assert_frame_equal(found_df, self.frame)
            tm.assert_frame_equal(found_df2, self.frame2)

    def test_roundtrip(self, merge_cells, engine, ext):
        self.frame['A'][:5] = nan

        self.frame.to_excel(self.path, 'test1')
        self.frame.to_excel(self.path, 'test1', columns=['A', 'B'])
        self.frame.to_excel(self.path, 'test1', header=False)
        self.frame.to_excel(self.path, 'test1', index=False)

        # test roundtrip
        self.frame.to_excel(self.path, 'test1')
        recons = read_excel(self.path, 'test1', index_col=0)
        tm.assert_frame_equal(self.frame, recons)

        self.frame.to_excel(self.path, 'test1', index=False)
        recons = read_excel(self.path, 'test1', index_col=None)
        recons.index = self.frame.index
        tm.assert_frame_equal(self.frame, recons)

        self.frame.to_excel(self.path, 'test1', na_rep='NA')
        recons = read_excel(self.path, 'test1', index_col=0, na_values=['NA'])
        tm.assert_frame_equal(self.frame, recons)

        # GH 3611
        self.frame.to_excel(self.path, 'test1', na_rep='88')
        recons = read_excel(self.path, 'test1', index_col=0, na_values=['88'])
        tm.assert_frame_equal(self.frame, recons)

        self.frame.to_excel(self.path, 'test1', na_rep='88')
        recons = read_excel(self.path, 'test1', index_col=0,
                            na_values=[88, 88.0])
        tm.assert_frame_equal(self.frame, recons)

        # GH 6573
        self.frame.to_excel(self.path, 'Sheet1')
        recons = read_excel(self.path, index_col=0)
        tm.assert_frame_equal(self.frame, recons)

        self.frame.to_excel(self.path, '0')
        recons = read_excel(self.path, index_col=0)
        tm.assert_frame_equal(self.frame, recons)

        # GH 8825 Pandas Series should provide to_excel method
        s = self.frame["A"]
        s.to_excel(self.path)
        recons = read_excel(self.path, index_col=0)
        tm.assert_frame_equal(s.to_frame(), recons)

    def test_mixed(self, merge_cells, engine, ext):
        self.mixed_frame.to_excel(self.path, 'test1')
        reader = ExcelFile(self.path)
        recons = read_excel(reader, 'test1', index_col=0)
        tm.assert_frame_equal(self.mixed_frame, recons)

    def test_tsframe(self, merge_cells, engine, ext):
        df = tm.makeTimeDataFrame()[:5]

        df.to_excel(self.path, 'test1')
        reader = ExcelFile(self.path)
        recons = read_excel(reader, 'test1')
        tm.assert_frame_equal(df, recons)

    def test_basics_with_nan(self, merge_cells, engine, ext):
        self.frame['A'][:5] = nan
        self.frame.to_excel(self.path, 'test1')
        self.frame.to_excel(self.path, 'test1', columns=['A', 'B'])
        self.frame.to_excel(self.path, 'test1', header=False)
        self.frame.to_excel(self.path, 'test1', index=False)

    @pytest.mark.parametrize("np_type", [
        np.int8, np.int16, np.int32, np.int64])
    def test_int_types(self, merge_cells, engine, ext, np_type):
        # Test np.int values read come back as int (rather than float
        # which is Excel's format).
        frame = DataFrame(np.random.randint(-10, 10, size=(10, 2)),
                          dtype=np_type)
        frame.to_excel(self.path, 'test1')
        reader = ExcelFile(self.path)
        recons = read_excel(reader, 'test1')
        int_frame = frame.astype(np.int64)
        tm.assert_frame_equal(int_frame, recons)
        recons2 = read_excel(self.path, 'test1')
        tm.assert_frame_equal(int_frame, recons2)

        # test with convert_float=False comes back as float
        float_frame = frame.astype(float)
        recons = read_excel(self.path, 'test1', convert_float=False)
        tm.assert_frame_equal(recons, float_frame,
                              check_index_type=False,
                              check_column_type=False)

    @pytest.mark.parametrize("np_type", [
        np.float16, np.float32, np.float64])
    def test_float_types(self, merge_cells, engine, ext, np_type):
        # Test np.float values read come back as float.
        frame = DataFrame(np.random.random_sample(10), dtype=np_type)
        frame.to_excel(self.path, 'test1')
        reader = ExcelFile(self.path)
        recons = read_excel(reader, 'test1').astype(np_type)
        tm.assert_frame_equal(frame, recons, check_dtype=False)

    @pytest.mark.parametrize("np_type", [np.bool8, np.bool_])
    def test_bool_types(self, merge_cells, engine, ext, np_type):
        # Test np.bool values read come back as float.
        frame = (DataFrame([1, 0, True, False], dtype=np_type))
        frame.to_excel(self.path, 'test1')
        reader = ExcelFile(self.path)
        recons = read_excel(reader, 'test1').astype(np_type)
        tm.assert_frame_equal(frame, recons)

    def test_inf_roundtrip(self, merge_cells, engine, ext):
        frame = DataFrame([(1, np.inf), (2, 3), (5, -np.inf)])
        frame.to_excel(self.path, 'test1')
        reader = ExcelFile(self.path)
        recons = read_excel(reader, 'test1')
        tm.assert_frame_equal(frame, recons)

    def test_sheets(self, merge_cells, engine, ext):
        self.frame['A'][:5] = nan

        self.frame.to_excel(self.path, 'test1')
        self.frame.to_excel(self.path, 'test1', columns=['A', 'B'])
        self.frame.to_excel(self.path, 'test1', header=False)
        self.frame.to_excel(self.path, 'test1', index=False)

        # Test writing to separate sheets
        writer = ExcelWriter(self.path)
        self.frame.to_excel(writer, 'test1')
        self.tsframe.to_excel(writer, 'test2')
        writer.save()
        reader = ExcelFile(self.path)
        recons = read_excel(reader, 'test1', index_col=0)
        tm.assert_frame_equal(self.frame, recons)
        recons = read_excel(reader, 'test2', index_col=0)
        tm.assert_frame_equal(self.tsframe, recons)
        assert 2 == len(reader.sheet_names)
        assert 'test1' == reader.sheet_names[0]
        assert 'test2' == reader.sheet_names[1]

    def test_colaliases(self, merge_cells, engine, ext):
        self.frame['A'][:5] = nan

        self.frame.to_excel(self.path, 'test1')
        self.frame.to_excel(self.path, 'test1', columns=['A', 'B'])
        self.frame.to_excel(self.path, 'test1', header=False)
        self.frame.to_excel(self.path, 'test1', index=False)

        # column aliases
        col_aliases = Index(['AA', 'X', 'Y', 'Z'])
        self.frame2.to_excel(self.path, 'test1', header=col_aliases)
        reader = ExcelFile(self.path)
        rs = read_excel(reader, 'test1', index_col=0)
        xp = self.frame2.copy()
        xp.columns = col_aliases
        tm.assert_frame_equal(xp, rs)

    def test_roundtrip_indexlabels(self, merge_cells, engine, ext):
        self.frame['A'][:5] = nan

        self.frame.to_excel(self.path, 'test1')
        self.frame.to_excel(self.path, 'test1', columns=['A', 'B'])
        self.frame.to_excel(self.path, 'test1', header=False)
        self.frame.to_excel(self.path, 'test1', index=False)

        # test index_label
        frame = (DataFrame(np.random.randn(10, 2)) >= 0)
        frame.to_excel(self.path, 'test1',
                       index_label=['test'],
                       merge_cells=merge_cells)
        reader = ExcelFile(self.path)
        recons = read_excel(reader, 'test1',
                            index_col=0,
                            ).astype(np.int64)
        frame.index.names = ['test']
        assert frame.index.names == recons.index.names

        frame = (DataFrame(np.random.randn(10, 2)) >= 0)
        frame.to_excel(self.path,
                       'test1',
                       index_label=['test', 'dummy', 'dummy2'],
                       merge_cells=merge_cells)
        reader = ExcelFile(self.path)
        recons = read_excel(reader, 'test1',
                            index_col=0,
                            ).astype(np.int64)
        frame.index.names = ['test']
        assert frame.index.names == recons.index.names

        frame = (DataFrame(np.random.randn(10, 2)) >= 0)
        frame.to_excel(self.path,
                       'test1',
                       index_label='test',
                       merge_cells=merge_cells)
        reader = ExcelFile(self.path)
        recons = read_excel(reader, 'test1',
                            index_col=0,
                            ).astype(np.int64)
        frame.index.names = ['test']
        tm.assert_frame_equal(frame, recons.astype(bool))

        self.frame.to_excel(self.path,
                            'test1',
                            columns=['A', 'B', 'C', 'D'],
                            index=False, merge_cells=merge_cells)
        # take 'A' and 'B' as indexes (same row as cols 'C', 'D')
        df = self.frame.copy()
        df = df.set_index(['A', 'B'])

        reader = ExcelFile(self.path)
        recons = read_excel(reader, 'test1', index_col=[0, 1])
        tm.assert_frame_equal(df, recons, check_less_precise=True)

    def test_excel_roundtrip_indexname(self, merge_cells, engine, ext):
        df = DataFrame(np.random.randn(10, 4))
        df.index.name = 'foo'

        df.to_excel(self.path, merge_cells=merge_cells)

        xf = ExcelFile(self.path)
        result = read_excel(xf, xf.sheet_names[0],
                            index_col=0)

        tm.assert_frame_equal(result, df)
        assert result.index.name == 'foo'

    def test_excel_roundtrip_datetime(self, merge_cells, engine, ext):
        # datetime.date, not sure what to test here exactly
        tsf = self.tsframe.copy()

        tsf.index = [x.date() for x in self.tsframe.index]
        tsf.to_excel(self.path, 'test1', merge_cells=merge_cells)
        reader = ExcelFile(self.path)
        recons = read_excel(reader, 'test1')
        tm.assert_frame_equal(self.tsframe, recons)

    # GH4133 - excel output format strings
    def test_excel_date_datetime_format(self, merge_cells, engine, ext):
        df = DataFrame([[date(2014, 1, 31),
                         date(1999, 9, 24)],
                        [datetime(1998, 5, 26, 23, 33, 4),
                         datetime(2014, 2, 28, 13, 5, 13)]],
                       index=['DATE', 'DATETIME'], columns=['X', 'Y'])
        df_expected = DataFrame([[datetime(2014, 1, 31),
                                  datetime(1999, 9, 24)],
                                 [datetime(1998, 5, 26, 23, 33, 4),
                                  datetime(2014, 2, 28, 13, 5, 13)]],
                                index=['DATE', 'DATETIME'], columns=['X', 'Y'])

        with ensure_clean(ext) as filename2:
            writer1 = ExcelWriter(self.path)
            writer2 = ExcelWriter(filename2,
                                  date_format='DD.MM.YYYY',
                                  datetime_format='DD.MM.YYYY HH-MM-SS')

            df.to_excel(writer1, 'test1')
            df.to_excel(writer2, 'test1')

            writer1.close()
            writer2.close()

            reader1 = ExcelFile(self.path)
            reader2 = ExcelFile(filename2)

            rs1 = read_excel(reader1, 'test1', index_col=None)
            rs2 = read_excel(reader2, 'test1', index_col=None)

            tm.assert_frame_equal(rs1, rs2)

            # since the reader returns a datetime object for dates, we need
            # to use df_expected to check the result
            tm.assert_frame_equal(rs2, df_expected)

    def test_to_excel_interval_no_labels(self, merge_cells, engine, ext):
        # GH19242 - test writing Interval without labels
        frame = DataFrame(np.random.randint(-10, 10, size=(20, 1)),
                          dtype=np.int64)
        expected = frame.copy()
        frame['new'] = pd.cut(frame[0], 10)
        expected['new'] = pd.cut(expected[0], 10).astype(str)
        frame.to_excel(self.path, 'test1')
        reader = ExcelFile(self.path)
        recons = read_excel(reader, 'test1')
        tm.assert_frame_equal(expected, recons)

    def test_to_excel_interval_labels(self, merge_cells, engine, ext):
        # GH19242 - test writing Interval with labels
        frame = DataFrame(np.random.randint(-10, 10, size=(20, 1)),
                          dtype=np.int64)
        expected = frame.copy()
        intervals = pd.cut(frame[0], 10, labels=['A', 'B', 'C', 'D', 'E',
                                                 'F', 'G', 'H', 'I', 'J'])
        frame['new'] = intervals
        expected['new'] = pd.Series(list(intervals))
        frame.to_excel(self.path, 'test1')
        reader = ExcelFile(self.path)
        recons = read_excel(reader, 'test1')
        tm.assert_frame_equal(expected, recons)

    def test_to_excel_timedelta(self, merge_cells, engine, ext):
        # GH 19242, GH9155 - test writing timedelta to xls
        frame = DataFrame(np.random.randint(-10, 10, size=(20, 1)),
                          columns=['A'],
                          dtype=np.int64
                          )
        expected = frame.copy()
        frame['new'] = frame['A'].apply(lambda x: timedelta(seconds=x))
        expected['new'] = expected['A'].apply(
            lambda x: timedelta(seconds=x).total_seconds() / float(86400))
        frame.to_excel(self.path, 'test1')
        reader = ExcelFile(self.path)
        recons = read_excel(reader, 'test1')
        tm.assert_frame_equal(expected, recons)

    def test_to_excel_periodindex(self, merge_cells, engine, ext):
        frame = self.tsframe
        xp = frame.resample('M', kind='period').mean()

        xp.to_excel(self.path, 'sht1')

        reader = ExcelFile(self.path)
        rs = read_excel(reader, 'sht1', index_col=0)
        tm.assert_frame_equal(xp, rs.to_period('M'))

    def test_to_excel_multiindex(self, merge_cells, engine, ext):
        frame = self.frame
        arrays = np.arange(len(frame.index) * 2).reshape(2, -1)
        new_index = MultiIndex.from_arrays(arrays,
                                           names=['first', 'second'])
        frame.index = new_index

        frame.to_excel(self.path, 'test1', header=False)
        frame.to_excel(self.path, 'test1', columns=['A', 'B'])

        # round trip
        frame.to_excel(self.path, 'test1', merge_cells=merge_cells)
        reader = ExcelFile(self.path)
        df = read_excel(reader, 'test1', index_col=[0, 1])
        tm.assert_frame_equal(frame, df)

    # GH13511
    def test_to_excel_multiindex_nan_label(self, merge_cells, engine, ext):
        frame = pd.DataFrame({'A': [None, 2, 3],
                              'B': [10, 20, 30],
                              'C': np.random.sample(3)})
        frame = frame.set_index(['A', 'B'])

        frame.to_excel(self.path, merge_cells=merge_cells)
        df = read_excel(self.path, index_col=[0, 1])
        tm.assert_frame_equal(frame, df)

    # Test for Issue 11328. If column indices are integers, make
    # sure they are handled correctly for either setting of
    # merge_cells
    def test_to_excel_multiindex_cols(self, merge_cells, engine, ext):
        frame = self.frame
        arrays = np.arange(len(frame.index) * 2).reshape(2, -1)
        new_index = MultiIndex.from_arrays(arrays,
                                           names=['first', 'second'])
        frame.index = new_index

        new_cols_index = MultiIndex.from_tuples([(40, 1), (40, 2),
                                                 (50, 1), (50, 2)])
        frame.columns = new_cols_index
        header = [0, 1]
        if not merge_cells:
            header = 0

        # round trip
        frame.to_excel(self.path, 'test1', merge_cells=merge_cells)
        reader = ExcelFile(self.path)
        df = read_excel(reader, 'test1', header=header,
                        index_col=[0, 1])
        if not merge_cells:
            fm = frame.columns.format(sparsify=False,
                                      adjoin=False, names=False)
            frame.columns = [".".join(map(str, q)) for q in zip(*fm)]
        tm.assert_frame_equal(frame, df)

    def test_to_excel_multiindex_dates(self, merge_cells, engine, ext):
        # try multiindex with dates
        tsframe = self.tsframe.copy()
        new_index = [tsframe.index, np.arange(len(tsframe.index))]
        tsframe.index = MultiIndex.from_arrays(new_index)

        tsframe.index.names = ['time', 'foo']
        tsframe.to_excel(self.path, 'test1', merge_cells=merge_cells)
        reader = ExcelFile(self.path)
        recons = read_excel(reader, 'test1',
                            index_col=[0, 1])

        tm.assert_frame_equal(tsframe, recons)
        assert recons.index.names == ('time', 'foo')

    def test_to_excel_multiindex_no_write_index(self, merge_cells, engine,
                                                ext):
        # Test writing and re-reading a MI witout the index. GH 5616.

        # Initial non-MI frame.
        frame1 = DataFrame({'a': [10, 20], 'b': [30, 40], 'c': [50, 60]})

        # Add a MI.
        frame2 = frame1.copy()
        multi_index = MultiIndex.from_tuples([(70, 80), (90, 100)])
        frame2.index = multi_index

        # Write out to Excel without the index.
        frame2.to_excel(self.path, 'test1', index=False)

        # Read it back in.
        reader = ExcelFile(self.path)
        frame3 = read_excel(reader, 'test1')

        # Test that it is the same as the initial frame.
        tm.assert_frame_equal(frame1, frame3)

    def test_to_excel_float_format(self, merge_cells, engine, ext):
        df = DataFrame([[0.123456, 0.234567, 0.567567],
                        [12.32112, 123123.2, 321321.2]],
                       index=['A', 'B'], columns=['X', 'Y', 'Z'])

        df.to_excel(self.path, 'test1', float_format='%.2f')

        reader = ExcelFile(self.path)
        rs = read_excel(reader, 'test1', index_col=None)
        xp = DataFrame([[0.12, 0.23, 0.57],
                        [12.32, 123123.20, 321321.20]],
                       index=['A', 'B'], columns=['X', 'Y', 'Z'])
        tm.assert_frame_equal(rs, xp)

    def test_to_excel_output_encoding(self, merge_cells, engine, ext):
        # avoid mixed inferred_type
        df = DataFrame([[u'\u0192', u'\u0193', u'\u0194'],
                        [u'\u0195', u'\u0196', u'\u0197']],
                       index=[u'A\u0192', u'B'],
                       columns=[u'X\u0193', u'Y', u'Z'])

        with ensure_clean('__tmp_to_excel_float_format__.' + ext) as filename:
            df.to_excel(filename, sheet_name='TestSheet', encoding='utf8')
            result = read_excel(filename, 'TestSheet', encoding='utf8')
            tm.assert_frame_equal(result, df)

    def test_to_excel_unicode_filename(self, merge_cells, engine, ext):
        with ensure_clean(u('\u0192u.') + ext) as filename:
            try:
                f = open(filename, 'wb')
            except UnicodeEncodeError:
                pytest.skip('no unicode file names on this system')
            else:
                f.close()

            df = DataFrame([[0.123456, 0.234567, 0.567567],
                            [12.32112, 123123.2, 321321.2]],
                           index=['A', 'B'], columns=['X', 'Y', 'Z'])

            df.to_excel(filename, 'test1', float_format='%.2f')

            reader = ExcelFile(filename)
            rs = read_excel(reader, 'test1', index_col=None)
            xp = DataFrame([[0.12, 0.23, 0.57],
                            [12.32, 123123.20, 321321.20]],
                           index=['A', 'B'], columns=['X', 'Y', 'Z'])
            tm.assert_frame_equal(rs, xp)

    # def test_to_excel_header_styling_xls(self, merge_cells, engine, ext):

    #     import StringIO
    #     s = StringIO(
    #     """Date,ticker,type,value
    #     2001-01-01,x,close,12.2
    #     2001-01-01,x,open ,12.1
    #     2001-01-01,y,close,12.2
    #     2001-01-01,y,open ,12.1
    #     2001-02-01,x,close,12.2
    #     2001-02-01,x,open ,12.1
    #     2001-02-01,y,close,12.2
    #     2001-02-01,y,open ,12.1
    #     2001-03-01,x,close,12.2
    #     2001-03-01,x,open ,12.1
    #     2001-03-01,y,close,12.2
    #     2001-03-01,y,open ,12.1""")
    #     df = read_csv(s, parse_dates=["Date"])
    #     pdf = df.pivot_table(values="value", rows=["ticker"],
    #                                          cols=["Date", "type"])

    #     try:
    #         import xlwt
    #         import xlrd
    #     except ImportError:
    #         pytest.skip

    #     filename = '__tmp_to_excel_header_styling_xls__.xls'
    #     pdf.to_excel(filename, 'test1')

    #     wbk = xlrd.open_workbook(filename,
    #                              formatting_info=True)
    #     assert ["test1"] == wbk.sheet_names()
    #     ws = wbk.sheet_by_name('test1')
    #     assert [(0, 1, 5, 7), (0, 1, 3, 5), (0, 1, 1, 3)] == ws.merged_cells
    #     for i in range(0, 2):
    #         for j in range(0, 7):
    #             xfx = ws.cell_xf_index(0, 0)
    #             cell_xf = wbk.xf_list[xfx]
    #             font = wbk.font_list
    #             assert 1 == font[cell_xf.font_index].bold
    #             assert 1 == cell_xf.border.top_line_style
    #             assert 1 == cell_xf.border.right_line_style
    #             assert 1 == cell_xf.border.bottom_line_style
    #             assert 1 == cell_xf.border.left_line_style
    #             assert 2 == cell_xf.alignment.hor_align
    #     os.remove(filename)
    # def test_to_excel_header_styling_xlsx(self, merge_cells, engine, ext):
    #     import StringIO
    #     s = StringIO(
    #     """Date,ticker,type,value
    #     2001-01-01,x,close,12.2
    #     2001-01-01,x,open ,12.1
    #     2001-01-01,y,close,12.2
    #     2001-01-01,y,open ,12.1
    #     2001-02-01,x,close,12.2
    #     2001-02-01,x,open ,12.1
    #     2001-02-01,y,close,12.2
    #     2001-02-01,y,open ,12.1
    #     2001-03-01,x,close,12.2
    #     2001-03-01,x,open ,12.1
    #     2001-03-01,y,close,12.2
    #     2001-03-01,y,open ,12.1""")
    #     df = read_csv(s, parse_dates=["Date"])
    #     pdf = df.pivot_table(values="value", rows=["ticker"],
    #                                          cols=["Date", "type"])
    #     try:
    #         import openpyxl
    #         from openpyxl.cell import get_column_letter
    #     except ImportError:
    #         pytest.skip
    #     if openpyxl.__version__ < '1.6.1':
    #         pytest.skip
    #     # test xlsx_styling
    #     filename = '__tmp_to_excel_header_styling_xlsx__.xlsx'
    #     pdf.to_excel(filename, 'test1')
    #     wbk = openpyxl.load_workbook(filename)
    #     assert ["test1"] == wbk.get_sheet_names()
    #     ws = wbk.get_sheet_by_name('test1')
    #     xlsaddrs = ["%s2" % chr(i) for i in range(ord('A'), ord('H'))]
    #     xlsaddrs += ["A%s" % i for i in range(1, 6)]
    #     xlsaddrs += ["B1", "D1", "F1"]
    #     for xlsaddr in xlsaddrs:
    #         cell = ws.cell(xlsaddr)
    #         assert cell.style.font.bold
    #         assert (openpyxl.style.Border.BORDER_THIN ==
    #                 cell.style.borders.top.border_style)
    #         assert (openpyxl.style.Border.BORDER_THIN ==
    #                 cell.style.borders.right.border_style)
    #         assert (openpyxl.style.Border.BORDER_THIN ==
    #                 cell.style.borders.bottom.border_style)
    #         assert (openpyxl.style.Border.BORDER_THIN ==
    #                 cell.style.borders.left.border_style)
    #         assert (openpyxl.style.Alignment.HORIZONTAL_CENTER ==
    #                 cell.style.alignment.horizontal)
    #     mergedcells_addrs = ["C1", "E1", "G1"]
    #     for maddr in mergedcells_addrs:
    #         assert ws.cell(maddr).merged
    #     os.remove(filename)

    def test_excel_010_hemstring(self, merge_cells, engine, ext):
        if merge_cells:
            pytest.skip('Skip tests for merged MI format.')

        from pandas.util.testing import makeCustomDataframe as mkdf
        # ensure limited functionality in 0.10
        # override of #2370 until sorted out in 0.11

        def roundtrip(df, header=True, parser_hdr=0, index=True):

            df.to_excel(self.path, header=header,
                        merge_cells=merge_cells, index=index)
            xf = ExcelFile(self.path)
            res = read_excel(xf, xf.sheet_names[0], header=parser_hdr)
            return res

        nrows = 5
        ncols = 3
        for use_headers in (True, False):
            for i in range(1, 4):  # row multindex up to nlevel=3
                for j in range(1, 4):  # col ""
                    df = mkdf(nrows, ncols, r_idx_nlevels=i, c_idx_nlevels=j)

                    # this if will be removed once multi column excel writing
                    # is implemented for now fixing #9794
                    if j > 1:
                        with pytest.raises(NotImplementedError):
                            res = roundtrip(df, use_headers, index=False)
                    else:
                        res = roundtrip(df, use_headers)

                    if use_headers:
                        assert res.shape == (nrows, ncols + i)
                    else:
                        # first row taken as columns
                        assert res.shape == (nrows - 1, ncols + i)

                    # no nans
                    for r in range(len(res.index)):
                        for c in range(len(res.columns)):
                            assert res.iloc[r, c] is not np.nan

        res = roundtrip(DataFrame([0]))
        assert res.shape == (1, 1)
        assert res.iloc[0, 0] is not np.nan

        res = roundtrip(DataFrame([0]), False, None)
        assert res.shape == (1, 2)
        assert res.iloc[0, 0] is not np.nan

    def test_excel_010_hemstring_raises_NotImplementedError(self, merge_cells,
                                                            engine, ext):
        # This test was failing only for j>1 and header=False,
        # So I reproduced a simple test.
        if merge_cells:
            pytest.skip('Skip tests for merged MI format.')

        from pandas.util.testing import makeCustomDataframe as mkdf
        # ensure limited functionality in 0.10
        # override of #2370 until sorted out in 0.11

        def roundtrip2(df, header=True, parser_hdr=0, index=True):

            df.to_excel(self.path, header=header,
                        merge_cells=merge_cells, index=index)
            xf = ExcelFile(self.path)
            res = read_excel(xf, xf.sheet_names[0], header=parser_hdr)
            return res

        nrows = 5
        ncols = 3
        j = 2
        i = 1
        df = mkdf(nrows, ncols, r_idx_nlevels=i, c_idx_nlevels=j)
        with pytest.raises(NotImplementedError):
            roundtrip2(df, header=False, index=False)

    def test_duplicated_columns(self, merge_cells, engine, ext):
        # Test for issue #5235
        write_frame = DataFrame([[1, 2, 3], [1, 2, 3], [1, 2, 3]])
        colnames = ['A', 'B', 'B']

        write_frame.columns = colnames
        write_frame.to_excel(self.path, 'test1')

        read_frame = read_excel(self.path, 'test1')
        read_frame.columns = colnames
        tm.assert_frame_equal(write_frame, read_frame)

        # 11007 / #10970
        write_frame = DataFrame([[1, 2, 3, 4], [5, 6, 7, 8]],
                                columns=['A', 'B', 'A', 'B'])
        write_frame.to_excel(self.path, 'test1')
        read_frame = read_excel(self.path, 'test1')
        read_frame.columns = ['A', 'B', 'A', 'B']
        tm.assert_frame_equal(write_frame, read_frame)

        # 10982
        write_frame.to_excel(self.path, 'test1', index=False, header=False)
        read_frame = read_excel(self.path, 'test1', header=None)
        write_frame.columns = [0, 1, 2, 3]
        tm.assert_frame_equal(write_frame, read_frame)

    def test_swapped_columns(self, merge_cells, engine, ext):
        # Test for issue #5427.
        write_frame = DataFrame({'A': [1, 1, 1],
                                 'B': [2, 2, 2]})
        write_frame.to_excel(self.path, 'test1', columns=['B', 'A'])

        read_frame = read_excel(self.path, 'test1', header=0)

        tm.assert_series_equal(write_frame['A'], read_frame['A'])
        tm.assert_series_equal(write_frame['B'], read_frame['B'])

    def test_invalid_columns(self, merge_cells, engine, ext):
        # 10982
        write_frame = DataFrame({'A': [1, 1, 1],
                                 'B': [2, 2, 2]})

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            write_frame.to_excel(self.path, 'test1', columns=['B', 'C'])
        expected = write_frame.reindex(columns=['B', 'C'])
        read_frame = read_excel(self.path, 'test1')
        tm.assert_frame_equal(expected, read_frame)

        with pytest.raises(KeyError):
            write_frame.to_excel(self.path, 'test1', columns=['C', 'D'])

    def test_comment_arg(self, merge_cells, engine, ext):
        # Re issue #18735
        # Test the comment argument functionality to read_excel

        # Create file to read in
        df = DataFrame({'A': ['one', '#one', 'one'],
                        'B': ['two', 'two', '#two']})
        df.to_excel(self.path, 'test_c')

        # Read file without comment arg
        result1 = read_excel(self.path, 'test_c')
        result1.iloc[1, 0] = None
        result1.iloc[1, 1] = None
        result1.iloc[2, 1] = None
        result2 = read_excel(self.path, 'test_c', comment='#')
        tm.assert_frame_equal(result1, result2)

    def test_comment_default(self, merge_cells, engine, ext):
        # Re issue #18735
        # Test the comment argument default to read_excel

        # Create file to read in
        df = DataFrame({'A': ['one', '#one', 'one'],
                        'B': ['two', 'two', '#two']})
        df.to_excel(self.path, 'test_c')

        # Read file with default and explicit comment=None
        result1 = read_excel(self.path, 'test_c')
        result2 = read_excel(self.path, 'test_c', comment=None)
        tm.assert_frame_equal(result1, result2)

    def test_comment_used(self, merge_cells, engine, ext):
        # Re issue #18735
        # Test the comment argument is working as expected when used

        # Create file to read in
        df = DataFrame({'A': ['one', '#one', 'one'],
                        'B': ['two', 'two', '#two']})
        df.to_excel(self.path, 'test_c')

        # Test read_frame_comment against manually produced expected output
        expected = DataFrame({'A': ['one', None, 'one'],
                              'B': ['two', None, None]})
        result = read_excel(self.path, 'test_c', comment='#')
        tm.assert_frame_equal(result, expected)

    def test_comment_emptyline(self, merge_cells, engine, ext):
        # Re issue #18735
        # Test that read_excel ignores commented lines at the end of file

        df = DataFrame({'a': ['1', '#2'], 'b': ['2', '3']})
        df.to_excel(self.path, index=False)

        # Test that all-comment lines at EoF are ignored
        expected = DataFrame({'a': [1], 'b': [2]})
        result = read_excel(self.path, comment='#')
        tm.assert_frame_equal(result, expected)

    def test_datetimes(self, merge_cells, engine, ext):

        # Test writing and reading datetimes. For issue #9139. (xref #9185)
        datetimes = [datetime(2013, 1, 13, 1, 2, 3),
                     datetime(2013, 1, 13, 2, 45, 56),
                     datetime(2013, 1, 13, 4, 29, 49),
                     datetime(2013, 1, 13, 6, 13, 42),
                     datetime(2013, 1, 13, 7, 57, 35),
                     datetime(2013, 1, 13, 9, 41, 28),
                     datetime(2013, 1, 13, 11, 25, 21),
                     datetime(2013, 1, 13, 13, 9, 14),
                     datetime(2013, 1, 13, 14, 53, 7),
                     datetime(2013, 1, 13, 16, 37, 0),
                     datetime(2013, 1, 13, 18, 20, 52)]

        write_frame = DataFrame({'A': datetimes})
        write_frame.to_excel(self.path, 'Sheet1')
        read_frame = read_excel(self.path, 'Sheet1', header=0)

        tm.assert_series_equal(write_frame['A'], read_frame['A'])

    # GH7074
    def test_bytes_io(self, merge_cells, engine, ext):
        bio = BytesIO()
        df = DataFrame(np.random.randn(10, 2))
        # pass engine explicitly as there is no file path to infer from
        writer = ExcelWriter(bio, engine=engine)
        df.to_excel(writer)
        writer.save()
        bio.seek(0)
        reread_df = read_excel(bio)
        tm.assert_frame_equal(df, reread_df)

    # GH8188
    def test_write_lists_dict(self, merge_cells, engine, ext):
        df = DataFrame({'mixed': ['a', ['b', 'c'], {'d': 'e', 'f': 2}],
                        'numeric': [1, 2, 3.0],
                        'str': ['apple', 'banana', 'cherry']})
        expected = df.copy()
        expected.mixed = expected.mixed.apply(str)
        expected.numeric = expected.numeric.astype('int64')

        df.to_excel(self.path, 'Sheet1')
        read = read_excel(self.path, 'Sheet1', header=0)
        tm.assert_frame_equal(read, expected)

    # GH13347
    def test_true_and_false_value_options(self, merge_cells, engine, ext):
        df = pd.DataFrame([['foo', 'bar']], columns=['col1', 'col2'])
        expected = df.replace({'foo': True,
                               'bar': False})

        df.to_excel(self.path)
        read_frame = read_excel(self.path, true_values=['foo'],
                                false_values=['bar'])
        tm.assert_frame_equal(read_frame, expected)

    def test_freeze_panes(self, merge_cells, engine, ext):
        # GH15160
        expected = DataFrame([[1, 2], [3, 4]], columns=['col1', 'col2'])
        expected.to_excel(self.path, "Sheet1", freeze_panes=(1, 1))
        result = read_excel(self.path)
        tm.assert_frame_equal(expected, result)

    def test_path_pathlib(self, merge_cells, engine, ext):
        df = tm.makeDataFrame()
        writer = partial(df.to_excel, engine=engine)
        reader = partial(pd.read_excel)
        result = tm.round_trip_pathlib(writer, reader,
                                       path="foo.{}".format(ext))
        tm.assert_frame_equal(df, result)

    def test_path_localpath(self, merge_cells, engine, ext):
        df = tm.makeDataFrame()
        writer = partial(df.to_excel, engine=engine)
        reader = partial(pd.read_excel)
        result = tm.round_trip_pathlib(writer, reader,
                                       path="foo.{}".format(ext))
        tm.assert_frame_equal(df, result)


@td.skip_if_no('openpyxl')
@pytest.mark.parametrize("merge_cells,ext,engine", [
    (None, '.xlsx', 'openpyxl')])
class TestOpenpyxlTests(_WriterBase):

    def test_to_excel_styleconverter(self, merge_cells, ext, engine):
        from openpyxl import styles

        hstyle = {
            "font": {
                "color": '00FF0000',
                "bold": True,
            },
            "borders": {
                "top": "thin",
                "right": "thin",
                "bottom": "thin",
                "left": "thin",
            },
            "alignment": {
                "horizontal": "center",
                "vertical": "top",
            },
            "fill": {
                "patternType": 'solid',
                'fgColor': {
                    'rgb': '006666FF',
                    'tint': 0.3,
                },
            },
            "number_format": {
                "format_code": "0.00"
            },
            "protection": {
                "locked": True,
                "hidden": False,
            },
        }

        font_color = styles.Color('00FF0000')
        font = styles.Font(bold=True, color=font_color)
        side = styles.Side(style=styles.borders.BORDER_THIN)
        border = styles.Border(top=side, right=side, bottom=side, left=side)
        alignment = styles.Alignment(horizontal='center', vertical='top')
        fill_color = styles.Color(rgb='006666FF', tint=0.3)
        fill = styles.PatternFill(patternType='solid', fgColor=fill_color)

        number_format = '0.00'

        protection = styles.Protection(locked=True, hidden=False)

        kw = _OpenpyxlWriter._convert_to_style_kwargs(hstyle)
        assert kw['font'] == font
        assert kw['border'] == border
        assert kw['alignment'] == alignment
        assert kw['fill'] == fill
        assert kw['number_format'] == number_format
        assert kw['protection'] == protection

    def test_write_cells_merge_styled(self, merge_cells, ext, engine):
        from pandas.io.formats.excel import ExcelCell

        sheet_name = 'merge_styled'

        sty_b1 = {'font': {'color': '00FF0000'}}
        sty_a2 = {'font': {'color': '0000FF00'}}

        initial_cells = [
            ExcelCell(col=1, row=0, val=42, style=sty_b1),
            ExcelCell(col=0, row=1, val=99, style=sty_a2),
        ]

        sty_merged = {'font': {'color': '000000FF', 'bold': True}}
        sty_kwargs = _OpenpyxlWriter._convert_to_style_kwargs(sty_merged)
        openpyxl_sty_merged = sty_kwargs['font']
        merge_cells = [
            ExcelCell(col=0, row=0, val='pandas',
                      mergestart=1, mergeend=1, style=sty_merged),
        ]

        with ensure_clean(ext) as path:
            writer = _OpenpyxlWriter(path)
            writer.write_cells(initial_cells, sheet_name=sheet_name)
            writer.write_cells(merge_cells, sheet_name=sheet_name)

            wks = writer.sheets[sheet_name]
            xcell_b1 = wks['B1']
            xcell_a2 = wks['A2']
            assert xcell_b1.font == openpyxl_sty_merged
            assert xcell_a2.font == openpyxl_sty_merged


@td.skip_if_no('xlwt')
@pytest.mark.parametrize("merge_cells,ext,engine", [
    (None, '.xls', 'xlwt')])
class TestXlwtTests(_WriterBase):

    def test_excel_raise_error_on_multiindex_columns_and_no_index(
            self, merge_cells, ext, engine):
        # MultiIndex as columns is not yet implemented 9794
        cols = MultiIndex.from_tuples([('site', ''),
                                       ('2014', 'height'),
                                       ('2014', 'weight')])
        df = DataFrame(np.random.randn(10, 3), columns=cols)
        with pytest.raises(NotImplementedError):
            with ensure_clean(ext) as path:
                df.to_excel(path, index=False)

    def test_excel_multiindex_columns_and_index_true(self, merge_cells, ext,
                                                     engine):
        cols = MultiIndex.from_tuples([('site', ''),
                                       ('2014', 'height'),
                                       ('2014', 'weight')])
        df = pd.DataFrame(np.random.randn(10, 3), columns=cols)
        with ensure_clean(ext) as path:
            df.to_excel(path, index=True)

    def test_excel_multiindex_index(self, merge_cells, ext, engine):
        # MultiIndex as index works so assert no error #9794
        cols = MultiIndex.from_tuples([('site', ''),
                                       ('2014', 'height'),
                                       ('2014', 'weight')])
        df = DataFrame(np.random.randn(3, 10), index=cols)
        with ensure_clean(ext) as path:
            df.to_excel(path, index=False)

    def test_to_excel_styleconverter(self, merge_cells, ext, engine):
        import xlwt

        hstyle = {"font": {"bold": True},
                  "borders": {"top": "thin",
                              "right": "thin",
                              "bottom": "thin",
                              "left": "thin"},
                  "alignment": {"horizontal": "center", "vertical": "top"}}

        xls_style = _XlwtWriter._convert_to_style(hstyle)
        assert xls_style.font.bold
        assert xlwt.Borders.THIN == xls_style.borders.top
        assert xlwt.Borders.THIN == xls_style.borders.right
        assert xlwt.Borders.THIN == xls_style.borders.bottom
        assert xlwt.Borders.THIN == xls_style.borders.left
        assert xlwt.Alignment.HORZ_CENTER == xls_style.alignment.horz
        assert xlwt.Alignment.VERT_TOP == xls_style.alignment.vert


@td.skip_if_no('xlsxwriter')
@pytest.mark.parametrize("merge_cells,ext,engine", [
    (None, '.xlsx', 'xlsxwriter')])
class TestXlsxWriterTests(_WriterBase):

    @td.skip_if_no('openpyxl')
    def test_column_format(self, merge_cells, ext, engine):
        # Test that column formats are applied to cells. Test for issue #9167.
        # Applicable to xlsxwriter only.
        with warnings.catch_warnings():
            # Ignore the openpyxl lxml warning.
            warnings.simplefilter("ignore")
            import openpyxl

        with ensure_clean(ext) as path:
            frame = DataFrame({'A': [123456, 123456],
                               'B': [123456, 123456]})

            writer = ExcelWriter(path)
            frame.to_excel(writer)

            # Add a number format to col B and ensure it is applied to cells.
            num_format = '#,##0'
            write_workbook = writer.book
            write_worksheet = write_workbook.worksheets()[0]
            col_format = write_workbook.add_format({'num_format': num_format})
            write_worksheet.set_column('B:B', None, col_format)
            writer.save()

            read_workbook = openpyxl.load_workbook(path)
            try:
                read_worksheet = read_workbook['Sheet1']
            except TypeError:
                # compat
                read_worksheet = read_workbook.get_sheet_by_name(name='Sheet1')

            # Get the number format from the cell.
            try:
                cell = read_worksheet['B2']
            except TypeError:
                # compat
                cell = read_worksheet.cell('B2')

            try:
                read_num_format = cell.number_format
            except Exception:
                read_num_format = cell.style.number_format._format_code

            assert read_num_format == num_format


class TestExcelWriterEngineTests(object):

    @pytest.mark.parametrize('klass,ext', [
        pytest.param(_XlsxWriter, '.xlsx', marks=pytest.mark.skipif(
            not td.safe_import('xlsxwriter'), reason='No xlsxwriter')),
        pytest.param(_OpenpyxlWriter, '.xlsx', marks=pytest.mark.skipif(
            not td.safe_import('openpyxl'), reason='No openpyxl')),
        pytest.param(_XlwtWriter, '.xls', marks=pytest.mark.skipif(
            not td.safe_import('xlwt'), reason='No xlwt'))
    ])
    def test_ExcelWriter_dispatch(self, klass, ext):
        with ensure_clean(ext) as path:
            writer = ExcelWriter(path)
            if ext == '.xlsx' and td.safe_import('xlsxwriter'):
                # xlsxwriter has preference over openpyxl if both installed
                assert isinstance(writer, _XlsxWriter)
            else:
                assert isinstance(writer, klass)

    def test_ExcelWriter_dispatch_raises(self):
        with tm.assert_raises_regex(ValueError, 'No engine'):
            ExcelWriter('nothing')

    def test_register_writer(self):
        # some awkward mocking to test out dispatch and such actually works
        called_save = []
        called_write_cells = []

        class DummyClass(ExcelWriter):
            called_save = False
            called_write_cells = False
            supported_extensions = ['test', 'xlsx', 'xls']
            engine = 'dummy'

            def save(self):
                called_save.append(True)

            def write_cells(self, *args, **kwargs):
                called_write_cells.append(True)

        def check_called(func):
            func()
            assert len(called_save) >= 1
            assert len(called_write_cells) >= 1
            del called_save[:]
            del called_write_cells[:]

        with pd.option_context('io.excel.xlsx.writer', 'dummy'):
            register_writer(DummyClass)
            writer = ExcelWriter('something.test')
            assert isinstance(writer, DummyClass)
            df = tm.makeCustomDataframe(1, 1)

            with catch_warnings(record=True):
                panel = tm.makePanel()
                func = lambda: df.to_excel('something.test')
                check_called(func)
                check_called(lambda: panel.to_excel('something.test'))
                check_called(lambda: df.to_excel('something.xlsx'))
                check_called(
                    lambda: df.to_excel(
                        'something.xls', engine='dummy'))


@pytest.mark.parametrize('engine', [
    pytest.param('xlwt',
                 marks=pytest.mark.xfail(reason='xlwt does not support '
                                                'openpyxl-compatible '
                                                'style dicts')),
    'xlsxwriter',
    'openpyxl',
])
def test_styler_to_excel(engine):
    def style(df):
        # XXX: RGB colors not supported in xlwt
        return DataFrame([['font-weight: bold', '', ''],
                          ['', 'color: blue', ''],
                          ['', '', 'text-decoration: underline'],
                          ['border-style: solid', '', ''],
                          ['', 'font-style: italic', ''],
                          ['', '', 'text-align: right'],
                          ['background-color: red', '', ''],
                          ['', '', ''],
                          ['', '', ''],
                          ['', '', '']],
                         index=df.index, columns=df.columns)

    def assert_equal_style(cell1, cell2):
        # XXX: should find a better way to check equality
        assert cell1.alignment.__dict__ == cell2.alignment.__dict__
        assert cell1.border.__dict__ == cell2.border.__dict__
        assert cell1.fill.__dict__ == cell2.fill.__dict__
        assert cell1.font.__dict__ == cell2.font.__dict__
        assert cell1.number_format == cell2.number_format
        assert cell1.protection.__dict__ == cell2.protection.__dict__

    def custom_converter(css):
        # use bold iff there is custom style attached to the cell
        if css.strip(' \n;'):
            return {'font': {'bold': True}}
        return {}

    pytest.importorskip('jinja2')
    pytest.importorskip(engine)

    # Prepare spreadsheets

    df = DataFrame(np.random.randn(10, 3))
    with ensure_clean('.xlsx' if engine != 'xlwt' else '.xls') as path:
        writer = ExcelWriter(path, engine=engine)
        df.to_excel(writer, sheet_name='frame')
        df.style.to_excel(writer, sheet_name='unstyled')
        styled = df.style.apply(style, axis=None)
        styled.to_excel(writer, sheet_name='styled')
        ExcelFormatter(styled, style_converter=custom_converter).write(
            writer, sheet_name='custom')
        writer.save()

        if engine not in ('openpyxl', 'xlsxwriter'):
            # For other engines, we only smoke test
            return
        openpyxl = pytest.importorskip('openpyxl')
        wb = openpyxl.load_workbook(path)

        # (1) compare DataFrame.to_excel and Styler.to_excel when unstyled
        n_cells = 0
        for col1, col2 in zip(wb['frame'].columns,
                              wb['unstyled'].columns):
            assert len(col1) == len(col2)
            for cell1, cell2 in zip(col1, col2):
                assert cell1.value == cell2.value
                assert_equal_style(cell1, cell2)
                n_cells += 1

        # ensure iteration actually happened:
        assert n_cells == (10 + 1) * (3 + 1)

        # (2) check styling with default converter

        # XXX: openpyxl (as at 2.4) prefixes colors with 00, xlsxwriter with FF
        alpha = '00' if engine == 'openpyxl' else 'FF'

        n_cells = 0
        for col1, col2 in zip(wb['frame'].columns,
                              wb['styled'].columns):
            assert len(col1) == len(col2)
            for cell1, cell2 in zip(col1, col2):
                ref = '%s%d' % (cell2.column, cell2.row)
                # XXX: this isn't as strong a test as ideal; we should
                #      confirm that differences are exclusive
                if ref == 'B2':
                    assert not cell1.font.bold
                    assert cell2.font.bold
                elif ref == 'C3':
                    assert cell1.font.color.rgb != cell2.font.color.rgb
                    assert cell2.font.color.rgb == alpha + '0000FF'
                elif ref == 'D4':
                    # This fails with engine=xlsxwriter due to
                    # https://bitbucket.org/openpyxl/openpyxl/issues/800
                    if engine == 'xlsxwriter' \
                       and (LooseVersion(openpyxl.__version__) <
                            LooseVersion('2.4.6')):
                        pass
                    else:
                        assert cell1.font.underline != cell2.font.underline
                        assert cell2.font.underline == 'single'
                elif ref == 'B5':
                    assert not cell1.border.left.style
                    assert (cell2.border.top.style ==
                            cell2.border.right.style ==
                            cell2.border.bottom.style ==
                            cell2.border.left.style ==
                            'medium')
                elif ref == 'C6':
                    assert not cell1.font.italic
                    assert cell2.font.italic
                elif ref == 'D7':
                    assert (cell1.alignment.horizontal !=
                            cell2.alignment.horizontal)
                    assert cell2.alignment.horizontal == 'right'
                elif ref == 'B8':
                    assert cell1.fill.fgColor.rgb != cell2.fill.fgColor.rgb
                    assert cell1.fill.patternType != cell2.fill.patternType
                    assert cell2.fill.fgColor.rgb == alpha + 'FF0000'
                    assert cell2.fill.patternType == 'solid'
                else:
                    assert_equal_style(cell1, cell2)

                assert cell1.value == cell2.value
                n_cells += 1

        assert n_cells == (10 + 1) * (3 + 1)

        # (3) check styling with custom converter
        n_cells = 0
        for col1, col2 in zip(wb['frame'].columns,
                              wb['custom'].columns):
            assert len(col1) == len(col2)
            for cell1, cell2 in zip(col1, col2):
                ref = '%s%d' % (cell2.column, cell2.row)
                if ref in ('B2', 'C3', 'D4', 'B5', 'C6', 'D7', 'B8'):
                    assert not cell1.font.bold
                    assert cell2.font.bold
                else:
                    assert_equal_style(cell1, cell2)

                assert cell1.value == cell2.value
                n_cells += 1

        assert n_cells == (10 + 1) * (3 + 1)


@td.skip_if_no('openpyxl')
@pytest.mark.skipif(not PY36, reason='requires fspath')
class TestFSPath(object):

    def test_excelfile_fspath(self):
        with tm.ensure_clean('foo.xlsx') as path:
            df = DataFrame({"A": [1, 2]})
            df.to_excel(path)
            xl = ExcelFile(path)
            result = os.fspath(xl)
            assert result == path

    def test_excelwriter_fspath(self):
        with tm.ensure_clean('foo.xlsx') as path:
            writer = ExcelWriter(path)
            assert os.fspath(writer) == str(path)
