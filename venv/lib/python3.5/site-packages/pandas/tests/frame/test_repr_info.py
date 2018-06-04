# -*- coding: utf-8 -*-

from __future__ import print_function

from datetime import datetime, timedelta
import re
import sys
import textwrap

from numpy import nan
import numpy as np
import pytest

from pandas import (DataFrame, Series, compat, option_context,
                    date_range, period_range, Categorical)
from pandas.compat import StringIO, lrange, u, PYPY
import pandas.io.formats.format as fmt
import pandas as pd

import pandas.util.testing as tm

from pandas.tests.frame.common import TestData


# Segregated collection of methods that require the BlockManager internal data
# structure


class TestDataFrameReprInfoEtc(TestData):

    def test_repr_empty(self):
        # empty
        foo = repr(self.empty)  # noqa

        # empty with index
        frame = DataFrame(index=np.arange(1000))
        foo = repr(frame)  # noqa

    def test_repr_mixed(self):
        buf = StringIO()

        # mixed
        foo = repr(self.mixed_frame)  # noqa
        self.mixed_frame.info(verbose=False, buf=buf)

    @pytest.mark.slow
    def test_repr_mixed_big(self):
        # big mixed
        biggie = DataFrame({'A': np.random.randn(200),
                            'B': tm.makeStringIndex(200)},
                           index=lrange(200))
        biggie.loc[:20, 'A'] = nan
        biggie.loc[:20, 'B'] = nan

        foo = repr(biggie)  # noqa

    def test_repr(self):
        buf = StringIO()

        # small one
        foo = repr(self.frame)
        self.frame.info(verbose=False, buf=buf)

        # even smaller
        self.frame.reindex(columns=['A']).info(verbose=False, buf=buf)
        self.frame.reindex(columns=['A', 'B']).info(verbose=False, buf=buf)

        # exhausting cases in DataFrame.info

        # columns but no index
        no_index = DataFrame(columns=[0, 1, 3])
        foo = repr(no_index)  # noqa

        # no columns or index
        self.empty.info(buf=buf)

        df = DataFrame(["a\n\r\tb"], columns=["a\n\r\td"], index=["a\n\r\tf"])
        assert "\t" not in repr(df)
        assert "\r" not in repr(df)
        assert "a\n" not in repr(df)

    def test_repr_dimensions(self):
        df = DataFrame([[1, 2, ], [3, 4]])
        with option_context('display.show_dimensions', True):
            assert "2 rows x 2 columns" in repr(df)

        with option_context('display.show_dimensions', False):
            assert "2 rows x 2 columns" not in repr(df)

        with option_context('display.show_dimensions', 'truncate'):
            assert "2 rows x 2 columns" not in repr(df)

    @pytest.mark.slow
    def test_repr_big(self):
        # big one
        biggie = DataFrame(np.zeros((200, 4)), columns=lrange(4),
                           index=lrange(200))
        repr(biggie)

    def test_repr_unsortable(self):
        # columns are not sortable
        import warnings
        warn_filters = warnings.filters
        warnings.filterwarnings('ignore',
                                category=FutureWarning,
                                module=".*format")

        unsortable = DataFrame({'foo': [1] * 50,
                                datetime.today(): [1] * 50,
                                'bar': ['bar'] * 50,
                                datetime.today() + timedelta(1): ['bar'] * 50},
                               index=np.arange(50))
        repr(unsortable)

        fmt.set_option('display.precision', 3, 'display.column_space', 10)
        repr(self.frame)

        fmt.set_option('display.max_rows', 10, 'display.max_columns', 2)
        repr(self.frame)

        fmt.set_option('display.max_rows', 1000, 'display.max_columns', 1000)
        repr(self.frame)

        tm.reset_display_options()

        warnings.filters = warn_filters

    def test_repr_unicode(self):
        uval = u('\u03c3\u03c3\u03c3\u03c3')

        # TODO(wesm): is this supposed to be used?
        bval = uval.encode('utf-8')  # noqa

        df = DataFrame({'A': [uval, uval]})

        result = repr(df)
        ex_top = '      A'
        assert result.split('\n')[0].rstrip() == ex_top

        df = DataFrame({'A': [uval, uval]})
        result = repr(df)
        assert result.split('\n')[0].rstrip() == ex_top

    def test_unicode_string_with_unicode(self):
        df = DataFrame({'A': [u("\u05d0")]})

        if compat.PY3:
            str(df)
        else:
            compat.text_type(df)

    def test_bytestring_with_unicode(self):
        df = DataFrame({'A': [u("\u05d0")]})
        if compat.PY3:
            bytes(df)
        else:
            str(df)

    def test_very_wide_info_repr(self):
        df = DataFrame(np.random.randn(10, 20),
                       columns=tm.rands_array(10, 20))
        repr(df)

    def test_repr_column_name_unicode_truncation_bug(self):
        # #1906
        df = DataFrame({'Id': [7117434],
                        'StringCol': ('Is it possible to modify drop plot code'
                                      ' so that the output graph is displayed '
                                      'in iphone simulator, Is it possible to '
                                      'modify drop plot code so that the '
                                      'output graph is \xe2\x80\xa8displayed '
                                      'in iphone simulator.Now we are adding '
                                      'the CSV file externally. I want to Call'
                                      ' the File through the code..')})

        with option_context('display.max_columns', 20):
            assert 'StringCol' in repr(df)

    def test_latex_repr(self):
        result = r"""\begin{tabular}{llll}
\toprule
{} &         0 &  1 &  2 \\
\midrule
0 &  $\alpha$ &  b &  c \\
1 &         1 &  2 &  3 \\
\bottomrule
\end{tabular}
"""
        with option_context("display.latex.escape", False,
                            'display.latex.repr', True):
            df = DataFrame([[r'$\alpha$', 'b', 'c'], [1, 2, 3]])
            assert result == df._repr_latex_()

        # GH 12182
        assert df._repr_latex_() is None

    @tm.capture_stdout
    def test_info(self):
        io = StringIO()
        self.frame.info(buf=io)
        self.tsframe.info(buf=io)

        frame = DataFrame(np.random.randn(5, 3))

        frame.info()
        frame.info(verbose=False)

    def test_info_memory(self):
        # https://github.com/pandas-dev/pandas/issues/21056
        df = pd.DataFrame({'a': pd.Series([1, 2], dtype='i8')})
        buf = StringIO()
        df.info(buf=buf)
        result = buf.getvalue()
        bytes = float(df.memory_usage().sum())

        expected = textwrap.dedent("""\
        <class 'pandas.core.frame.DataFrame'>
        RangeIndex: 2 entries, 0 to 1
        Data columns (total 1 columns):
        a    2 non-null int64
        dtypes: int64(1)
        memory usage: {} bytes
        """.format(bytes))

        assert result == expected

    def test_info_wide(self):
        from pandas import set_option, reset_option
        io = StringIO()
        df = DataFrame(np.random.randn(5, 101))
        df.info(buf=io)

        io = StringIO()
        df.info(buf=io, max_cols=101)
        rs = io.getvalue()
        assert len(rs.splitlines()) > 100
        xp = rs

        set_option('display.max_info_columns', 101)
        io = StringIO()
        df.info(buf=io)
        assert rs == xp
        reset_option('display.max_info_columns')

    def test_info_duplicate_columns(self):
        io = StringIO()

        # it works!
        frame = DataFrame(np.random.randn(1500, 4),
                          columns=['a', 'a', 'b', 'b'])
        frame.info(buf=io)

    def test_info_duplicate_columns_shows_correct_dtypes(self):
        # GH11761
        io = StringIO()

        frame = DataFrame([[1, 2.0]],
                          columns=['a', 'a'])
        frame.info(buf=io)
        io.seek(0)
        lines = io.readlines()
        assert 'a    1 non-null int64\n' == lines[3]
        assert 'a    1 non-null float64\n' == lines[4]

    def test_info_shows_column_dtypes(self):
        dtypes = ['int64', 'float64', 'datetime64[ns]', 'timedelta64[ns]',
                  'complex128', 'object', 'bool']
        data = {}
        n = 10
        for i, dtype in enumerate(dtypes):
            data[i] = np.random.randint(2, size=n).astype(dtype)
        df = DataFrame(data)
        buf = StringIO()
        df.info(buf=buf)
        res = buf.getvalue()
        for i, dtype in enumerate(dtypes):
            name = '%d    %d non-null %s' % (i, n, dtype)
            assert name in res

    def test_info_max_cols(self):
        df = DataFrame(np.random.randn(10, 5))
        for len_, verbose in [(5, None), (5, False), (10, True)]:
            # For verbose always      ^ setting  ^ summarize ^ full output
            with option_context('max_info_columns', 4):
                buf = StringIO()
                df.info(buf=buf, verbose=verbose)
                res = buf.getvalue()
                assert len(res.strip().split('\n')) == len_

        for len_, verbose in [(10, None), (5, False), (10, True)]:

            # max_cols no exceeded
            with option_context('max_info_columns', 5):
                buf = StringIO()
                df.info(buf=buf, verbose=verbose)
                res = buf.getvalue()
                assert len(res.strip().split('\n')) == len_

        for len_, max_cols in [(10, 5), (5, 4)]:
            # setting truncates
            with option_context('max_info_columns', 4):
                buf = StringIO()
                df.info(buf=buf, max_cols=max_cols)
                res = buf.getvalue()
                assert len(res.strip().split('\n')) == len_

            # setting wouldn't truncate
            with option_context('max_info_columns', 5):
                buf = StringIO()
                df.info(buf=buf, max_cols=max_cols)
                res = buf.getvalue()
                assert len(res.strip().split('\n')) == len_

    def test_info_memory_usage(self):
        # Ensure memory usage is displayed, when asserted, on the last line
        dtypes = ['int64', 'float64', 'datetime64[ns]', 'timedelta64[ns]',
                  'complex128', 'object', 'bool']
        data = {}
        n = 10
        for i, dtype in enumerate(dtypes):
            data[i] = np.random.randint(2, size=n).astype(dtype)
        df = DataFrame(data)
        buf = StringIO()

        # display memory usage case
        df.info(buf=buf, memory_usage=True)
        res = buf.getvalue().splitlines()
        assert "memory usage: " in res[-1]

        # do not display memory usage case
        df.info(buf=buf, memory_usage=False)
        res = buf.getvalue().splitlines()
        assert "memory usage: " not in res[-1]

        df.info(buf=buf, memory_usage=True)
        res = buf.getvalue().splitlines()

        # memory usage is a lower bound, so print it as XYZ+ MB
        assert re.match(r"memory usage: [^+]+\+", res[-1])

        df.iloc[:, :5].info(buf=buf, memory_usage=True)
        res = buf.getvalue().splitlines()

        # excluded column with object dtype, so estimate is accurate
        assert not re.match(r"memory usage: [^+]+\+", res[-1])

        # Test a DataFrame with duplicate columns
        dtypes = ['int64', 'int64', 'int64', 'float64']
        data = {}
        n = 100
        for i, dtype in enumerate(dtypes):
            data[i] = np.random.randint(2, size=n).astype(dtype)
        df = DataFrame(data)
        df.columns = dtypes

        df_with_object_index = pd.DataFrame({'a': [1]}, index=['foo'])
        df_with_object_index.info(buf=buf, memory_usage=True)
        res = buf.getvalue().splitlines()
        assert re.match(r"memory usage: [^+]+\+", res[-1])

        df_with_object_index.info(buf=buf, memory_usage='deep')
        res = buf.getvalue().splitlines()
        assert re.match(r"memory usage: [^+]+$", res[-1])

        # Ensure df size is as expected
        # (cols * rows * bytes) + index size
        df_size = df.memory_usage().sum()
        exp_size = len(dtypes) * n * 8 + df.index.nbytes
        assert df_size == exp_size

        # Ensure number of cols in memory_usage is the same as df
        size_df = np.size(df.columns.values) + 1  # index=True; default
        assert size_df == np.size(df.memory_usage())

        # assert deep works only on object
        assert df.memory_usage().sum() == df.memory_usage(deep=True).sum()

        # test for validity
        DataFrame(1, index=['a'], columns=['A']
                  ).memory_usage(index=True)
        DataFrame(1, index=['a'], columns=['A']
                  ).index.nbytes
        df = DataFrame(
            data=1,
            index=pd.MultiIndex.from_product(
                [['a'], range(1000)]),
            columns=['A']
        )
        df.index.nbytes
        df.memory_usage(index=True)
        df.index.values.nbytes

        mem = df.memory_usage(deep=True).sum()
        assert mem > 0

    @pytest.mark.skipif(PYPY,
                        reason="on PyPy deep=True doesn't change result")
    def test_info_memory_usage_deep_not_pypy(self):
        df_with_object_index = pd.DataFrame({'a': [1]}, index=['foo'])
        assert (df_with_object_index.memory_usage(
                index=True, deep=True).sum() >
                df_with_object_index.memory_usage(
                    index=True).sum())

        df_object = pd.DataFrame({'a': ['a']})
        assert (df_object.memory_usage(deep=True).sum() >
                df_object.memory_usage().sum())

    @pytest.mark.skipif(not PYPY,
                        reason="on PyPy deep=True does not change result")
    def test_info_memory_usage_deep_pypy(self):
        df_with_object_index = pd.DataFrame({'a': [1]}, index=['foo'])
        assert (df_with_object_index.memory_usage(
                index=True, deep=True).sum() ==
                df_with_object_index.memory_usage(
                    index=True).sum())

        df_object = pd.DataFrame({'a': ['a']})
        assert (df_object.memory_usage(deep=True).sum() ==
                df_object.memory_usage().sum())

    @pytest.mark.skipif(PYPY, reason="PyPy getsizeof() fails by design")
    def test_usage_via_getsizeof(self):
        df = DataFrame(
            data=1,
            index=pd.MultiIndex.from_product(
                [['a'], range(1000)]),
            columns=['A']
        )
        mem = df.memory_usage(deep=True).sum()
        # sys.getsizeof will call the .memory_usage with
        # deep=True, and add on some GC overhead
        diff = mem - sys.getsizeof(df)
        assert abs(diff) < 100

    def test_info_memory_usage_qualified(self):

        buf = StringIO()
        df = DataFrame(1, columns=list('ab'),
                       index=[1, 2, 3])
        df.info(buf=buf)
        assert '+' not in buf.getvalue()

        buf = StringIO()
        df = DataFrame(1, columns=list('ab'),
                       index=list('ABC'))
        df.info(buf=buf)
        assert '+' in buf.getvalue()

        buf = StringIO()
        df = DataFrame(1, columns=list('ab'),
                       index=pd.MultiIndex.from_product(
                           [range(3), range(3)]))
        df.info(buf=buf)
        assert '+' not in buf.getvalue()

        buf = StringIO()
        df = DataFrame(1, columns=list('ab'),
                       index=pd.MultiIndex.from_product(
                           [range(3), ['foo', 'bar']]))
        df.info(buf=buf)
        assert '+' in buf.getvalue()

    def test_info_memory_usage_bug_on_multiindex(self):
        # GH 14308
        # memory usage introspection should not materialize .values

        from string import ascii_uppercase as uppercase

        def memory_usage(f):
            return f.memory_usage(deep=True).sum()

        N = 100
        M = len(uppercase)
        index = pd.MultiIndex.from_product([list(uppercase),
                                            pd.date_range('20160101',
                                                          periods=N)],
                                           names=['id', 'date'])
        df = DataFrame({'value': np.random.randn(N * M)}, index=index)

        unstacked = df.unstack('id')
        assert df.values.nbytes == unstacked.values.nbytes
        assert memory_usage(df) > memory_usage(unstacked)

        # high upper bound
        assert memory_usage(unstacked) - memory_usage(df) < 2000

    def test_info_categorical(self):
        # GH14298
        idx = pd.CategoricalIndex(['a', 'b'])
        df = pd.DataFrame(np.zeros((2, 2)), index=idx, columns=idx)

        buf = StringIO()
        df.info(buf=buf)

    def test_info_categorical_column(self):

        # make sure it works
        n = 2500
        df = DataFrame({'int64': np.random.randint(100, size=n)})
        df['category'] = Series(np.array(list('abcdefghij')).take(
            np.random.randint(0, 10, size=n))).astype('category')
        df.isna()
        buf = StringIO()
        df.info(buf=buf)

        df2 = df[df['category'] == 'd']
        buf = compat.StringIO()
        df2.info(buf=buf)

    def test_repr_categorical_dates_periods(self):
        # normal DataFrame
        dt = date_range('2011-01-01 09:00', freq='H', periods=5,
                        tz='US/Eastern')
        p = period_range('2011-01', freq='M', periods=5)
        df = DataFrame({'dt': dt, 'p': p})
        exp = """                         dt       p
0 2011-01-01 09:00:00-05:00 2011-01
1 2011-01-01 10:00:00-05:00 2011-02
2 2011-01-01 11:00:00-05:00 2011-03
3 2011-01-01 12:00:00-05:00 2011-04
4 2011-01-01 13:00:00-05:00 2011-05"""

        df = DataFrame({'dt': Categorical(dt), 'p': Categorical(p)})
        assert repr(df) == exp
