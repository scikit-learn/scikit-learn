# -*- coding: utf-8 -*-

"""
Tests dtype specification during parsing
for all of the parsers defined in parsers.py
"""

import pytest

import numpy as np
import pandas as pd
import pandas.util.testing as tm

from pandas import DataFrame, Series, Index, MultiIndex, Categorical
from pandas.compat import StringIO
from pandas.core.dtypes.dtypes import CategoricalDtype
from pandas.errors import ParserWarning


class DtypeTests(object):

    def test_passing_dtype(self):
        # see gh-6607
        df = DataFrame(np.random.rand(5, 2).round(4), columns=list(
            'AB'), index=['1A', '1B', '1C', '1D', '1E'])

        with tm.ensure_clean('__passing_str_as_dtype__.csv') as path:
            df.to_csv(path)

            # see gh-3795: passing 'str' as the dtype
            result = self.read_csv(path, dtype=str, index_col=0)
            expected = df.astype(str)
            tm.assert_frame_equal(result, expected)

            # for parsing, interpret object as str
            result = self.read_csv(path, dtype=object, index_col=0)
            tm.assert_frame_equal(result, expected)

            # we expect all object columns, so need to
            # convert to test for equivalence
            result = result.astype(float)
            tm.assert_frame_equal(result, df)

            # invalid dtype
            pytest.raises(TypeError, self.read_csv, path,
                          dtype={'A': 'foo', 'B': 'float64'},
                          index_col=0)

        # see gh-12048: empty frame
        actual = self.read_csv(StringIO('A,B'), dtype=str)
        expected = DataFrame({'A': [], 'B': []}, index=[], dtype=str)
        tm.assert_frame_equal(actual, expected)

    def test_pass_dtype(self):
        data = """\
one,two
1,2.5
2,3.5
3,4.5
4,5.5"""

        result = self.read_csv(StringIO(data), dtype={'one': 'u1', 1: 'S1'})
        assert result['one'].dtype == 'u1'
        assert result['two'].dtype == 'object'

    def test_categorical_dtype(self):
        # GH 10153
        data = """a,b,c
1,a,3.4
1,a,3.4
2,b,4.5"""
        expected = pd.DataFrame({'a': Categorical(['1', '1', '2']),
                                 'b': Categorical(['a', 'a', 'b']),
                                 'c': Categorical(['3.4', '3.4', '4.5'])})
        actual = self.read_csv(StringIO(data), dtype='category')
        tm.assert_frame_equal(actual, expected)

        actual = self.read_csv(StringIO(data), dtype=CategoricalDtype())
        tm.assert_frame_equal(actual, expected)

        actual = self.read_csv(StringIO(data), dtype={'a': 'category',
                                                      'b': 'category',
                                                      'c': CategoricalDtype()})
        tm.assert_frame_equal(actual, expected)

        actual = self.read_csv(StringIO(data), dtype={'b': 'category'})
        expected = pd.DataFrame({'a': [1, 1, 2],
                                 'b': Categorical(['a', 'a', 'b']),
                                 'c': [3.4, 3.4, 4.5]})
        tm.assert_frame_equal(actual, expected)

        actual = self.read_csv(StringIO(data), dtype={1: 'category'})
        tm.assert_frame_equal(actual, expected)

        # unsorted
        data = """a,b,c
1,b,3.4
1,b,3.4
2,a,4.5"""
        expected = pd.DataFrame({'a': Categorical(['1', '1', '2']),
                                 'b': Categorical(['b', 'b', 'a']),
                                 'c': Categorical(['3.4', '3.4', '4.5'])})
        actual = self.read_csv(StringIO(data), dtype='category')
        tm.assert_frame_equal(actual, expected)

        # missing
        data = """a,b,c
1,b,3.4
1,nan,3.4
2,a,4.5"""
        expected = pd.DataFrame({'a': Categorical(['1', '1', '2']),
                                 'b': Categorical(['b', np.nan, 'a']),
                                 'c': Categorical(['3.4', '3.4', '4.5'])})
        actual = self.read_csv(StringIO(data), dtype='category')
        tm.assert_frame_equal(actual, expected)

    @pytest.mark.slow
    def test_categorical_dtype_high_cardinality_numeric(self):
        # GH 18186
        data = np.sort([str(i) for i in range(524289)])
        expected = DataFrame({'a': Categorical(data, ordered=True)})
        actual = self.read_csv(StringIO('a\n' + '\n'.join(data)),
                               dtype='category')
        actual["a"] = actual["a"].cat.reorder_categories(
            np.sort(actual.a.cat.categories), ordered=True)
        tm.assert_frame_equal(actual, expected)

    def test_categorical_dtype_encoding(self):
        # GH 10153
        pth = tm.get_data_path('unicode_series.csv')
        encoding = 'latin-1'
        expected = self.read_csv(pth, header=None, encoding=encoding)
        expected[1] = Categorical(expected[1])
        actual = self.read_csv(pth, header=None, encoding=encoding,
                               dtype={1: 'category'})
        tm.assert_frame_equal(actual, expected)

        pth = tm.get_data_path('utf16_ex.txt')
        encoding = 'utf-16'
        expected = self.read_table(pth, encoding=encoding)
        expected = expected.apply(Categorical)
        actual = self.read_table(pth, encoding=encoding, dtype='category')
        tm.assert_frame_equal(actual, expected)

    def test_categorical_dtype_chunksize(self):
        # GH 10153
        data = """a,b
1,a
1,b
1,b
2,c"""
        expecteds = [pd.DataFrame({'a': [1, 1],
                                   'b': Categorical(['a', 'b'])}),
                     pd.DataFrame({'a': [1, 2],
                                   'b': Categorical(['b', 'c'])},
                                  index=[2, 3])]
        actuals = self.read_csv(StringIO(data), dtype={'b': 'category'},
                                chunksize=2)

        for actual, expected in zip(actuals, expecteds):
            tm.assert_frame_equal(actual, expected)

    @pytest.mark.parametrize('ordered', [False, True])
    @pytest.mark.parametrize('categories', [
        ['a', 'b', 'c'],
        ['a', 'c', 'b'],
        ['a', 'b', 'c', 'd'],
        ['c', 'b', 'a'],
    ])
    def test_categorical_categoricaldtype(self, categories, ordered):
        data = """a,b
1,a
1,b
1,b
2,c"""
        expected = pd.DataFrame({
            "a": [1, 1, 1, 2],
            "b": Categorical(['a', 'b', 'b', 'c'],
                             categories=categories,
                             ordered=ordered)
        })
        dtype = {"b": CategoricalDtype(categories=categories,
                                       ordered=ordered)}
        result = self.read_csv(StringIO(data), dtype=dtype)
        tm.assert_frame_equal(result, expected)

    def test_categorical_categoricaldtype_unsorted(self):
        data = """a,b
1,a
1,b
1,b
2,c"""
        dtype = CategoricalDtype(['c', 'b', 'a'])
        expected = pd.DataFrame({
            'a': [1, 1, 1, 2],
            'b': Categorical(['a', 'b', 'b', 'c'], categories=['c', 'b', 'a'])
        })
        result = self.read_csv(StringIO(data), dtype={'b': dtype})
        tm.assert_frame_equal(result, expected)

    def test_categoricaldtype_coerces_numeric(self):
        dtype = {'b': CategoricalDtype([1, 2, 3])}
        data = "b\n1\n1\n2\n3"
        expected = pd.DataFrame({'b': Categorical([1, 1, 2, 3])})
        result = self.read_csv(StringIO(data), dtype=dtype)
        tm.assert_frame_equal(result, expected)

    def test_categoricaldtype_coerces_datetime(self):
        dtype = {
            'b': CategoricalDtype(pd.date_range('2017', '2019', freq='AS'))
        }
        data = "b\n2017-01-01\n2018-01-01\n2019-01-01"
        expected = pd.DataFrame({'b': Categorical(dtype['b'].categories)})
        result = self.read_csv(StringIO(data), dtype=dtype)
        tm.assert_frame_equal(result, expected)

        dtype = {
            'b': CategoricalDtype([pd.Timestamp("2014")])
        }
        data = "b\n2014-01-01\n2014-01-01T00:00:00"
        expected = pd.DataFrame({'b': Categorical([pd.Timestamp('2014')] * 2)})
        result = self.read_csv(StringIO(data), dtype=dtype)
        tm.assert_frame_equal(result, expected)

    def test_categoricaldtype_coerces_timedelta(self):
        dtype = {'b': CategoricalDtype(pd.to_timedelta(['1H', '2H', '3H']))}
        data = "b\n1H\n2H\n3H"
        expected = pd.DataFrame({'b': Categorical(dtype['b'].categories)})
        result = self.read_csv(StringIO(data), dtype=dtype)
        tm.assert_frame_equal(result, expected)

    def test_categoricaldtype_unexpected_categories(self):
        dtype = {'b': CategoricalDtype(['a', 'b', 'd', 'e'])}
        data = "b\nd\na\nc\nd"  # Unexpected c
        expected = pd.DataFrame({"b": Categorical(list('dacd'),
                                                  dtype=dtype['b'])})
        result = self.read_csv(StringIO(data), dtype=dtype)
        tm.assert_frame_equal(result, expected)

    def test_categorical_categoricaldtype_chunksize(self):
        # GH 10153
        data = """a,b
1,a
1,b
1,b
2,c"""
        cats = ['a', 'b', 'c']
        expecteds = [pd.DataFrame({'a': [1, 1],
                                   'b': Categorical(['a', 'b'],
                                                    categories=cats)}),
                     pd.DataFrame({'a': [1, 2],
                                   'b': Categorical(['b', 'c'],
                                                    categories=cats)},
                                  index=[2, 3])]
        dtype = CategoricalDtype(cats)
        actuals = self.read_csv(StringIO(data), dtype={'b': dtype},
                                chunksize=2)

        for actual, expected in zip(actuals, expecteds):
            tm.assert_frame_equal(actual, expected)

    def test_empty_pass_dtype(self):
        data = 'one,two'
        result = self.read_csv(StringIO(data), dtype={'one': 'u1'})

        expected = DataFrame({'one': np.empty(0, dtype='u1'),
                              'two': np.empty(0, dtype=np.object)})
        tm.assert_frame_equal(result, expected, check_index_type=False)

    def test_empty_with_index_pass_dtype(self):
        data = 'one,two'
        result = self.read_csv(StringIO(data), index_col=['one'],
                               dtype={'one': 'u1', 1: 'f'})

        expected = DataFrame({'two': np.empty(0, dtype='f')},
                             index=Index([], dtype='u1', name='one'))
        tm.assert_frame_equal(result, expected, check_index_type=False)

    def test_empty_with_multiindex_pass_dtype(self):
        data = 'one,two,three'
        result = self.read_csv(StringIO(data), index_col=['one', 'two'],
                               dtype={'one': 'u1', 1: 'f8'})

        exp_idx = MultiIndex.from_arrays([np.empty(0, dtype='u1'),
                                          np.empty(0, dtype='O')],
                                         names=['one', 'two'])
        expected = DataFrame(
            {'three': np.empty(0, dtype=np.object)}, index=exp_idx)
        tm.assert_frame_equal(result, expected, check_index_type=False)

    def test_empty_with_mangled_column_pass_dtype_by_names(self):
        data = 'one,one'
        result = self.read_csv(StringIO(data), dtype={
            'one': 'u1', 'one.1': 'f'})

        expected = DataFrame(
            {'one': np.empty(0, dtype='u1'), 'one.1': np.empty(0, dtype='f')})
        tm.assert_frame_equal(result, expected, check_index_type=False)

    def test_empty_with_mangled_column_pass_dtype_by_indexes(self):
        data = 'one,one'
        result = self.read_csv(StringIO(data), dtype={0: 'u1', 1: 'f'})

        expected = DataFrame(
            {'one': np.empty(0, dtype='u1'), 'one.1': np.empty(0, dtype='f')})
        tm.assert_frame_equal(result, expected, check_index_type=False)

    def test_empty_with_dup_column_pass_dtype_by_indexes(self):
        # see gh-9424
        expected = pd.concat([Series([], name='one', dtype='u1'),
                              Series([], name='one.1', dtype='f')], axis=1)

        data = 'one,one'
        result = self.read_csv(StringIO(data), dtype={0: 'u1', 1: 'f'})
        tm.assert_frame_equal(result, expected, check_index_type=False)

        with tm.assert_produces_warning(UserWarning, check_stacklevel=False):
            data = ''
            result = self.read_csv(StringIO(data), names=['one', 'one'],
                                   dtype={0: 'u1', 1: 'f'})
            tm.assert_frame_equal(result, expected, check_index_type=False)

    def test_raise_on_passed_int_dtype_with_nas(self):
        # see gh-2631
        data = """YEAR, DOY, a
2001,106380451,10
2001,,11
2001,106380451,67"""
        pytest.raises(ValueError, self.read_csv, StringIO(data),
                      sep=",", skipinitialspace=True,
                      dtype={'DOY': np.int64})

    def test_dtype_with_converter(self):
        data = """a,b
1.1,2.2
1.2,2.3"""
        # dtype spec ignored if converted specified
        with tm.assert_produces_warning(ParserWarning):
            result = self.read_csv(StringIO(data), dtype={'a': 'i8'},
                                   converters={'a': lambda x: str(x)})
        expected = DataFrame({'a': ['1.1', '1.2'], 'b': [2.2, 2.3]})
        tm.assert_frame_equal(result, expected)

    def test_empty_dtype(self):
        # see gh-14712
        data = 'a,b'

        expected = pd.DataFrame(columns=['a', 'b'], dtype=np.float64)
        result = self.read_csv(StringIO(data), header=0, dtype=np.float64)
        tm.assert_frame_equal(result, expected)

        expected = pd.DataFrame({'a': pd.Categorical([]),
                                 'b': pd.Categorical([])},
                                index=[])
        result = self.read_csv(StringIO(data), header=0,
                               dtype='category')
        tm.assert_frame_equal(result, expected)
        result = self.read_csv(StringIO(data), header=0,
                               dtype={'a': 'category', 'b': 'category'})
        tm.assert_frame_equal(result, expected)

        expected = pd.DataFrame(columns=['a', 'b'], dtype='datetime64[ns]')
        result = self.read_csv(StringIO(data), header=0,
                               dtype='datetime64[ns]')
        tm.assert_frame_equal(result, expected)

        expected = pd.DataFrame({'a': pd.Series([], dtype='timedelta64[ns]'),
                                 'b': pd.Series([], dtype='timedelta64[ns]')},
                                index=[])
        result = self.read_csv(StringIO(data), header=0,
                               dtype='timedelta64[ns]')
        tm.assert_frame_equal(result, expected)

        expected = pd.DataFrame(columns=['a', 'b'])
        expected['a'] = expected['a'].astype(np.float64)
        result = self.read_csv(StringIO(data), header=0,
                               dtype={'a': np.float64})
        tm.assert_frame_equal(result, expected)

        expected = pd.DataFrame(columns=['a', 'b'])
        expected['a'] = expected['a'].astype(np.float64)
        result = self.read_csv(StringIO(data), header=0,
                               dtype={0: np.float64})
        tm.assert_frame_equal(result, expected)

        expected = pd.DataFrame(columns=['a', 'b'])
        expected['a'] = expected['a'].astype(np.int32)
        expected['b'] = expected['b'].astype(np.float64)
        result = self.read_csv(StringIO(data), header=0,
                               dtype={'a': np.int32, 1: np.float64})
        tm.assert_frame_equal(result, expected)

    def test_numeric_dtype(self):
        data = '0\n1'

        for dt in np.typecodes['AllInteger'] + np.typecodes['Float']:
            expected = pd.DataFrame([0, 1], dtype=dt)
            result = self.read_csv(StringIO(data), header=None, dtype=dt)
            tm.assert_frame_equal(expected, result)
