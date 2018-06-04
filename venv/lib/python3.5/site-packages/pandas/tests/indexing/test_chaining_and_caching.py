from warnings import catch_warnings

import pytest

import numpy as np
import pandas as pd
from pandas.core import common as com
from pandas import (compat, DataFrame, option_context,
                    Series, MultiIndex, date_range, Timestamp)
from pandas.util import testing as tm


class TestCaching(object):

    def test_slice_consolidate_invalidate_item_cache(self):

        # this is chained assignment, but will 'work'
        with option_context('chained_assignment', None):

            # #3970
            df = DataFrame({"aa": compat.lrange(5), "bb": [2.2] * 5})

            # Creates a second float block
            df["cc"] = 0.0

            # caches a reference to the 'bb' series
            df["bb"]

            # repr machinery triggers consolidation
            repr(df)

            # Assignment to wrong series
            df['bb'].iloc[0] = 0.17
            df._clear_item_cache()
            tm.assert_almost_equal(df['bb'][0], 0.17)

    def test_setitem_cache_updating(self):
        # GH 5424
        cont = ['one', 'two', 'three', 'four', 'five', 'six', 'seven']

        for do_ref in [False, False]:
            df = DataFrame({'a': cont,
                            "b": cont[3:] + cont[:3],
                            'c': np.arange(7)})

            # ref the cache
            if do_ref:
                df.loc[0, "c"]

            # set it
            df.loc[7, 'c'] = 1

            assert df.loc[0, 'c'] == 0.0
            assert df.loc[7, 'c'] == 1.0

        # GH 7084
        # not updating cache on series setting with slices
        expected = DataFrame({'A': [600, 600, 600]},
                             index=date_range('5/7/2014', '5/9/2014'))
        out = DataFrame({'A': [0, 0, 0]},
                        index=date_range('5/7/2014', '5/9/2014'))
        df = DataFrame({'C': ['A', 'A', 'A'], 'D': [100, 200, 300]})

        # loop through df to update out
        six = Timestamp('5/7/2014')
        eix = Timestamp('5/9/2014')
        for ix, row in df.iterrows():
            out.loc[six:eix, row['C']] = out.loc[six:eix, row['C']] + row['D']

        tm.assert_frame_equal(out, expected)
        tm.assert_series_equal(out['A'], expected['A'])

        # try via a chain indexing
        # this actually works
        out = DataFrame({'A': [0, 0, 0]},
                        index=date_range('5/7/2014', '5/9/2014'))
        for ix, row in df.iterrows():
            v = out[row['C']][six:eix] + row['D']
            out[row['C']][six:eix] = v

        tm.assert_frame_equal(out, expected)
        tm.assert_series_equal(out['A'], expected['A'])

        out = DataFrame({'A': [0, 0, 0]},
                        index=date_range('5/7/2014', '5/9/2014'))
        for ix, row in df.iterrows():
            out.loc[six:eix, row['C']] += row['D']

        tm.assert_frame_equal(out, expected)
        tm.assert_series_equal(out['A'], expected['A'])


class TestChaining(object):

    def test_setitem_chained_setfault(self):

        # GH6026
        # setfaults under numpy 1.7.1 (ok on 1.8)
        data = ['right', 'left', 'left', 'left', 'right', 'left', 'timeout']
        mdata = ['right', 'left', 'left', 'left', 'right', 'left', 'none']

        df = DataFrame({'response': np.array(data)})
        mask = df.response == 'timeout'
        df.response[mask] = 'none'
        tm.assert_frame_equal(df, DataFrame({'response': mdata}))

        recarray = np.rec.fromarrays([data], names=['response'])
        df = DataFrame(recarray)
        mask = df.response == 'timeout'
        df.response[mask] = 'none'
        tm.assert_frame_equal(df, DataFrame({'response': mdata}))

        df = DataFrame({'response': data, 'response1': data})
        mask = df.response == 'timeout'
        df.response[mask] = 'none'
        tm.assert_frame_equal(df, DataFrame({'response': mdata,
                                             'response1': data}))

        # GH 6056
        expected = DataFrame(dict(A=[np.nan, 'bar', 'bah', 'foo', 'bar']))
        df = DataFrame(dict(A=np.array(['foo', 'bar', 'bah', 'foo', 'bar'])))
        df['A'].iloc[0] = np.nan
        result = df.head()
        tm.assert_frame_equal(result, expected)

        df = DataFrame(dict(A=np.array(['foo', 'bar', 'bah', 'foo', 'bar'])))
        df.A.iloc[0] = np.nan
        result = df.head()
        tm.assert_frame_equal(result, expected)

    def test_detect_chained_assignment(self):

        pd.set_option('chained_assignment', 'raise')

        # work with the chain
        expected = DataFrame([[-5, 1], [-6, 3]], columns=list('AB'))
        df = DataFrame(np.arange(4).reshape(2, 2),
                       columns=list('AB'), dtype='int64')
        assert df._is_copy is None

        df['A'][0] = -5
        df['A'][1] = -6
        tm.assert_frame_equal(df, expected)

        # test with the chaining
        df = DataFrame({'A': Series(range(2), dtype='int64'),
                        'B': np.array(np.arange(2, 4), dtype=np.float64)})
        assert df._is_copy is None

        with pytest.raises(com.SettingWithCopyError):
            df['A'][0] = -5

        with pytest.raises(com.SettingWithCopyError):
            df['A'][1] = np.nan

        assert df['A']._is_copy is None

        # Using a copy (the chain), fails
        df = DataFrame({'A': Series(range(2), dtype='int64'),
                        'B': np.array(np.arange(2, 4), dtype=np.float64)})

        with pytest.raises(com.SettingWithCopyError):
            df.loc[0]['A'] = -5

        # Doc example
        df = DataFrame({'a': ['one', 'one', 'two', 'three',
                              'two', 'one', 'six'],
                        'c': Series(range(7), dtype='int64')})
        assert df._is_copy is None

        with pytest.raises(com.SettingWithCopyError):
            indexer = df.a.str.startswith('o')
            df[indexer]['c'] = 42

        expected = DataFrame({'A': [111, 'bbb', 'ccc'], 'B': [1, 2, 3]})
        df = DataFrame({'A': ['aaa', 'bbb', 'ccc'], 'B': [1, 2, 3]})

        with pytest.raises(com.SettingWithCopyError):
            df['A'][0] = 111

        with pytest.raises(com.SettingWithCopyError):
            df.loc[0]['A'] = 111

        df.loc[0, 'A'] = 111
        tm.assert_frame_equal(df, expected)

        # gh-5475: Make sure that is_copy is picked up reconstruction
        df = DataFrame({"A": [1, 2]})
        assert df._is_copy is None

        with tm.ensure_clean('__tmp__pickle') as path:
            df.to_pickle(path)
            df2 = pd.read_pickle(path)
            df2["B"] = df2["A"]
            df2["B"] = df2["A"]

        # gh-5597: a spurious raise as we are setting the entire column here
        from string import ascii_letters as letters

        def random_text(nobs=100):
            df = []
            for i in range(nobs):
                idx = np.random.randint(len(letters), size=2)
                idx.sort()

                df.append([letters[idx[0]:idx[1]]])

            return DataFrame(df, columns=['letters'])

        df = random_text(100000)

        # Always a copy
        x = df.iloc[[0, 1, 2]]
        assert x._is_copy is not None

        x = df.iloc[[0, 1, 2, 4]]
        assert x._is_copy is not None

        # Explicitly copy
        indexer = df.letters.apply(lambda x: len(x) > 10)
        df = df.loc[indexer].copy()

        assert df._is_copy is None
        df['letters'] = df['letters'].apply(str.lower)

        # Implicitly take
        df = random_text(100000)
        indexer = df.letters.apply(lambda x: len(x) > 10)
        df = df.loc[indexer]

        assert df._is_copy is not None
        df['letters'] = df['letters'].apply(str.lower)

        # Implicitly take 2
        df = random_text(100000)
        indexer = df.letters.apply(lambda x: len(x) > 10)

        df = df.loc[indexer]
        assert df._is_copy is not None
        df.loc[:, 'letters'] = df['letters'].apply(str.lower)

        # Should be ok even though it's a copy!
        assert df._is_copy is None

        df['letters'] = df['letters'].apply(str.lower)
        assert df._is_copy is None

        df = random_text(100000)
        indexer = df.letters.apply(lambda x: len(x) > 10)
        df.loc[indexer, 'letters'] = (
            df.loc[indexer, 'letters'].apply(str.lower))

        # an identical take, so no copy
        df = DataFrame({'a': [1]}).dropna()
        assert df._is_copy is None
        df['a'] += 1

        # Inplace ops, originally from:
        # http://stackoverflow.com/questions/20508968/series-fillna-in-a-multiindex-dataframe-does-not-fill-is-this-a-bug
        a = [12, 23]
        b = [123, None]
        c = [1234, 2345]
        d = [12345, 23456]
        tuples = [('eyes', 'left'), ('eyes', 'right'), ('ears', 'left'),
                  ('ears', 'right')]
        events = {('eyes', 'left'): a,
                  ('eyes', 'right'): b,
                  ('ears', 'left'): c,
                  ('ears', 'right'): d}
        multiind = MultiIndex.from_tuples(tuples, names=['part', 'side'])
        zed = DataFrame(events, index=['a', 'b'], columns=multiind)

        with pytest.raises(com.SettingWithCopyError):
            zed['eyes']['right'].fillna(value=555, inplace=True)

        df = DataFrame(np.random.randn(10, 4))
        s = df.iloc[:, 0].sort_values()

        tm.assert_series_equal(s, df.iloc[:, 0].sort_values())
        tm.assert_series_equal(s, df[0].sort_values())

        # see gh-6025: false positives
        df = DataFrame({'column1': ['a', 'a', 'a'], 'column2': [4, 8, 9]})
        str(df)

        df['column1'] = df['column1'] + 'b'
        str(df)

        df = df[df['column2'] != 8]
        str(df)

        df['column1'] = df['column1'] + 'c'
        str(df)

        # from SO:
        # http://stackoverflow.com/questions/24054495/potential-bug-setting-value-for-undefined-column-using-iloc
        df = DataFrame(np.arange(0, 9), columns=['count'])
        df['group'] = 'b'

        with pytest.raises(com.SettingWithCopyError):
            df.iloc[0:5]['group'] = 'a'

        # Mixed type setting but same dtype & changing dtype
        df = DataFrame(dict(A=date_range('20130101', periods=5),
                            B=np.random.randn(5),
                            C=np.arange(5, dtype='int64'),
                            D=list('abcde')))

        with pytest.raises(com.SettingWithCopyError):
            df.loc[2]['D'] = 'foo'

        with pytest.raises(com.SettingWithCopyError):
            df.loc[2]['C'] = 'foo'

        with pytest.raises(com.SettingWithCopyError):
            df['C'][2] = 'foo'

    def test_setting_with_copy_bug(self):

        # operating on a copy
        df = DataFrame({'a': list(range(4)),
                        'b': list('ab..'),
                        'c': ['a', 'b', np.nan, 'd']})
        mask = pd.isna(df.c)

        def f():
            df[['c']][mask] = df[['b']][mask]

        pytest.raises(com.SettingWithCopyError, f)

        # invalid warning as we are returning a new object
        # GH 8730
        df1 = DataFrame({'x': Series(['a', 'b', 'c']),
                         'y': Series(['d', 'e', 'f'])})
        df2 = df1[['x']]

        # this should not raise
        df2['y'] = ['g', 'h', 'i']

    def test_detect_chained_assignment_warnings(self):

        # warnings
        with option_context('chained_assignment', 'warn'):
            df = DataFrame({'A': ['aaa', 'bbb', 'ccc'], 'B': [1, 2, 3]})
            with tm.assert_produces_warning(
                    expected_warning=com.SettingWithCopyWarning):
                df.loc[0]['A'] = 111

    def test_chained_getitem_with_lists(self):

        # GH6394
        # Regression in chained getitem indexing with embedded list-like from
        # 0.12
        def check(result, expected):
            tm.assert_numpy_array_equal(result, expected)
            assert isinstance(result, np.ndarray)

        df = DataFrame({'A': 5 * [np.zeros(3)], 'B': 5 * [np.ones(3)]})
        expected = df['A'].iloc[2]
        result = df.loc[2, 'A']
        check(result, expected)
        result2 = df.iloc[2]['A']
        check(result2, expected)
        result3 = df['A'].loc[2]
        check(result3, expected)
        result4 = df['A'].iloc[2]
        check(result4, expected)

    def test_cache_updating(self):
        # GH 4939, make sure to update the cache on setitem

        df = tm.makeDataFrame()
        df['A']  # cache series
        with catch_warnings(record=True):
            df.ix["Hello Friend"] = df.ix[0]
        assert "Hello Friend" in df['A'].index
        assert "Hello Friend" in df['B'].index

        with catch_warnings(record=True):
            panel = tm.makePanel()
            panel.ix[0]  # get first item into cache
            panel.ix[:, :, 'A+1'] = panel.ix[:, :, 'A'] + 1
            assert "A+1" in panel.ix[0].columns
            assert "A+1" in panel.ix[1].columns

        # 5216
        # make sure that we don't try to set a dead cache
        a = np.random.rand(10, 3)
        df = DataFrame(a, columns=['x', 'y', 'z'])
        tuples = [(i, j) for i in range(5) for j in range(2)]
        index = MultiIndex.from_tuples(tuples)
        df.index = index

        # setting via chained assignment
        # but actually works, since everything is a view
        df.loc[0]['z'].iloc[0] = 1.
        result = df.loc[(0, 0), 'z']
        assert result == 1

        # correct setting
        df.loc[(0, 0), 'z'] = 2
        result = df.loc[(0, 0), 'z']
        assert result == 2

        # 10264
        df = DataFrame(np.zeros((5, 5), dtype='int64'), columns=[
                       'a', 'b', 'c', 'd', 'e'], index=range(5))
        df['f'] = 0
        df.f.values[3] = 1

        # TODO(wesm): unused?
        # y = df.iloc[np.arange(2, len(df))]

        df.f.values[3] = 2
        expected = DataFrame(np.zeros((5, 6), dtype='int64'), columns=[
                             'a', 'b', 'c', 'd', 'e', 'f'], index=range(5))
        expected.at[3, 'f'] = 2
        tm.assert_frame_equal(df, expected)
        expected = Series([0, 0, 0, 2, 0], name='f')
        tm.assert_series_equal(df.f, expected)

    def test_deprecate_is_copy(self):
        # GH18801
        df = DataFrame({"A": [1, 2, 3]})
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            # getter
            df.is_copy

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            # setter
            df.is_copy = "test deprecated is_copy"
