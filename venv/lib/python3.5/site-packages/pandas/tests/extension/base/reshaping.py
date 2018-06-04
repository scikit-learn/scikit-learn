import pytest
import numpy as np

import pandas as pd
from pandas.core.internals import ExtensionBlock

from .base import BaseExtensionTests


class BaseReshapingTests(BaseExtensionTests):
    """Tests for reshaping and concatenation."""
    @pytest.mark.parametrize('in_frame', [True, False])
    def test_concat(self, data, in_frame):
        wrapped = pd.Series(data)
        if in_frame:
            wrapped = pd.DataFrame(wrapped)
        result = pd.concat([wrapped, wrapped], ignore_index=True)

        assert len(result) == len(data) * 2

        if in_frame:
            dtype = result.dtypes[0]
        else:
            dtype = result.dtype

        assert dtype == data.dtype
        assert isinstance(result._data.blocks[0], ExtensionBlock)

    @pytest.mark.parametrize('in_frame', [True, False])
    def test_concat_all_na_block(self, data_missing, in_frame):
        valid_block = pd.Series(data_missing.take([1, 1]), index=[0, 1])
        na_block = pd.Series(data_missing.take([0, 0]), index=[2, 3])
        if in_frame:
            valid_block = pd.DataFrame({"a": valid_block})
            na_block = pd.DataFrame({"a": na_block})
        result = pd.concat([valid_block, na_block])
        if in_frame:
            expected = pd.DataFrame({"a": data_missing.take([1, 1, 0, 0])})
            self.assert_frame_equal(result, expected)
        else:
            expected = pd.Series(data_missing.take([1, 1, 0, 0]))
            self.assert_series_equal(result, expected)

    def test_concat_mixed_dtypes(self, data):
        # https://github.com/pandas-dev/pandas/issues/20762
        df1 = pd.DataFrame({'A': data[:3]})
        df2 = pd.DataFrame({"A": [1, 2, 3]})
        df3 = pd.DataFrame({"A": ['a', 'b', 'c']}).astype('category')
        df4 = pd.DataFrame({"A": pd.SparseArray([1, 2, 3])})
        dfs = [df1, df2, df3, df4]

        # dataframes
        result = pd.concat(dfs)
        expected = pd.concat([x.astype(object) for x in dfs])
        self.assert_frame_equal(result, expected)

        # series
        result = pd.concat([x['A'] for x in dfs])
        expected = pd.concat([x['A'].astype(object) for x in dfs])
        self.assert_series_equal(result, expected)

        # simple test for just EA and one other
        result = pd.concat([df1, df2])
        expected = pd.concat([df1.astype('object'), df2.astype('object')])
        self.assert_frame_equal(result, expected)

        result = pd.concat([df1['A'], df2['A']])
        expected = pd.concat([df1['A'].astype('object'),
                              df2['A'].astype('object')])
        self.assert_series_equal(result, expected)

    def test_concat_columns(self, data, na_value):
        df1 = pd.DataFrame({'A': data[:3]})
        df2 = pd.DataFrame({'B': [1, 2, 3]})

        expected = pd.DataFrame({'A': data[:3], 'B': [1, 2, 3]})
        result = pd.concat([df1, df2], axis=1)
        self.assert_frame_equal(result, expected)
        result = pd.concat([df1['A'], df2['B']], axis=1)
        self.assert_frame_equal(result, expected)

        # non-aligned
        df2 = pd.DataFrame({'B': [1, 2, 3]}, index=[1, 2, 3])
        expected = pd.DataFrame({
            'A': data._from_sequence(list(data[:3]) + [na_value]),
            'B': [np.nan, 1, 2, 3]})
        result = pd.concat([df1, df2], axis=1)
        self.assert_frame_equal(result, expected)
        result = pd.concat([df1['A'], df2['B']], axis=1)
        self.assert_frame_equal(result, expected)

    def test_align(self, data, na_value):
        a = data[:3]
        b = data[2:5]
        r1, r2 = pd.Series(a).align(pd.Series(b, index=[1, 2, 3]))

        # Assumes that the ctor can take a list of scalars of the type
        e1 = pd.Series(data._from_sequence(list(a) + [na_value]))
        e2 = pd.Series(data._from_sequence([na_value] + list(b)))
        self.assert_series_equal(r1, e1)
        self.assert_series_equal(r2, e2)

    def test_align_frame(self, data, na_value):
        a = data[:3]
        b = data[2:5]
        r1, r2 = pd.DataFrame({'A': a}).align(
            pd.DataFrame({'A': b}, index=[1, 2, 3])
        )

        # Assumes that the ctor can take a list of scalars of the type
        e1 = pd.DataFrame({'A': data._from_sequence(list(a) + [na_value])})
        e2 = pd.DataFrame({'A': data._from_sequence([na_value] + list(b))})
        self.assert_frame_equal(r1, e1)
        self.assert_frame_equal(r2, e2)

    def test_align_series_frame(self, data, na_value):
        # https://github.com/pandas-dev/pandas/issues/20576
        ser = pd.Series(data, name='a')
        df = pd.DataFrame({"col": np.arange(len(ser) + 1)})
        r1, r2 = ser.align(df)

        e1 = pd.Series(data._from_sequence(list(data) + [na_value]),
                       name=ser.name)

        self.assert_series_equal(r1, e1)
        self.assert_frame_equal(r2, df)

    def test_set_frame_expand_regular_with_extension(self, data):
        df = pd.DataFrame({"A": [1] * len(data)})
        df['B'] = data
        expected = pd.DataFrame({"A": [1] * len(data), "B": data})
        self.assert_frame_equal(df, expected)

    def test_set_frame_expand_extension_with_regular(self, data):
        df = pd.DataFrame({'A': data})
        df['B'] = [1] * len(data)
        expected = pd.DataFrame({"A": data, "B": [1] * len(data)})
        self.assert_frame_equal(df, expected)

    def test_set_frame_overwrite_object(self, data):
        # https://github.com/pandas-dev/pandas/issues/20555
        df = pd.DataFrame({"A": [1] * len(data)}, dtype=object)
        df['A'] = data
        assert df.dtypes['A'] == data.dtype

    def test_merge(self, data, na_value):
        # GH-20743
        df1 = pd.DataFrame({'ext': data[:3], 'int1': [1, 2, 3],
                            'key': [0, 1, 2]})
        df2 = pd.DataFrame({'int2': [1, 2, 3, 4], 'key': [0, 0, 1, 3]})

        res = pd.merge(df1, df2)
        exp = pd.DataFrame(
            {'int1': [1, 1, 2], 'int2': [1, 2, 3], 'key': [0, 0, 1],
             'ext': data._from_sequence([data[0], data[0], data[1]])})
        self.assert_frame_equal(res, exp[['ext', 'int1', 'key', 'int2']])

        res = pd.merge(df1, df2, how='outer')
        exp = pd.DataFrame(
            {'int1': [1, 1, 2, 3, np.nan], 'int2': [1, 2, 3, np.nan, 4],
             'key': [0, 0, 1, 2, 3],
             'ext': data._from_sequence(
                 [data[0], data[0], data[1], data[2], na_value])})
        self.assert_frame_equal(res, exp[['ext', 'int1', 'key', 'int2']])
