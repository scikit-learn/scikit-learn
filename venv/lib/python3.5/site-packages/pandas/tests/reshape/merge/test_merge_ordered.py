import pandas as pd
from pandas import DataFrame, merge_ordered
from pandas.util import testing as tm
from pandas.util.testing import assert_frame_equal

from numpy import nan


class TestMergeOrdered(object):

    def setup_method(self, method):
        self.left = DataFrame({'key': ['a', 'c', 'e'],
                               'lvalue': [1, 2., 3]})

        self.right = DataFrame({'key': ['b', 'c', 'd', 'f'],
                                'rvalue': [1, 2, 3., 4]})

    def test_basic(self):
        result = merge_ordered(self.left, self.right, on='key')
        expected = DataFrame({'key': ['a', 'b', 'c', 'd', 'e', 'f'],
                              'lvalue': [1, nan, 2, nan, 3, nan],
                              'rvalue': [nan, 1, 2, 3, nan, 4]})

        assert_frame_equal(result, expected)

    def test_ffill(self):
        result = merge_ordered(
            self.left, self.right, on='key', fill_method='ffill')
        expected = DataFrame({'key': ['a', 'b', 'c', 'd', 'e', 'f'],
                              'lvalue': [1., 1, 2, 2, 3, 3.],
                              'rvalue': [nan, 1, 2, 3, 3, 4]})
        assert_frame_equal(result, expected)

    def test_multigroup(self):
        left = pd.concat([self.left, self.left], ignore_index=True)

        left['group'] = ['a'] * 3 + ['b'] * 3

        result = merge_ordered(left, self.right, on='key', left_by='group',
                               fill_method='ffill')
        expected = DataFrame({'key': ['a', 'b', 'c', 'd', 'e', 'f'] * 2,
                              'lvalue': [1., 1, 2, 2, 3, 3.] * 2,
                              'rvalue': [nan, 1, 2, 3, 3, 4] * 2})
        expected['group'] = ['a'] * 6 + ['b'] * 6

        assert_frame_equal(result, expected.loc[:, result.columns])

        result2 = merge_ordered(self.right, left, on='key', right_by='group',
                                fill_method='ffill')
        assert_frame_equal(result, result2.loc[:, result.columns])

        result = merge_ordered(left, self.right, on='key', left_by='group')
        assert result['group'].notna().all()

    def test_merge_type(self):
        class NotADataFrame(DataFrame):

            @property
            def _constructor(self):
                return NotADataFrame

        nad = NotADataFrame(self.left)
        result = nad.merge(self.right, on='key')

        assert isinstance(result, NotADataFrame)

    def test_empty_sequence_concat(self):
        # GH 9157
        empty_pat = "[Nn]o objects"
        none_pat = "objects.*None"
        test_cases = [
            ((), empty_pat),
            ([], empty_pat),
            ({}, empty_pat),
            ([None], none_pat),
            ([None, None], none_pat)
        ]
        for df_seq, pattern in test_cases:
            tm.assert_raises_regex(ValueError, pattern, pd.concat, df_seq)

        pd.concat([pd.DataFrame()])
        pd.concat([None, pd.DataFrame()])
        pd.concat([pd.DataFrame(), None])

    def test_doc_example(self):
        left = DataFrame({'group': list('aaabbb'),
                          'key': ['a', 'c', 'e', 'a', 'c', 'e'],
                          'lvalue': [1, 2, 3] * 2,
                          })

        right = DataFrame({'key': ['b', 'c', 'd'],
                           'rvalue': [1, 2, 3]})

        result = merge_ordered(left, right, fill_method='ffill',
                               left_by='group')

        expected = DataFrame({'group': list('aaaaabbbbb'),
                              'key': ['a', 'b', 'c', 'd', 'e'] * 2,
                              'lvalue': [1, 1, 2, 2, 3] * 2,
                              'rvalue': [nan, 1, 2, 3, 3] * 2})

        assert_frame_equal(result, expected)
