# -*- coding: utf-8 -*-

import warnings

import numpy as np
import pandas as pd
from pandas.core.api import Series, DataFrame, MultiIndex
import pandas.util.testing as tm
import pytest


class TestIndexingSlow(object):

    @pytest.mark.slow
    def test_multiindex_get_loc(self):  # GH7724, GH2646

        with warnings.catch_warnings(record=True):

            # test indexing into a multi-index before & past the lexsort depth
            from numpy.random import randint, choice, randn
            cols = ['jim', 'joe', 'jolie', 'joline', 'jolia']

            def validate(mi, df, key):
                mask = np.ones(len(df)).astype('bool')

                # test for all partials of this key
                for i, k in enumerate(key):
                    mask &= df.iloc[:, i] == k

                    if not mask.any():
                        assert key[:i + 1] not in mi.index
                        continue

                    assert key[:i + 1] in mi.index
                    right = df[mask].copy()

                    if i + 1 != len(key):  # partial key
                        right.drop(cols[:i + 1], axis=1, inplace=True)
                        right.set_index(cols[i + 1:-1], inplace=True)
                        tm.assert_frame_equal(mi.loc[key[:i + 1]], right)

                    else:  # full key
                        right.set_index(cols[:-1], inplace=True)
                        if len(right) == 1:  # single hit
                            right = Series(right['jolia'].values,
                                           name=right.index[0],
                                           index=['jolia'])
                            tm.assert_series_equal(mi.loc[key[:i + 1]], right)
                        else:  # multi hit
                            tm.assert_frame_equal(mi.loc[key[:i + 1]], right)

            def loop(mi, df, keys):
                for key in keys:
                    validate(mi, df, key)

            n, m = 1000, 50

            vals = [randint(0, 10, n), choice(
                list('abcdefghij'), n), choice(
                    pd.date_range('20141009', periods=10).tolist(), n), choice(
                        list('ZYXWVUTSRQ'), n), randn(n)]
            vals = list(map(tuple, zip(*vals)))

            # bunch of keys for testing
            keys = [randint(0, 11, m), choice(
                list('abcdefghijk'), m), choice(
                    pd.date_range('20141009', periods=11).tolist(), m), choice(
                        list('ZYXWVUTSRQP'), m)]
            keys = list(map(tuple, zip(*keys)))
            keys += list(map(lambda t: t[:-1], vals[::n // m]))

            # covers both unique index and non-unique index
            df = DataFrame(vals, columns=cols)
            a, b = pd.concat([df, df]), df.drop_duplicates(subset=cols[:-1])

            for frame in a, b:
                for i in range(5):  # lexsort depth
                    df = frame.copy() if i == 0 else frame.sort_values(
                        by=cols[:i])
                    mi = df.set_index(cols[:-1])
                    assert not mi.index.lexsort_depth < i
                    loop(mi, df, keys)

    @pytest.mark.slow
    def test_large_dataframe_indexing(self):
        # GH10692
        result = DataFrame({'x': range(10 ** 6)}, dtype='int64')
        result.loc[len(result)] = len(result) + 1
        expected = DataFrame({'x': range(10 ** 6 + 1)}, dtype='int64')
        tm.assert_frame_equal(result, expected)

    @pytest.mark.slow
    def test_large_mi_dataframe_indexing(self):
        # GH10645
        result = MultiIndex.from_arrays([range(10 ** 6), range(10 ** 6)])
        assert (not (10 ** 6, 0) in result)
