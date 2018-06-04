# -*- coding: utf-8 -*-
# pylint: disable-msg=W0612,E1101

import numpy as np
import pandas as pd
import pandas.util.testing as tm


class TestIndexingCallable(object):

    def test_frame_loc_ix_callable(self):
        # GH 11485
        df = pd.DataFrame({'A': [1, 2, 3, 4], 'B': list('aabb'),
                           'C': [1, 2, 3, 4]})
        # iloc cannot use boolean Series (see GH3635)

        # return bool indexer
        res = df.loc[lambda x: x.A > 2]
        tm.assert_frame_equal(res, df.loc[df.A > 2])

        res = df.loc[lambda x: x.A > 2]
        tm.assert_frame_equal(res, df.loc[df.A > 2])

        res = df.loc[lambda x: x.A > 2, ]
        tm.assert_frame_equal(res, df.loc[df.A > 2, ])

        res = df.loc[lambda x: x.A > 2, ]
        tm.assert_frame_equal(res, df.loc[df.A > 2, ])

        res = df.loc[lambda x: x.B == 'b', :]
        tm.assert_frame_equal(res, df.loc[df.B == 'b', :])

        res = df.loc[lambda x: x.B == 'b', :]
        tm.assert_frame_equal(res, df.loc[df.B == 'b', :])

        res = df.loc[lambda x: x.A > 2, lambda x: x.columns == 'B']
        tm.assert_frame_equal(res, df.loc[df.A > 2, [False, True, False]])

        res = df.loc[lambda x: x.A > 2, lambda x: x.columns == 'B']
        tm.assert_frame_equal(res, df.loc[df.A > 2, [False, True, False]])

        res = df.loc[lambda x: x.A > 2, lambda x: 'B']
        tm.assert_series_equal(res, df.loc[df.A > 2, 'B'])

        res = df.loc[lambda x: x.A > 2, lambda x: 'B']
        tm.assert_series_equal(res, df.loc[df.A > 2, 'B'])

        res = df.loc[lambda x: x.A > 2, lambda x: ['A', 'B']]
        tm.assert_frame_equal(res, df.loc[df.A > 2, ['A', 'B']])

        res = df.loc[lambda x: x.A > 2, lambda x: ['A', 'B']]
        tm.assert_frame_equal(res, df.loc[df.A > 2, ['A', 'B']])

        res = df.loc[lambda x: x.A == 2, lambda x: ['A', 'B']]
        tm.assert_frame_equal(res, df.loc[df.A == 2, ['A', 'B']])

        res = df.loc[lambda x: x.A == 2, lambda x: ['A', 'B']]
        tm.assert_frame_equal(res, df.loc[df.A == 2, ['A', 'B']])

        # scalar
        res = df.loc[lambda x: 1, lambda x: 'A']
        assert res == df.loc[1, 'A']

        res = df.loc[lambda x: 1, lambda x: 'A']
        assert res == df.loc[1, 'A']

    def test_frame_loc_ix_callable_mixture(self):
        # GH 11485
        df = pd.DataFrame({'A': [1, 2, 3, 4], 'B': list('aabb'),
                           'C': [1, 2, 3, 4]})

        res = df.loc[lambda x: x.A > 2, ['A', 'B']]
        tm.assert_frame_equal(res, df.loc[df.A > 2, ['A', 'B']])

        res = df.loc[lambda x: x.A > 2, ['A', 'B']]
        tm.assert_frame_equal(res, df.loc[df.A > 2, ['A', 'B']])

        res = df.loc[[2, 3], lambda x: ['A', 'B']]
        tm.assert_frame_equal(res, df.loc[[2, 3], ['A', 'B']])

        res = df.loc[[2, 3], lambda x: ['A', 'B']]
        tm.assert_frame_equal(res, df.loc[[2, 3], ['A', 'B']])

        res = df.loc[3, lambda x: ['A', 'B']]
        tm.assert_series_equal(res, df.loc[3, ['A', 'B']])

        res = df.loc[3, lambda x: ['A', 'B']]
        tm.assert_series_equal(res, df.loc[3, ['A', 'B']])

    def test_frame_loc_callable(self):
        # GH 11485
        df = pd.DataFrame({'X': [1, 2, 3, 4],
                           'Y': list('aabb')},
                          index=list('ABCD'))

        # return label
        res = df.loc[lambda x: ['A', 'C']]
        tm.assert_frame_equal(res, df.loc[['A', 'C']])

        res = df.loc[lambda x: ['A', 'C'], ]
        tm.assert_frame_equal(res, df.loc[['A', 'C'], ])

        res = df.loc[lambda x: ['A', 'C'], :]
        tm.assert_frame_equal(res, df.loc[['A', 'C'], :])

        res = df.loc[lambda x: ['A', 'C'], lambda x: 'X']
        tm.assert_series_equal(res, df.loc[['A', 'C'], 'X'])

        res = df.loc[lambda x: ['A', 'C'], lambda x: ['X']]
        tm.assert_frame_equal(res, df.loc[['A', 'C'], ['X']])

        # mixture
        res = df.loc[['A', 'C'], lambda x: 'X']
        tm.assert_series_equal(res, df.loc[['A', 'C'], 'X'])

        res = df.loc[['A', 'C'], lambda x: ['X']]
        tm.assert_frame_equal(res, df.loc[['A', 'C'], ['X']])

        res = df.loc[lambda x: ['A', 'C'], 'X']
        tm.assert_series_equal(res, df.loc[['A', 'C'], 'X'])

        res = df.loc[lambda x: ['A', 'C'], ['X']]
        tm.assert_frame_equal(res, df.loc[['A', 'C'], ['X']])

    def test_frame_loc_callable_setitem(self):
        # GH 11485
        df = pd.DataFrame({'X': [1, 2, 3, 4],
                           'Y': list('aabb')},
                          index=list('ABCD'))

        # return label
        res = df.copy()
        res.loc[lambda x: ['A', 'C']] = -20
        exp = df.copy()
        exp.loc[['A', 'C']] = -20
        tm.assert_frame_equal(res, exp)

        res = df.copy()
        res.loc[lambda x: ['A', 'C'], :] = 20
        exp = df.copy()
        exp.loc[['A', 'C'], :] = 20
        tm.assert_frame_equal(res, exp)

        res = df.copy()
        res.loc[lambda x: ['A', 'C'], lambda x: 'X'] = -1
        exp = df.copy()
        exp.loc[['A', 'C'], 'X'] = -1
        tm.assert_frame_equal(res, exp)

        res = df.copy()
        res.loc[lambda x: ['A', 'C'], lambda x: ['X']] = [5, 10]
        exp = df.copy()
        exp.loc[['A', 'C'], ['X']] = [5, 10]
        tm.assert_frame_equal(res, exp)

        # mixture
        res = df.copy()
        res.loc[['A', 'C'], lambda x: 'X'] = np.array([-1, -2])
        exp = df.copy()
        exp.loc[['A', 'C'], 'X'] = np.array([-1, -2])
        tm.assert_frame_equal(res, exp)

        res = df.copy()
        res.loc[['A', 'C'], lambda x: ['X']] = 10
        exp = df.copy()
        exp.loc[['A', 'C'], ['X']] = 10
        tm.assert_frame_equal(res, exp)

        res = df.copy()
        res.loc[lambda x: ['A', 'C'], 'X'] = -2
        exp = df.copy()
        exp.loc[['A', 'C'], 'X'] = -2
        tm.assert_frame_equal(res, exp)

        res = df.copy()
        res.loc[lambda x: ['A', 'C'], ['X']] = -4
        exp = df.copy()
        exp.loc[['A', 'C'], ['X']] = -4
        tm.assert_frame_equal(res, exp)

    def test_frame_iloc_callable(self):
        # GH 11485
        df = pd.DataFrame({'X': [1, 2, 3, 4],
                           'Y': list('aabb')},
                          index=list('ABCD'))

        # return location
        res = df.iloc[lambda x: [1, 3]]
        tm.assert_frame_equal(res, df.iloc[[1, 3]])

        res = df.iloc[lambda x: [1, 3], :]
        tm.assert_frame_equal(res, df.iloc[[1, 3], :])

        res = df.iloc[lambda x: [1, 3], lambda x: 0]
        tm.assert_series_equal(res, df.iloc[[1, 3], 0])

        res = df.iloc[lambda x: [1, 3], lambda x: [0]]
        tm.assert_frame_equal(res, df.iloc[[1, 3], [0]])

        # mixture
        res = df.iloc[[1, 3], lambda x: 0]
        tm.assert_series_equal(res, df.iloc[[1, 3], 0])

        res = df.iloc[[1, 3], lambda x: [0]]
        tm.assert_frame_equal(res, df.iloc[[1, 3], [0]])

        res = df.iloc[lambda x: [1, 3], 0]
        tm.assert_series_equal(res, df.iloc[[1, 3], 0])

        res = df.iloc[lambda x: [1, 3], [0]]
        tm.assert_frame_equal(res, df.iloc[[1, 3], [0]])

    def test_frame_iloc_callable_setitem(self):
        # GH 11485
        df = pd.DataFrame({'X': [1, 2, 3, 4],
                           'Y': list('aabb')},
                          index=list('ABCD'))

        # return location
        res = df.copy()
        res.iloc[lambda x: [1, 3]] = 0
        exp = df.copy()
        exp.iloc[[1, 3]] = 0
        tm.assert_frame_equal(res, exp)

        res = df.copy()
        res.iloc[lambda x: [1, 3], :] = -1
        exp = df.copy()
        exp.iloc[[1, 3], :] = -1
        tm.assert_frame_equal(res, exp)

        res = df.copy()
        res.iloc[lambda x: [1, 3], lambda x: 0] = 5
        exp = df.copy()
        exp.iloc[[1, 3], 0] = 5
        tm.assert_frame_equal(res, exp)

        res = df.copy()
        res.iloc[lambda x: [1, 3], lambda x: [0]] = 25
        exp = df.copy()
        exp.iloc[[1, 3], [0]] = 25
        tm.assert_frame_equal(res, exp)

        # mixture
        res = df.copy()
        res.iloc[[1, 3], lambda x: 0] = -3
        exp = df.copy()
        exp.iloc[[1, 3], 0] = -3
        tm.assert_frame_equal(res, exp)

        res = df.copy()
        res.iloc[[1, 3], lambda x: [0]] = -5
        exp = df.copy()
        exp.iloc[[1, 3], [0]] = -5
        tm.assert_frame_equal(res, exp)

        res = df.copy()
        res.iloc[lambda x: [1, 3], 0] = 10
        exp = df.copy()
        exp.iloc[[1, 3], 0] = 10
        tm.assert_frame_equal(res, exp)

        res = df.copy()
        res.iloc[lambda x: [1, 3], [0]] = [-5, -5]
        exp = df.copy()
        exp.iloc[[1, 3], [0]] = [-5, -5]
        tm.assert_frame_equal(res, exp)
