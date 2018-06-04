import pytest
from warnings import catch_warnings

import numpy as np
from pandas.util import testing as tm
from pandas import Panel, date_range, DataFrame


class TestPanel(object):

    def test_iloc_getitem_panel(self):

        with catch_warnings(record=True):
            # GH 7189
            p = Panel(np.arange(4 * 3 * 2).reshape(4, 3, 2),
                      items=['A', 'B', 'C', 'D'],
                      major_axis=['a', 'b', 'c'],
                      minor_axis=['one', 'two'])

            result = p.iloc[1]
            expected = p.loc['B']
            tm.assert_frame_equal(result, expected)

            result = p.iloc[1, 1]
            expected = p.loc['B', 'b']
            tm.assert_series_equal(result, expected)

            result = p.iloc[1, 1, 1]
            expected = p.loc['B', 'b', 'two']
            assert result == expected

            # slice
            result = p.iloc[1:3]
            expected = p.loc[['B', 'C']]
            tm.assert_panel_equal(result, expected)

            result = p.iloc[:, 0:2]
            expected = p.loc[:, ['a', 'b']]
            tm.assert_panel_equal(result, expected)

            # list of integers
            result = p.iloc[[0, 2]]
            expected = p.loc[['A', 'C']]
            tm.assert_panel_equal(result, expected)

            # neg indicies
            result = p.iloc[[-1, 1], [-1, 1]]
            expected = p.loc[['D', 'B'], ['c', 'b']]
            tm.assert_panel_equal(result, expected)

            # dups indicies
            result = p.iloc[[-1, -1, 1], [-1, 1]]
            expected = p.loc[['D', 'D', 'B'], ['c', 'b']]
            tm.assert_panel_equal(result, expected)

            # combined
            result = p.iloc[0, [True, True], [0, 1]]
            expected = p.loc['A', ['a', 'b'], ['one', 'two']]
            tm.assert_frame_equal(result, expected)

            # out-of-bounds exception
            with pytest.raises(IndexError):
                p.iloc[tuple([10, 5])]

            def f():
                p.iloc[0, [True, True], [0, 1, 2]]

            pytest.raises(IndexError, f)

            # trying to use a label
            with pytest.raises(ValueError):
                p.iloc[tuple(['j', 'D'])]

            # GH
            p = Panel(
                np.random.rand(4, 3, 2), items=['A', 'B', 'C', 'D'],
                major_axis=['U', 'V', 'W'], minor_axis=['X', 'Y'])
            expected = p['A']

            result = p.iloc[0, :, :]
            tm.assert_frame_equal(result, expected)

            result = p.iloc[0, [True, True, True], :]
            tm.assert_frame_equal(result, expected)

            result = p.iloc[0, [True, True, True], [0, 1]]
            tm.assert_frame_equal(result, expected)

            def f():
                p.iloc[0, [True, True, True], [0, 1, 2]]

            pytest.raises(IndexError, f)

            def f():
                p.iloc[0, [True, True, True], [2]]

            pytest.raises(IndexError, f)

    def test_iloc_panel_issue(self):

        with catch_warnings(record=True):
            # see gh-3617
            p = Panel(np.random.randn(4, 4, 4))

            assert p.iloc[:3, :3, :3].shape == (3, 3, 3)
            assert p.iloc[1, :3, :3].shape == (3, 3)
            assert p.iloc[:3, 1, :3].shape == (3, 3)
            assert p.iloc[:3, :3, 1].shape == (3, 3)
            assert p.iloc[1, 1, :3].shape == (3, )
            assert p.iloc[1, :3, 1].shape == (3, )
            assert p.iloc[:3, 1, 1].shape == (3, )

    def test_panel_getitem(self):

        with catch_warnings(record=True):
            # GH4016, date selection returns a frame when a partial string
            # selection
            ind = date_range(start="2000", freq="D", periods=1000)
            df = DataFrame(
                np.random.randn(
                    len(ind), 5), index=ind, columns=list('ABCDE'))
            panel = Panel({'frame_' + c: df for c in list('ABC')})

            test2 = panel.loc[:, "2002":"2002-12-31"]
            test1 = panel.loc[:, "2002"]
            tm.assert_panel_equal(test1, test2)

            # GH8710
            # multi-element getting with a list
            panel = tm.makePanel()

            expected = panel.iloc[[0, 1]]

            result = panel.loc[['ItemA', 'ItemB']]
            tm.assert_panel_equal(result, expected)

            result = panel.loc[['ItemA', 'ItemB'], :, :]
            tm.assert_panel_equal(result, expected)

            result = panel[['ItemA', 'ItemB']]
            tm.assert_panel_equal(result, expected)

            result = panel.loc['ItemA':'ItemB']
            tm.assert_panel_equal(result, expected)

            with catch_warnings(record=True):
                result = panel.ix[['ItemA', 'ItemB']]
            tm.assert_panel_equal(result, expected)

            # with an object-like
            # GH 9140
            class TestObject(object):

                def __str__(self):
                    return "TestObject"

            obj = TestObject()

            p = Panel(np.random.randn(1, 5, 4), items=[obj],
                      major_axis=date_range('1/1/2000', periods=5),
                      minor_axis=['A', 'B', 'C', 'D'])

            expected = p.iloc[0]
            result = p[obj]
            tm.assert_frame_equal(result, expected)

    def test_panel_setitem(self):

        with catch_warnings(record=True):
            # GH 7763
            # loc and setitem have setting differences
            np.random.seed(0)
            index = range(3)
            columns = list('abc')

            panel = Panel({'A': DataFrame(np.random.randn(3, 3),
                                          index=index, columns=columns),
                           'B': DataFrame(np.random.randn(3, 3),
                                          index=index, columns=columns),
                           'C': DataFrame(np.random.randn(3, 3),
                                          index=index, columns=columns)})

            replace = DataFrame(np.eye(3, 3), index=range(3), columns=columns)
            expected = Panel({'A': replace, 'B': replace, 'C': replace})

            p = panel.copy()
            for idx in list('ABC'):
                p[idx] = replace
            tm.assert_panel_equal(p, expected)

            p = panel.copy()
            for idx in list('ABC'):
                p.loc[idx, :, :] = replace
            tm.assert_panel_equal(p, expected)

    def test_panel_assignment(self):

        with catch_warnings(record=True):
            # GH3777
            wp = Panel(np.random.randn(2, 5, 4), items=['Item1', 'Item2'],
                       major_axis=date_range('1/1/2000', periods=5),
                       minor_axis=['A', 'B', 'C', 'D'])
            wp2 = Panel(np.random.randn(2, 5, 4), items=['Item1', 'Item2'],
                        major_axis=date_range('1/1/2000', periods=5),
                        minor_axis=['A', 'B', 'C', 'D'])

            # TODO: unused?
            # expected = wp.loc[['Item1', 'Item2'], :, ['A', 'B']]

            def f():
                wp.loc[['Item1', 'Item2'], :, ['A', 'B']] = wp2.loc[
                    ['Item1', 'Item2'], :, ['A', 'B']]

            pytest.raises(NotImplementedError, f)

            # to_assign = wp2.loc[['Item1', 'Item2'], :, ['A', 'B']]
            # wp.loc[['Item1', 'Item2'], :, ['A', 'B']] = to_assign
            # result = wp.loc[['Item1', 'Item2'], :, ['A', 'B']]
            # tm.assert_panel_equal(result,expected)
