# -*- coding: utf-8 -*-
"""Catch all for categorical functions"""
from __future__ import absolute_import, division, print_function

import pytest
import numpy as np

from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import matplotlib.category as cat

# Python2/3 text handling
_to_str = cat.StrCategoryFormatter._text


class TestUnitData(object):
    test_cases = [('single', (["hello world"], [0])),
                  ('unicode', (["Здравствуйте мир"], [0])),
                  ('mixed', (['A', "np.nan", 'B', "3.14", "мир"],
                             [0, 1, 2, 3, 4]))]
    ids, data = zip(*test_cases)

    @pytest.mark.parametrize("data, locs", data, ids=ids)
    def test_unit(self, data, locs):
        unit = cat.UnitData(data)
        assert list(unit._mapping.keys()) == data
        assert list(unit._mapping.values()) == locs

    def test_update(self):
        data = ['a', 'd']
        locs = [0, 1]

        data_update = ['b', 'd', 'e']
        unique_data = ['a', 'd', 'b', 'e']
        updated_locs = [0, 1, 2, 3]

        unit = cat.UnitData(data)
        assert list(unit._mapping.keys()) == data
        assert list(unit._mapping.values()) == locs

        unit.update(data_update)
        assert list(unit._mapping.keys()) == unique_data
        assert list(unit._mapping.values()) == updated_locs

    failing_test_cases = [("number", 3.14), ("nan", np.nan),
                          ("list", [3.14, 12]), ("mixed type", ["A", 2])]

    fids, fdata = zip(*test_cases)

    @pytest.mark.parametrize("fdata", fdata, ids=fids)
    def test_non_string_fails(self, fdata):
        with pytest.raises(TypeError):
            cat.UnitData(fdata)

    @pytest.mark.parametrize("fdata", fdata, ids=fids)
    def test_non_string_update_fails(self, fdata):
        unitdata = cat.UnitData()
        with pytest.raises(TypeError):
            unitdata.update(fdata)


class FakeAxis(object):
    def __init__(self, units):
        self.units = units


class TestStrCategoryConverter(object):
    """Based on the pandas conversion and factorization tests:

    ref: /pandas/tseries/tests/test_converter.py
         /pandas/tests/test_algos.py:TestFactorize
    """
    test_cases = [("unicode", ["Здравствуйте мир"]),
                  ("ascii", ["hello world"]),
                  ("single", ['a', 'b', 'c']),
                  ("integer string", ["1", "2"]),
                  ("single + values>10", ["A", "B", "C", "D", "E", "F", "G",
                                          "H", "I", "J", "K", "L", "M", "N",
                                          "O", "P", "Q", "R", "S", "T", "U",
                                          "V", "W", "X", "Y", "Z"])]

    ids, values = zip(*test_cases)

    failing_test_cases = [("mixed", [3.14, 'A', np.inf]),
                          ("string integer", ['42', 42])]

    fids, fvalues = zip(*failing_test_cases)

    @pytest.fixture(autouse=True)
    def mock_axis(self, request):
        self.cc = cat.StrCategoryConverter()
        # self.unit should be probably be replaced with real mock unit
        self.unit = cat.UnitData()
        self.ax = FakeAxis(self.unit)

    @pytest.mark.parametrize("vals", values, ids=ids)
    def test_convert(self, vals):
        np.testing.assert_allclose(self.cc.convert(vals, self.ax.units,
                                                   self.ax),
                                   range(len(vals)))

    @pytest.mark.parametrize("value", ["hi", "мир"], ids=["ascii", "unicode"])
    def test_convert_one_string(self, value):
        assert self.cc.convert(value, self.unit, self.ax) == 0

    def test_convert_one_number(self):
        actual = self.cc.convert(0.0, self.unit, self.ax)
        np.testing.assert_allclose(actual, np.array([0.]))

    def test_convert_float_array(self):
        data = np.array([1, 2, 3], dtype=float)
        actual = self.cc.convert(data, self.unit, self.ax)
        np.testing.assert_allclose(actual, np.array([1., 2., 3.]))

    @pytest.mark.parametrize("fvals", fvalues, ids=fids)
    def test_convert_fail(self, fvals):
        with pytest.raises(TypeError):
            self.cc.convert(fvals, self.unit, self.ax)

    def test_axisinfo(self):
        axis = self.cc.axisinfo(self.unit, self.ax)
        assert isinstance(axis.majloc, cat.StrCategoryLocator)
        assert isinstance(axis.majfmt, cat.StrCategoryFormatter)

    def test_default_units(self):
        assert isinstance(self.cc.default_units(["a"], self.ax), cat.UnitData)


@pytest.fixture
def ax():
    return plt.figure().subplots()


PLOT_LIST = [Axes.scatter, Axes.plot, Axes.bar]
PLOT_IDS = ["scatter", "plot", "bar"]


class TestStrCategoryLocator(object):
    def test_StrCategoryLocator(self):
        locs = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        unit = cat.UnitData([str(j) for j in locs])
        ticks = cat.StrCategoryLocator(unit._mapping)
        np.testing.assert_array_equal(ticks.tick_values(None, None), locs)

    @pytest.mark.parametrize("plotter", PLOT_LIST, ids=PLOT_IDS)
    def test_StrCategoryLocatorPlot(self, ax, plotter):
        ax.plot(["a", "b", "c"])
        np.testing.assert_array_equal(ax.yaxis.major.locator(), range(3))


class TestStrCategoryFormatter(object):
    test_cases = [("ascii", ["hello", "world", "hi"]),
                  ("unicode", ["Здравствуйте", "привет"])]

    ids, cases = zip(*test_cases)

    @pytest.mark.parametrize("ydata", cases, ids=ids)
    def test_StrCategoryFormatter(self, ax, ydata):
        unit = cat.UnitData(ydata)
        labels = cat.StrCategoryFormatter(unit._mapping)
        for i, d in enumerate(ydata):
            assert labels(i, i) == _to_str(d)

    @pytest.mark.parametrize("ydata", cases, ids=ids)
    @pytest.mark.parametrize("plotter", PLOT_LIST, ids=PLOT_IDS)
    def test_StrCategoryFormatterPlot(self, ax, ydata, plotter):
        plotter(ax, range(len(ydata)), ydata)
        for i, d in enumerate(ydata):
            assert ax.yaxis.major.formatter(i, i) == _to_str(d)
        assert ax.yaxis.major.formatter(i+1, i+1) == ""
        assert ax.yaxis.major.formatter(0, None) == ""


def axis_test(axis, labels):
    ticks = list(range(len(labels)))
    np.testing.assert_array_equal(axis.get_majorticklocs(), ticks)
    graph_labels = [axis.major.formatter(i, i) for i in ticks]
    assert graph_labels == [_to_str(l) for l in labels]
    assert list(axis.units._mapping.keys()) == [l for l in labels]
    assert list(axis.units._mapping.values()) == ticks


class TestPlotBytes(object):
    bytes_cases = [('string list', ['a', 'b', 'c']),
                   ('bytes list', [b'a', b'b', b'c']),
                   ('bytes ndarray', np.array([b'a', b'b', b'c']))]

    bytes_ids, bytes_data = zip(*bytes_cases)

    @pytest.mark.parametrize("plotter", PLOT_LIST, ids=PLOT_IDS)
    @pytest.mark.parametrize("bdata", bytes_data, ids=bytes_ids)
    def test_plot_bytes(self, ax, plotter, bdata):
        counts = np.array([4, 6, 5])
        plotter(ax, bdata, counts)
        axis_test(ax.xaxis, bdata)


class TestPlotNumlike(object):
    numlike_cases = [('string list', ['1', '11', '3']),
                     ('string ndarray', np.array(['1', '11', '3'])),
                     ('bytes list', [b'1', b'11', b'3']),
                     ('bytes ndarray', np.array([b'1', b'11', b'3']))]
    numlike_ids, numlike_data = zip(*numlike_cases)

    @pytest.mark.parametrize("plotter", PLOT_LIST, ids=PLOT_IDS)
    @pytest.mark.parametrize("ndata", numlike_data, ids=numlike_ids)
    def test_plot_numlike(self, ax, plotter, ndata):
        counts = np.array([4, 6, 5])
        plotter(ax, ndata, counts)
        axis_test(ax.xaxis, ndata)


class TestPlotTypes(object):
    @pytest.mark.parametrize("plotter", PLOT_LIST, ids=PLOT_IDS)
    def test_plot_unicode(self, ax, plotter):
        words = ['Здравствуйте', 'привет']
        plotter(ax, words, [0, 1])
        axis_test(ax.xaxis, words)

    @pytest.fixture
    def test_data(self):
        self.x = ["hello", "happy", "world"]
        self.xy = [2, 6, 3]
        self.y = ["Python", "is", "fun"]
        self.yx = [3, 4, 5]

    @pytest.mark.usefixtures("test_data")
    @pytest.mark.parametrize("plotter", PLOT_LIST, ids=PLOT_IDS)
    def test_plot_xaxis(self, ax, test_data, plotter):
        plotter(ax, self.x, self.xy)
        axis_test(ax.xaxis, self.x)

    @pytest.mark.usefixtures("test_data")
    @pytest.mark.parametrize("plotter", PLOT_LIST, ids=PLOT_IDS)
    def test_plot_yaxis(self, ax, test_data, plotter):
        plotter(ax, self.yx, self.y)
        axis_test(ax.yaxis, self.y)

    @pytest.mark.usefixtures("test_data")
    @pytest.mark.parametrize("plotter", PLOT_LIST, ids=PLOT_IDS)
    def test_plot_xyaxis(self, ax, test_data, plotter):
        plotter(ax, self.x, self.y)
        axis_test(ax.xaxis, self.x)
        axis_test(ax.yaxis, self.y)

    @pytest.mark.parametrize("plotter", PLOT_LIST, ids=PLOT_IDS)
    def test_update_plot(self, ax, plotter):
        plotter(ax, ['a', 'b'], ['e', 'g'])
        plotter(ax, ['a', 'b', 'd'], ['f', 'a', 'b'])
        plotter(ax, ['b', 'c', 'd'], ['g', 'e', 'd'])
        axis_test(ax.xaxis, ['a', 'b', 'd', 'c'])
        axis_test(ax.yaxis, ['e', 'g', 'f', 'a', 'b', 'd'])

    failing_test_cases = [("mixed", ['A', 3.14]),
                          ("number integer", ['1', 1]),
                          ("string integer", ['42', 42]),
                          ("missing", ['12', np.nan])]

    fids, fvalues = zip(*failing_test_cases)

    PLOT_BROKEN_LIST = [Axes.scatter,
                        pytest.param(Axes.plot, marks=pytest.mark.xfail),
                        pytest.param(Axes.bar, marks=pytest.mark.xfail)]

    PLOT_BROKEN_IDS = ["scatter", "plot", "bar"]

    @pytest.mark.parametrize("plotter", PLOT_BROKEN_LIST, ids=PLOT_BROKEN_IDS)
    @pytest.mark.parametrize("xdata", fvalues, ids=fids)
    def test_mixed_type_exception(self, ax, plotter, xdata):
        with pytest.raises(TypeError):
            plotter(ax, xdata, [1, 2])

    @pytest.mark.parametrize("plotter", PLOT_BROKEN_LIST, ids=PLOT_BROKEN_IDS)
    @pytest.mark.parametrize("xdata", fvalues, ids=fids)
    def test_mixed_type_update_exception(self, ax, plotter, xdata):
        with pytest.raises(TypeError):
            plotter(ax, [0, 3], [1, 3])
            plotter(ax, xdata, [1, 2])
