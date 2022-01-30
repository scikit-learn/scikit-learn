import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import pytest
import numpy.testing as npt
from numpy.testing import assert_array_equal
try:
    import pandas.testing as tm
except ImportError:
    import pandas.util.testing as tm

from .._core import categorical_order
from .. import rcmod
from ..palettes import color_palette
from ..relational import scatterplot
from ..distributions import histplot, kdeplot, distplot
from ..categorical import pointplot
from .. import axisgrid as ag
from .._testing import (
    assert_plots_equal,
    assert_colors_equal,
)

rs = np.random.RandomState(0)


class TestFacetGrid:

    df = pd.DataFrame(dict(x=rs.normal(size=60),
                           y=rs.gamma(4, size=60),
                           a=np.repeat(list("abc"), 20),
                           b=np.tile(list("mn"), 30),
                           c=np.tile(list("tuv"), 20),
                           d=np.tile(list("abcdefghijkl"), 5)))

    def test_self_data(self):

        g = ag.FacetGrid(self.df)
        assert g.data is self.df

    def test_self_figure(self):

        g = ag.FacetGrid(self.df)
        assert isinstance(g.figure, plt.Figure)
        assert g.figure is g._figure

    def test_self_axes(self):

        g = ag.FacetGrid(self.df, row="a", col="b", hue="c")
        for ax in g.axes.flat:
            assert isinstance(ax, plt.Axes)

    def test_axes_array_size(self):

        g = ag.FacetGrid(self.df)
        assert g.axes.shape == (1, 1)

        g = ag.FacetGrid(self.df, row="a")
        assert g.axes.shape == (3, 1)

        g = ag.FacetGrid(self.df, col="b")
        assert g.axes.shape == (1, 2)

        g = ag.FacetGrid(self.df, hue="c")
        assert g.axes.shape == (1, 1)

        g = ag.FacetGrid(self.df, row="a", col="b", hue="c")
        assert g.axes.shape == (3, 2)
        for ax in g.axes.flat:
            assert isinstance(ax, plt.Axes)

    def test_single_axes(self):

        g = ag.FacetGrid(self.df)
        assert isinstance(g.ax, plt.Axes)

        g = ag.FacetGrid(self.df, row="a")
        with pytest.raises(AttributeError):
            g.ax

        g = ag.FacetGrid(self.df, col="a")
        with pytest.raises(AttributeError):
            g.ax

        g = ag.FacetGrid(self.df, col="a", row="b")
        with pytest.raises(AttributeError):
            g.ax

    def test_col_wrap(self):

        n = len(self.df.d.unique())

        g = ag.FacetGrid(self.df, col="d")
        assert g.axes.shape == (1, n)
        assert g.facet_axis(0, 8) is g.axes[0, 8]

        g_wrap = ag.FacetGrid(self.df, col="d", col_wrap=4)
        assert g_wrap.axes.shape == (n,)
        assert g_wrap.facet_axis(0, 8) is g_wrap.axes[8]
        assert g_wrap._ncol == 4
        assert g_wrap._nrow == (n / 4)

        with pytest.raises(ValueError):
            g = ag.FacetGrid(self.df, row="b", col="d", col_wrap=4)

        df = self.df.copy()
        df.loc[df.d == "j"] = np.nan
        g_missing = ag.FacetGrid(df, col="d")
        assert g_missing.axes.shape == (1, n - 1)

        g_missing_wrap = ag.FacetGrid(df, col="d", col_wrap=4)
        assert g_missing_wrap.axes.shape == (n - 1,)

        g = ag.FacetGrid(self.df, col="d", col_wrap=1)
        assert len(list(g.facet_data())) == n

    def test_normal_axes(self):

        null = np.empty(0, object).flat

        g = ag.FacetGrid(self.df)
        npt.assert_array_equal(g._bottom_axes, g.axes.flat)
        npt.assert_array_equal(g._not_bottom_axes, null)
        npt.assert_array_equal(g._left_axes, g.axes.flat)
        npt.assert_array_equal(g._not_left_axes, null)
        npt.assert_array_equal(g._inner_axes, null)

        g = ag.FacetGrid(self.df, col="c")
        npt.assert_array_equal(g._bottom_axes, g.axes.flat)
        npt.assert_array_equal(g._not_bottom_axes, null)
        npt.assert_array_equal(g._left_axes, g.axes[:, 0].flat)
        npt.assert_array_equal(g._not_left_axes, g.axes[:, 1:].flat)
        npt.assert_array_equal(g._inner_axes, null)

        g = ag.FacetGrid(self.df, row="c")
        npt.assert_array_equal(g._bottom_axes, g.axes[-1, :].flat)
        npt.assert_array_equal(g._not_bottom_axes, g.axes[:-1, :].flat)
        npt.assert_array_equal(g._left_axes, g.axes.flat)
        npt.assert_array_equal(g._not_left_axes, null)
        npt.assert_array_equal(g._inner_axes, null)

        g = ag.FacetGrid(self.df, col="a", row="c")
        npt.assert_array_equal(g._bottom_axes, g.axes[-1, :].flat)
        npt.assert_array_equal(g._not_bottom_axes, g.axes[:-1, :].flat)
        npt.assert_array_equal(g._left_axes, g.axes[:, 0].flat)
        npt.assert_array_equal(g._not_left_axes, g.axes[:, 1:].flat)
        npt.assert_array_equal(g._inner_axes, g.axes[:-1, 1:].flat)

    def test_wrapped_axes(self):

        null = np.empty(0, object).flat

        g = ag.FacetGrid(self.df, col="a", col_wrap=2)
        npt.assert_array_equal(g._bottom_axes,
                               g.axes[np.array([1, 2])].flat)
        npt.assert_array_equal(g._not_bottom_axes, g.axes[:1].flat)
        npt.assert_array_equal(g._left_axes, g.axes[np.array([0, 2])].flat)
        npt.assert_array_equal(g._not_left_axes, g.axes[np.array([1])].flat)
        npt.assert_array_equal(g._inner_axes, null)

    def test_axes_dict(self):

        g = ag.FacetGrid(self.df)
        assert isinstance(g.axes_dict, dict)
        assert not g.axes_dict

        g = ag.FacetGrid(self.df, row="c")
        assert list(g.axes_dict.keys()) == g.row_names
        for (name, ax) in zip(g.row_names, g.axes.flat):
            assert g.axes_dict[name] is ax

        g = ag.FacetGrid(self.df, col="c")
        assert list(g.axes_dict.keys()) == g.col_names
        for (name, ax) in zip(g.col_names, g.axes.flat):
            assert g.axes_dict[name] is ax

        g = ag.FacetGrid(self.df, col="a", col_wrap=2)
        assert list(g.axes_dict.keys()) == g.col_names
        for (name, ax) in zip(g.col_names, g.axes.flat):
            assert g.axes_dict[name] is ax

        g = ag.FacetGrid(self.df, row="a", col="c")
        for (row_var, col_var), ax in g.axes_dict.items():
            i = g.row_names.index(row_var)
            j = g.col_names.index(col_var)
            assert g.axes[i, j] is ax

    def test_figure_size(self):

        g = ag.FacetGrid(self.df, row="a", col="b")
        npt.assert_array_equal(g.figure.get_size_inches(), (6, 9))

        g = ag.FacetGrid(self.df, row="a", col="b", height=6)
        npt.assert_array_equal(g.figure.get_size_inches(), (12, 18))

        g = ag.FacetGrid(self.df, col="c", height=4, aspect=.5)
        npt.assert_array_equal(g.figure.get_size_inches(), (6, 4))

    def test_figure_size_with_legend(self):

        g = ag.FacetGrid(self.df, col="a", hue="c", height=4, aspect=.5)
        npt.assert_array_equal(g.figure.get_size_inches(), (6, 4))
        g.add_legend()
        assert g.figure.get_size_inches()[0] > 6

        g = ag.FacetGrid(self.df, col="a", hue="c", height=4, aspect=.5,
                         legend_out=False)
        npt.assert_array_equal(g.figure.get_size_inches(), (6, 4))
        g.add_legend()
        npt.assert_array_equal(g.figure.get_size_inches(), (6, 4))

    def test_legend_data(self):

        g = ag.FacetGrid(self.df, hue="a")
        g.map(plt.plot, "x", "y")
        g.add_legend()
        palette = color_palette(n_colors=3)

        assert g._legend.get_title().get_text() == "a"

        a_levels = sorted(self.df.a.unique())

        lines = g._legend.get_lines()
        assert len(lines) == len(a_levels)

        for line, hue in zip(lines, palette):
            assert line.get_color() == hue

        labels = g._legend.get_texts()
        assert len(labels) == len(a_levels)

        for label, level in zip(labels, a_levels):
            assert label.get_text() == level

    def test_legend_data_missing_level(self):

        g = ag.FacetGrid(self.df, hue="a", hue_order=list("azbc"))
        g.map(plt.plot, "x", "y")
        g.add_legend()

        c1, c2, c3, c4 = color_palette(n_colors=4)
        palette = [c1, c3, c4]

        assert g._legend.get_title().get_text() == "a"

        a_levels = sorted(self.df.a.unique())

        lines = g._legend.get_lines()
        assert len(lines) == len(a_levels)

        for line, hue in zip(lines, palette):
            assert line.get_color() == hue

        labels = g._legend.get_texts()
        assert len(labels) == 4

        for label, level in zip(labels, list("azbc")):
            assert label.get_text() == level

    def test_get_boolean_legend_data(self):

        self.df["b_bool"] = self.df.b == "m"
        g = ag.FacetGrid(self.df, hue="b_bool")
        g.map(plt.plot, "x", "y")
        g.add_legend()
        palette = color_palette(n_colors=2)

        assert g._legend.get_title().get_text() == "b_bool"

        b_levels = list(map(str, categorical_order(self.df.b_bool)))

        lines = g._legend.get_lines()
        assert len(lines) == len(b_levels)

        for line, hue in zip(lines, palette):
            assert line.get_color() == hue

        labels = g._legend.get_texts()
        assert len(labels) == len(b_levels)

        for label, level in zip(labels, b_levels):
            assert label.get_text() == level

    def test_legend_tuples(self):

        g = ag.FacetGrid(self.df, hue="a")
        g.map(plt.plot, "x", "y")

        handles, labels = g.ax.get_legend_handles_labels()
        label_tuples = [("", l) for l in labels]
        legend_data = dict(zip(label_tuples, handles))
        g.add_legend(legend_data, label_tuples)
        for entry, label in zip(g._legend.get_texts(), labels):
            assert entry.get_text() == label

    def test_legend_options(self):

        g = ag.FacetGrid(self.df, hue="b")
        g.map(plt.plot, "x", "y")
        g.add_legend()

        g1 = ag.FacetGrid(self.df, hue="b", legend_out=False)
        g1.add_legend(adjust_subtitles=True)

        g1 = ag.FacetGrid(self.df, hue="b", legend_out=False)
        g1.add_legend(adjust_subtitles=False)

    def test_legendout_with_colwrap(self):

        g = ag.FacetGrid(self.df, col="d", hue='b',
                         col_wrap=4, legend_out=False)
        g.map(plt.plot, "x", "y", linewidth=3)
        g.add_legend()

    def test_legend_tight_layout(self):

        g = ag.FacetGrid(self.df, hue='b')
        g.map(plt.plot, "x", "y", linewidth=3)
        g.add_legend()
        g.tight_layout()

        axes_right_edge = g.ax.get_window_extent().xmax
        legend_left_edge = g._legend.get_window_extent().xmin

        assert axes_right_edge < legend_left_edge

    def test_subplot_kws(self):

        g = ag.FacetGrid(self.df, despine=False,
                         subplot_kws=dict(projection="polar"))
        for ax in g.axes.flat:
            assert "PolarAxesSubplot" in str(type(ax))

    def test_gridspec_kws(self):
        ratios = [3, 1, 2]

        gskws = dict(width_ratios=ratios)
        g = ag.FacetGrid(self.df, col='c', row='a', gridspec_kws=gskws)

        for ax in g.axes.flat:
            ax.set_xticks([])
            ax.set_yticks([])

        g.figure.tight_layout()

        for (l, m, r) in g.axes:
            assert l.get_position().width > m.get_position().width
            assert r.get_position().width > m.get_position().width

    def test_gridspec_kws_col_wrap(self):
        ratios = [3, 1, 2, 1, 1]

        gskws = dict(width_ratios=ratios)
        with pytest.warns(UserWarning):
            ag.FacetGrid(self.df, col='d', col_wrap=5, gridspec_kws=gskws)

    def test_data_generator(self):

        g = ag.FacetGrid(self.df, row="a")
        d = list(g.facet_data())
        assert len(d) == 3

        tup, data = d[0]
        assert tup == (0, 0, 0)
        assert (data["a"] == "a").all()

        tup, data = d[1]
        assert tup == (1, 0, 0)
        assert (data["a"] == "b").all()

        g = ag.FacetGrid(self.df, row="a", col="b")
        d = list(g.facet_data())
        assert len(d) == 6

        tup, data = d[0]
        assert tup == (0, 0, 0)
        assert (data["a"] == "a").all()
        assert (data["b"] == "m").all()

        tup, data = d[1]
        assert tup == (0, 1, 0)
        assert (data["a"] == "a").all()
        assert (data["b"] == "n").all()

        tup, data = d[2]
        assert tup == (1, 0, 0)
        assert (data["a"] == "b").all()
        assert (data["b"] == "m").all()

        g = ag.FacetGrid(self.df, hue="c")
        d = list(g.facet_data())
        assert len(d) == 3
        tup, data = d[1]
        assert tup == (0, 0, 1)
        assert (data["c"] == "u").all()

    def test_map(self):

        g = ag.FacetGrid(self.df, row="a", col="b", hue="c")
        g.map(plt.plot, "x", "y", linewidth=3)

        lines = g.axes[0, 0].lines
        assert len(lines) == 3

        line1, _, _ = lines
        assert line1.get_linewidth() == 3
        x, y = line1.get_data()
        mask = (self.df.a == "a") & (self.df.b == "m") & (self.df.c == "t")
        npt.assert_array_equal(x, self.df.x[mask])
        npt.assert_array_equal(y, self.df.y[mask])

    def test_map_dataframe(self):

        g = ag.FacetGrid(self.df, row="a", col="b", hue="c")

        def plot(x, y, data=None, **kws):
            plt.plot(data[x], data[y], **kws)
        # Modify __module__ so this doesn't look like a seaborn function
        plot.__module__ = "test"

        g.map_dataframe(plot, "x", "y", linestyle="--")

        lines = g.axes[0, 0].lines
        assert len(g.axes[0, 0].lines) == 3

        line1, _, _ = lines
        assert line1.get_linestyle() == "--"
        x, y = line1.get_data()
        mask = (self.df.a == "a") & (self.df.b == "m") & (self.df.c == "t")
        npt.assert_array_equal(x, self.df.x[mask])
        npt.assert_array_equal(y, self.df.y[mask])

    def test_set(self):

        g = ag.FacetGrid(self.df, row="a", col="b")
        xlim = (-2, 5)
        ylim = (3, 6)
        xticks = [-2, 0, 3, 5]
        yticks = [3, 4.5, 6]
        g.set(xlim=xlim, ylim=ylim, xticks=xticks, yticks=yticks)
        for ax in g.axes.flat:
            npt.assert_array_equal(ax.get_xlim(), xlim)
            npt.assert_array_equal(ax.get_ylim(), ylim)
            npt.assert_array_equal(ax.get_xticks(), xticks)
            npt.assert_array_equal(ax.get_yticks(), yticks)

    def test_set_titles(self):

        g = ag.FacetGrid(self.df, row="a", col="b")
        g.map(plt.plot, "x", "y")

        # Test the default titles
        assert g.axes[0, 0].get_title() == "a = a | b = m"
        assert g.axes[0, 1].get_title() == "a = a | b = n"
        assert g.axes[1, 0].get_title() == "a = b | b = m"

        # Test a provided title
        g.set_titles("{row_var} == {row_name} \\/ {col_var} == {col_name}")
        assert g.axes[0, 0].get_title() == "a == a \\/ b == m"
        assert g.axes[0, 1].get_title() == "a == a \\/ b == n"
        assert g.axes[1, 0].get_title() == "a == b \\/ b == m"

        # Test a single row
        g = ag.FacetGrid(self.df, col="b")
        g.map(plt.plot, "x", "y")

        # Test the default titles
        assert g.axes[0, 0].get_title() == "b = m"
        assert g.axes[0, 1].get_title() == "b = n"

        # test with dropna=False
        g = ag.FacetGrid(self.df, col="b", hue="b", dropna=False)
        g.map(plt.plot, 'x', 'y')

    def test_set_titles_margin_titles(self):

        g = ag.FacetGrid(self.df, row="a", col="b", margin_titles=True)
        g.map(plt.plot, "x", "y")

        # Test the default titles
        assert g.axes[0, 0].get_title() == "b = m"
        assert g.axes[0, 1].get_title() == "b = n"
        assert g.axes[1, 0].get_title() == ""

        # Test the row "titles"
        assert g.axes[0, 1].texts[0].get_text() == "a = a"
        assert g.axes[1, 1].texts[0].get_text() == "a = b"
        assert g.axes[0, 1].texts[0] is g._margin_titles_texts[0]

        # Test provided titles
        g.set_titles(col_template="{col_name}", row_template="{row_name}")
        assert g.axes[0, 0].get_title() == "m"
        assert g.axes[0, 1].get_title() == "n"
        assert g.axes[1, 0].get_title() == ""

        assert len(g.axes[1, 1].texts) == 1
        assert g.axes[1, 1].texts[0].get_text() == "b"

    def test_set_ticklabels(self):

        g = ag.FacetGrid(self.df, row="a", col="b")
        g.map(plt.plot, "x", "y")

        ax = g.axes[-1, 0]
        xlab = [l.get_text() + "h" for l in ax.get_xticklabels()]
        ylab = [l.get_text() + "i" for l in ax.get_yticklabels()]

        g.set_xticklabels(xlab)
        g.set_yticklabels(ylab)
        got_x = [l.get_text() for l in g.axes[-1, 1].get_xticklabels()]
        got_y = [l.get_text() for l in g.axes[0, 0].get_yticklabels()]
        npt.assert_array_equal(got_x, xlab)
        npt.assert_array_equal(got_y, ylab)

        x, y = np.arange(10), np.arange(10)
        df = pd.DataFrame(np.c_[x, y], columns=["x", "y"])
        g = ag.FacetGrid(df).map_dataframe(pointplot, x="x", y="y", order=x)
        g.set_xticklabels(step=2)
        got_x = [int(l.get_text()) for l in g.axes[0, 0].get_xticklabels()]
        npt.assert_array_equal(x[::2], got_x)

        g = ag.FacetGrid(self.df, col="d", col_wrap=5)
        g.map(plt.plot, "x", "y")
        g.set_xticklabels(rotation=45)
        g.set_yticklabels(rotation=75)
        for ax in g._bottom_axes:
            for l in ax.get_xticklabels():
                assert l.get_rotation() == 45
        for ax in g._left_axes:
            for l in ax.get_yticklabels():
                assert l.get_rotation() == 75

    def test_set_axis_labels(self):

        g = ag.FacetGrid(self.df, row="a", col="b")
        g.map(plt.plot, "x", "y")
        xlab = 'xx'
        ylab = 'yy'

        g.set_axis_labels(xlab, ylab)

        got_x = [ax.get_xlabel() for ax in g.axes[-1, :]]
        got_y = [ax.get_ylabel() for ax in g.axes[:, 0]]
        npt.assert_array_equal(got_x, xlab)
        npt.assert_array_equal(got_y, ylab)

        for ax in g.axes.flat:
            ax.set(xlabel="x", ylabel="y")

        g.set_axis_labels(xlab, ylab)
        for ax in g._not_bottom_axes:
            assert not ax.get_xlabel()
        for ax in g._not_left_axes:
            assert not ax.get_ylabel()

    def test_axis_lims(self):

        g = ag.FacetGrid(self.df, row="a", col="b", xlim=(0, 4), ylim=(-2, 3))
        assert g.axes[0, 0].get_xlim() == (0, 4)
        assert g.axes[0, 0].get_ylim() == (-2, 3)

    def test_data_orders(self):

        g = ag.FacetGrid(self.df, row="a", col="b", hue="c")

        assert g.row_names == list("abc")
        assert g.col_names == list("mn")
        assert g.hue_names == list("tuv")
        assert g.axes.shape == (3, 2)

        g = ag.FacetGrid(self.df, row="a", col="b", hue="c",
                         row_order=list("bca"),
                         col_order=list("nm"),
                         hue_order=list("vtu"))

        assert g.row_names == list("bca")
        assert g.col_names == list("nm")
        assert g.hue_names == list("vtu")
        assert g.axes.shape == (3, 2)

        g = ag.FacetGrid(self.df, row="a", col="b", hue="c",
                         row_order=list("bcda"),
                         col_order=list("nom"),
                         hue_order=list("qvtu"))

        assert g.row_names == list("bcda")
        assert g.col_names == list("nom")
        assert g.hue_names == list("qvtu")
        assert g.axes.shape == (4, 3)

    def test_palette(self):

        rcmod.set()

        g = ag.FacetGrid(self.df, hue="c")
        assert g._colors == color_palette(n_colors=len(self.df.c.unique()))

        g = ag.FacetGrid(self.df, hue="d")
        assert g._colors == color_palette("husl", len(self.df.d.unique()))

        g = ag.FacetGrid(self.df, hue="c", palette="Set2")
        assert g._colors == color_palette("Set2", len(self.df.c.unique()))

        dict_pal = dict(t="red", u="green", v="blue")
        list_pal = color_palette(["red", "green", "blue"], 3)
        g = ag.FacetGrid(self.df, hue="c", palette=dict_pal)
        assert g._colors == list_pal

        list_pal = color_palette(["green", "blue", "red"], 3)
        g = ag.FacetGrid(self.df, hue="c", hue_order=list("uvt"),
                         palette=dict_pal)
        assert g._colors == list_pal

    def test_hue_kws(self):

        kws = dict(marker=["o", "s", "D"])
        g = ag.FacetGrid(self.df, hue="c", hue_kws=kws)
        g.map(plt.plot, "x", "y")

        for line, marker in zip(g.axes[0, 0].lines, kws["marker"]):
            assert line.get_marker() == marker

    def test_dropna(self):

        df = self.df.copy()
        hasna = pd.Series(np.tile(np.arange(6), 10), dtype=float)
        hasna[hasna == 5] = np.nan
        df["hasna"] = hasna
        g = ag.FacetGrid(df, dropna=False, row="hasna")
        assert g._not_na.sum() == 60

        g = ag.FacetGrid(df, dropna=True, row="hasna")
        assert g._not_na.sum() == 50

    def test_categorical_column_missing_categories(self):

        df = self.df.copy()
        df['a'] = df['a'].astype('category')

        g = ag.FacetGrid(df[df['a'] == 'a'], col="a", col_wrap=1)

        assert g.axes.shape == (len(df['a'].cat.categories),)

    def test_categorical_warning(self):

        g = ag.FacetGrid(self.df, col="b")
        with pytest.warns(UserWarning):
            g.map(pointplot, "b", "x")

    def test_refline(self):

        g = ag.FacetGrid(self.df, row="a", col="b")
        g.refline()
        for ax in g.axes.ravel():
            assert not ax.lines

        refx = refy = 0.5
        hline = np.array([[0, refy], [1, refy]])
        vline = np.array([[refx, 0], [refx, 1]])
        g.refline(x=refx, y=refy)
        for ax in g.axes.ravel():
            assert ax.lines[0].get_color() == '.5'
            assert ax.lines[0].get_linestyle() == '--'
            assert len(ax.lines) == 2
            npt.assert_array_equal(ax.lines[0].get_xydata(), vline)
            npt.assert_array_equal(ax.lines[1].get_xydata(), hline)

        color, linestyle = 'red', '-'
        g.refline(x=refx, color=color, linestyle=linestyle)
        npt.assert_array_equal(g.axes[0, 0].lines[-1].get_xydata(), vline)
        assert g.axes[0, 0].lines[-1].get_color() == color
        assert g.axes[0, 0].lines[-1].get_linestyle() == linestyle


class TestPairGrid:

    rs = np.random.RandomState(sum(map(ord, "PairGrid")))
    df = pd.DataFrame(dict(x=rs.normal(size=60),
                           y=rs.randint(0, 4, size=(60)),
                           z=rs.gamma(3, size=60),
                           a=np.repeat(list("abc"), 20),
                           b=np.repeat(list("abcdefghijkl"), 5)))

    def test_self_data(self):

        g = ag.PairGrid(self.df)
        assert g.data is self.df

    def test_ignore_datelike_data(self):

        df = self.df.copy()
        df['date'] = pd.date_range('2010-01-01', periods=len(df), freq='d')
        result = ag.PairGrid(self.df).data
        expected = df.drop('date', axis=1)
        tm.assert_frame_equal(result, expected)

    def test_self_figure(self):

        g = ag.PairGrid(self.df)
        assert isinstance(g.figure, plt.Figure)
        assert g.figure is g._figure

    def test_self_axes(self):

        g = ag.PairGrid(self.df)
        for ax in g.axes.flat:
            assert isinstance(ax, plt.Axes)

    def test_default_axes(self):

        g = ag.PairGrid(self.df)
        assert g.axes.shape == (3, 3)
        assert g.x_vars == ["x", "y", "z"]
        assert g.y_vars == ["x", "y", "z"]
        assert g.square_grid

    @pytest.mark.parametrize("vars", [["z", "x"], np.array(["z", "x"])])
    def test_specific_square_axes(self, vars):

        g = ag.PairGrid(self.df, vars=vars)
        assert g.axes.shape == (len(vars), len(vars))
        assert g.x_vars == list(vars)
        assert g.y_vars == list(vars)
        assert g.square_grid

    def test_remove_hue_from_default(self):

        hue = "z"
        g = ag.PairGrid(self.df, hue=hue)
        assert hue not in g.x_vars
        assert hue not in g.y_vars

        vars = ["x", "y", "z"]
        g = ag.PairGrid(self.df, hue=hue, vars=vars)
        assert hue in g.x_vars
        assert hue in g.y_vars

    @pytest.mark.parametrize(
        "x_vars, y_vars",
        [
            (["x", "y"], ["z", "y", "x"]),
            (["x", "y"], "z"),
            (np.array(["x", "y"]), np.array(["z", "y", "x"])),
        ],
    )
    def test_specific_nonsquare_axes(self, x_vars, y_vars):

        g = ag.PairGrid(self.df, x_vars=x_vars, y_vars=y_vars)
        assert g.axes.shape == (len(y_vars), len(x_vars))
        assert g.x_vars == list(x_vars)
        assert g.y_vars == list(y_vars)
        assert not g.square_grid

    def test_corner(self):

        plot_vars = ["x", "y", "z"]
        g = ag.PairGrid(self.df, vars=plot_vars, corner=True)
        corner_size = sum([i + 1 for i in range(len(plot_vars))])
        assert len(g.figure.axes) == corner_size

        g.map_diag(plt.hist)
        assert len(g.figure.axes) == (corner_size + len(plot_vars))

        for ax in np.diag(g.axes):
            assert not ax.yaxis.get_visible()
            assert not g.axes[0, 0].get_ylabel()

        plot_vars = ["x", "y", "z"]
        g = ag.PairGrid(self.df, vars=plot_vars, corner=True)
        g.map(scatterplot)
        assert len(g.figure.axes) == corner_size

    def test_size(self):

        g1 = ag.PairGrid(self.df, height=3)
        npt.assert_array_equal(g1.fig.get_size_inches(), (9, 9))

        g2 = ag.PairGrid(self.df, height=4, aspect=.5)
        npt.assert_array_equal(g2.fig.get_size_inches(), (6, 12))

        g3 = ag.PairGrid(self.df, y_vars=["z"], x_vars=["x", "y"],
                         height=2, aspect=2)
        npt.assert_array_equal(g3.fig.get_size_inches(), (8, 2))

    def test_empty_grid(self):

        with pytest.raises(ValueError, match="No variables found"):
            ag.PairGrid(self.df[["a", "b"]])

    def test_map(self):

        vars = ["x", "y", "z"]
        g1 = ag.PairGrid(self.df)
        g1.map(plt.scatter)

        for i, axes_i in enumerate(g1.axes):
            for j, ax in enumerate(axes_i):
                x_in = self.df[vars[j]]
                y_in = self.df[vars[i]]
                x_out, y_out = ax.collections[0].get_offsets().T
                npt.assert_array_equal(x_in, x_out)
                npt.assert_array_equal(y_in, y_out)

        g2 = ag.PairGrid(self.df, hue="a")
        g2.map(plt.scatter)

        for i, axes_i in enumerate(g2.axes):
            for j, ax in enumerate(axes_i):
                x_in = self.df[vars[j]]
                y_in = self.df[vars[i]]
                for k, k_level in enumerate(self.df.a.unique()):
                    x_in_k = x_in[self.df.a == k_level]
                    y_in_k = y_in[self.df.a == k_level]
                    x_out, y_out = ax.collections[k].get_offsets().T
                npt.assert_array_equal(x_in_k, x_out)
                npt.assert_array_equal(y_in_k, y_out)

    def test_map_nonsquare(self):

        x_vars = ["x"]
        y_vars = ["y", "z"]
        g = ag.PairGrid(self.df, x_vars=x_vars, y_vars=y_vars)
        g.map(plt.scatter)

        x_in = self.df.x
        for i, i_var in enumerate(y_vars):
            ax = g.axes[i, 0]
            y_in = self.df[i_var]
            x_out, y_out = ax.collections[0].get_offsets().T
            npt.assert_array_equal(x_in, x_out)
            npt.assert_array_equal(y_in, y_out)

    def test_map_lower(self):

        vars = ["x", "y", "z"]
        g = ag.PairGrid(self.df)
        g.map_lower(plt.scatter)

        for i, j in zip(*np.tril_indices_from(g.axes, -1)):
            ax = g.axes[i, j]
            x_in = self.df[vars[j]]
            y_in = self.df[vars[i]]
            x_out, y_out = ax.collections[0].get_offsets().T
            npt.assert_array_equal(x_in, x_out)
            npt.assert_array_equal(y_in, y_out)

        for i, j in zip(*np.triu_indices_from(g.axes)):
            ax = g.axes[i, j]
            assert len(ax.collections) == 0

    def test_map_upper(self):

        vars = ["x", "y", "z"]
        g = ag.PairGrid(self.df)
        g.map_upper(plt.scatter)

        for i, j in zip(*np.triu_indices_from(g.axes, 1)):
            ax = g.axes[i, j]
            x_in = self.df[vars[j]]
            y_in = self.df[vars[i]]
            x_out, y_out = ax.collections[0].get_offsets().T
            npt.assert_array_equal(x_in, x_out)
            npt.assert_array_equal(y_in, y_out)

        for i, j in zip(*np.tril_indices_from(g.axes)):
            ax = g.axes[i, j]
            assert len(ax.collections) == 0

    def test_map_mixed_funcsig(self):

        vars = ["x", "y", "z"]
        g = ag.PairGrid(self.df, vars=vars)
        g.map_lower(scatterplot)
        g.map_upper(plt.scatter)

        for i, j in zip(*np.triu_indices_from(g.axes, 1)):
            ax = g.axes[i, j]
            x_in = self.df[vars[j]]
            y_in = self.df[vars[i]]
            x_out, y_out = ax.collections[0].get_offsets().T
            npt.assert_array_equal(x_in, x_out)
            npt.assert_array_equal(y_in, y_out)

    def test_map_diag(self):

        g = ag.PairGrid(self.df)
        g.map_diag(plt.hist)

        for var, ax in zip(g.diag_vars, g.diag_axes):
            assert len(ax.patches) == 10
            assert pytest.approx(ax.patches[0].get_x()) == self.df[var].min()

        g = ag.PairGrid(self.df, hue="a")
        g.map_diag(plt.hist)

        for ax in g.diag_axes:
            assert len(ax.patches) == 30

        g = ag.PairGrid(self.df, hue="a")
        g.map_diag(plt.hist, histtype='step')

        for ax in g.diag_axes:
            for ptch in ax.patches:
                assert not ptch.fill

    def test_map_diag_rectangular(self):

        x_vars = ["x", "y"]
        y_vars = ["x", "z", "y"]
        g1 = ag.PairGrid(self.df, x_vars=x_vars, y_vars=y_vars)
        g1.map_diag(plt.hist)
        g1.map_offdiag(plt.scatter)

        assert set(g1.diag_vars) == (set(x_vars) & set(y_vars))

        for var, ax in zip(g1.diag_vars, g1.diag_axes):
            assert len(ax.patches) == 10
            assert pytest.approx(ax.patches[0].get_x()) == self.df[var].min()

        for j, x_var in enumerate(x_vars):
            for i, y_var in enumerate(y_vars):

                ax = g1.axes[i, j]
                if x_var == y_var:
                    diag_ax = g1.diag_axes[j]  # because fewer x than y vars
                    assert ax.bbox.bounds == diag_ax.bbox.bounds

                else:
                    x, y = ax.collections[0].get_offsets().T
                    assert_array_equal(x, self.df[x_var])
                    assert_array_equal(y, self.df[y_var])

        g2 = ag.PairGrid(self.df, x_vars=x_vars, y_vars=y_vars, hue="a")
        g2.map_diag(plt.hist)
        g2.map_offdiag(plt.scatter)

        assert set(g2.diag_vars) == (set(x_vars) & set(y_vars))

        for ax in g2.diag_axes:
            assert len(ax.patches) == 30

        x_vars = ["x", "y", "z"]
        y_vars = ["x", "z"]
        g3 = ag.PairGrid(self.df, x_vars=x_vars, y_vars=y_vars)
        g3.map_diag(plt.hist)
        g3.map_offdiag(plt.scatter)

        assert set(g3.diag_vars) == (set(x_vars) & set(y_vars))

        for var, ax in zip(g3.diag_vars, g3.diag_axes):
            assert len(ax.patches) == 10
            assert pytest.approx(ax.patches[0].get_x()) == self.df[var].min()

        for j, x_var in enumerate(x_vars):
            for i, y_var in enumerate(y_vars):

                ax = g3.axes[i, j]
                if x_var == y_var:
                    diag_ax = g3.diag_axes[i]  # because fewer y than x vars
                    assert ax.bbox.bounds == diag_ax.bbox.bounds
                else:
                    x, y = ax.collections[0].get_offsets().T
                    assert_array_equal(x, self.df[x_var])
                    assert_array_equal(y, self.df[y_var])

    def test_map_diag_color(self):

        color = "red"

        g1 = ag.PairGrid(self.df)
        g1.map_diag(plt.hist, color=color)

        for ax in g1.diag_axes:
            for patch in ax.patches:
                assert_colors_equal(patch.get_facecolor(), color)

        g2 = ag.PairGrid(self.df)
        g2.map_diag(kdeplot, color='red')

        for ax in g2.diag_axes:
            for line in ax.lines:
                assert_colors_equal(line.get_color(), color)

    def test_map_diag_palette(self):

        palette = "muted"
        pal = color_palette(palette, n_colors=len(self.df.a.unique()))
        g = ag.PairGrid(self.df, hue="a", palette=palette)
        g.map_diag(kdeplot)

        for ax in g.diag_axes:
            for line, color in zip(ax.lines[::-1], pal):
                assert_colors_equal(line.get_color(), color)

    def test_map_diag_and_offdiag(self):

        vars = ["x", "y", "z"]
        g = ag.PairGrid(self.df)
        g.map_offdiag(plt.scatter)
        g.map_diag(plt.hist)

        for ax in g.diag_axes:
            assert len(ax.patches) == 10

        for i, j in zip(*np.triu_indices_from(g.axes, 1)):
            ax = g.axes[i, j]
            x_in = self.df[vars[j]]
            y_in = self.df[vars[i]]
            x_out, y_out = ax.collections[0].get_offsets().T
            npt.assert_array_equal(x_in, x_out)
            npt.assert_array_equal(y_in, y_out)

        for i, j in zip(*np.tril_indices_from(g.axes, -1)):
            ax = g.axes[i, j]
            x_in = self.df[vars[j]]
            y_in = self.df[vars[i]]
            x_out, y_out = ax.collections[0].get_offsets().T
            npt.assert_array_equal(x_in, x_out)
            npt.assert_array_equal(y_in, y_out)

        for i, j in zip(*np.diag_indices_from(g.axes)):
            ax = g.axes[i, j]
            assert len(ax.collections) == 0

    def test_diag_sharey(self):

        g = ag.PairGrid(self.df, diag_sharey=True)
        g.map_diag(kdeplot)
        for ax in g.diag_axes[1:]:
            assert ax.get_ylim() == g.diag_axes[0].get_ylim()

    def test_map_diag_matplotlib(self):

        bins = 10
        g = ag.PairGrid(self.df)
        g.map_diag(plt.hist, bins=bins)
        for ax in g.diag_axes:
            assert len(ax.patches) == bins

        levels = len(self.df["a"].unique())
        g = ag.PairGrid(self.df, hue="a")
        g.map_diag(plt.hist, bins=bins)
        for ax in g.diag_axes:
            assert len(ax.patches) == (bins * levels)

    def test_palette(self):

        rcmod.set()

        g = ag.PairGrid(self.df, hue="a")
        assert g.palette == color_palette(n_colors=len(self.df.a.unique()))

        g = ag.PairGrid(self.df, hue="b")
        assert g.palette == color_palette("husl", len(self.df.b.unique()))

        g = ag.PairGrid(self.df, hue="a", palette="Set2")
        assert g.palette == color_palette("Set2", len(self.df.a.unique()))

        dict_pal = dict(a="red", b="green", c="blue")
        list_pal = color_palette(["red", "green", "blue"])
        g = ag.PairGrid(self.df, hue="a", palette=dict_pal)
        assert g.palette == list_pal

        list_pal = color_palette(["blue", "red", "green"])
        g = ag.PairGrid(self.df, hue="a", hue_order=list("cab"),
                        palette=dict_pal)
        assert g.palette == list_pal

    def test_hue_kws(self):

        kws = dict(marker=["o", "s", "d", "+"])
        g = ag.PairGrid(self.df, hue="a", hue_kws=kws)
        g.map(plt.plot)

        for line, marker in zip(g.axes[0, 0].lines, kws["marker"]):
            assert line.get_marker() == marker

        g = ag.PairGrid(self.df, hue="a", hue_kws=kws,
                        hue_order=list("dcab"))
        g.map(plt.plot)

        for line, marker in zip(g.axes[0, 0].lines, kws["marker"]):
            assert line.get_marker() == marker

    def test_hue_order(self):

        order = list("dcab")
        g = ag.PairGrid(self.df, hue="a", hue_order=order)
        g.map(plt.plot)

        for line, level in zip(g.axes[1, 0].lines, order):
            x, y = line.get_xydata().T
            npt.assert_array_equal(x, self.df.loc[self.df.a == level, "x"])
            npt.assert_array_equal(y, self.df.loc[self.df.a == level, "y"])

        plt.close("all")

        g = ag.PairGrid(self.df, hue="a", hue_order=order)
        g.map_diag(plt.plot)

        for line, level in zip(g.axes[0, 0].lines, order):
            x, y = line.get_xydata().T
            npt.assert_array_equal(x, self.df.loc[self.df.a == level, "x"])
            npt.assert_array_equal(y, self.df.loc[self.df.a == level, "x"])

        plt.close("all")

        g = ag.PairGrid(self.df, hue="a", hue_order=order)
        g.map_lower(plt.plot)

        for line, level in zip(g.axes[1, 0].lines, order):
            x, y = line.get_xydata().T
            npt.assert_array_equal(x, self.df.loc[self.df.a == level, "x"])
            npt.assert_array_equal(y, self.df.loc[self.df.a == level, "y"])

        plt.close("all")

        g = ag.PairGrid(self.df, hue="a", hue_order=order)
        g.map_upper(plt.plot)

        for line, level in zip(g.axes[0, 1].lines, order):
            x, y = line.get_xydata().T
            npt.assert_array_equal(x, self.df.loc[self.df.a == level, "y"])
            npt.assert_array_equal(y, self.df.loc[self.df.a == level, "x"])

        plt.close("all")

    def test_hue_order_missing_level(self):

        order = list("dcaeb")
        g = ag.PairGrid(self.df, hue="a", hue_order=order)
        g.map(plt.plot)

        for line, level in zip(g.axes[1, 0].lines, order):
            x, y = line.get_xydata().T
            npt.assert_array_equal(x, self.df.loc[self.df.a == level, "x"])
            npt.assert_array_equal(y, self.df.loc[self.df.a == level, "y"])

        plt.close("all")

        g = ag.PairGrid(self.df, hue="a", hue_order=order)
        g.map_diag(plt.plot)

        for line, level in zip(g.axes[0, 0].lines, order):
            x, y = line.get_xydata().T
            npt.assert_array_equal(x, self.df.loc[self.df.a == level, "x"])
            npt.assert_array_equal(y, self.df.loc[self.df.a == level, "x"])

        plt.close("all")

        g = ag.PairGrid(self.df, hue="a", hue_order=order)
        g.map_lower(plt.plot)

        for line, level in zip(g.axes[1, 0].lines, order):
            x, y = line.get_xydata().T
            npt.assert_array_equal(x, self.df.loc[self.df.a == level, "x"])
            npt.assert_array_equal(y, self.df.loc[self.df.a == level, "y"])

        plt.close("all")

        g = ag.PairGrid(self.df, hue="a", hue_order=order)
        g.map_upper(plt.plot)

        for line, level in zip(g.axes[0, 1].lines, order):
            x, y = line.get_xydata().T
            npt.assert_array_equal(x, self.df.loc[self.df.a == level, "y"])
            npt.assert_array_equal(y, self.df.loc[self.df.a == level, "x"])

        plt.close("all")

    def test_nondefault_index(self):

        df = self.df.copy().set_index("b")

        plot_vars = ["x", "y", "z"]
        g1 = ag.PairGrid(df)
        g1.map(plt.scatter)

        for i, axes_i in enumerate(g1.axes):
            for j, ax in enumerate(axes_i):
                x_in = self.df[plot_vars[j]]
                y_in = self.df[plot_vars[i]]
                x_out, y_out = ax.collections[0].get_offsets().T
                npt.assert_array_equal(x_in, x_out)
                npt.assert_array_equal(y_in, y_out)

        g2 = ag.PairGrid(df, hue="a")
        g2.map(plt.scatter)

        for i, axes_i in enumerate(g2.axes):
            for j, ax in enumerate(axes_i):
                x_in = self.df[plot_vars[j]]
                y_in = self.df[plot_vars[i]]
                for k, k_level in enumerate(self.df.a.unique()):
                    x_in_k = x_in[self.df.a == k_level]
                    y_in_k = y_in[self.df.a == k_level]
                    x_out, y_out = ax.collections[k].get_offsets().T
                    npt.assert_array_equal(x_in_k, x_out)
                    npt.assert_array_equal(y_in_k, y_out)

    @pytest.mark.parametrize("func", [scatterplot, plt.scatter])
    def test_dropna(self, func):

        df = self.df.copy()
        n_null = 20
        df.loc[np.arange(n_null), "x"] = np.nan

        plot_vars = ["x", "y", "z"]

        g1 = ag.PairGrid(df, vars=plot_vars, dropna=True)
        g1.map(func)

        for i, axes_i in enumerate(g1.axes):
            for j, ax in enumerate(axes_i):
                x_in = df[plot_vars[j]]
                y_in = df[plot_vars[i]]
                x_out, y_out = ax.collections[0].get_offsets().T

                n_valid = (x_in * y_in).notnull().sum()

                assert n_valid == len(x_out)
                assert n_valid == len(y_out)

        g1.map_diag(histplot)
        for i, ax in enumerate(g1.diag_axes):
            var = plot_vars[i]
            count = sum([p.get_height() for p in ax.patches])
            assert count == df[var].notna().sum()

    def test_histplot_legend(self):

        # Tests _extract_legend_handles
        g = ag.PairGrid(self.df, vars=["x", "y"], hue="a")
        g.map_offdiag(histplot)
        g.add_legend()

        assert len(g._legend.legendHandles) == len(self.df["a"].unique())

    def test_pairplot(self):

        vars = ["x", "y", "z"]
        g = ag.pairplot(self.df)

        for ax in g.diag_axes:
            assert len(ax.patches) > 1

        for i, j in zip(*np.triu_indices_from(g.axes, 1)):
            ax = g.axes[i, j]
            x_in = self.df[vars[j]]
            y_in = self.df[vars[i]]
            x_out, y_out = ax.collections[0].get_offsets().T
            npt.assert_array_equal(x_in, x_out)
            npt.assert_array_equal(y_in, y_out)

        for i, j in zip(*np.tril_indices_from(g.axes, -1)):
            ax = g.axes[i, j]
            x_in = self.df[vars[j]]
            y_in = self.df[vars[i]]
            x_out, y_out = ax.collections[0].get_offsets().T
            npt.assert_array_equal(x_in, x_out)
            npt.assert_array_equal(y_in, y_out)

        for i, j in zip(*np.diag_indices_from(g.axes)):
            ax = g.axes[i, j]
            assert len(ax.collections) == 0

        g = ag.pairplot(self.df, hue="a")
        n = len(self.df.a.unique())

        for ax in g.diag_axes:
            assert len(ax.collections) == n

    def test_pairplot_reg(self):

        vars = ["x", "y", "z"]
        g = ag.pairplot(self.df, diag_kind="hist", kind="reg")

        for ax in g.diag_axes:
            assert len(ax.patches)

        for i, j in zip(*np.triu_indices_from(g.axes, 1)):
            ax = g.axes[i, j]
            x_in = self.df[vars[j]]
            y_in = self.df[vars[i]]
            x_out, y_out = ax.collections[0].get_offsets().T
            npt.assert_array_equal(x_in, x_out)
            npt.assert_array_equal(y_in, y_out)

            assert len(ax.lines) == 1
            assert len(ax.collections) == 2

        for i, j in zip(*np.tril_indices_from(g.axes, -1)):
            ax = g.axes[i, j]
            x_in = self.df[vars[j]]
            y_in = self.df[vars[i]]
            x_out, y_out = ax.collections[0].get_offsets().T
            npt.assert_array_equal(x_in, x_out)
            npt.assert_array_equal(y_in, y_out)

            assert len(ax.lines) == 1
            assert len(ax.collections) == 2

        for i, j in zip(*np.diag_indices_from(g.axes)):
            ax = g.axes[i, j]
            assert len(ax.collections) == 0

    def test_pairplot_reg_hue(self):

        markers = ["o", "s", "d"]
        g = ag.pairplot(self.df, kind="reg", hue="a", markers=markers)

        ax = g.axes[-1, 0]
        c1 = ax.collections[0]
        c2 = ax.collections[2]

        assert not np.array_equal(c1.get_facecolor(), c2.get_facecolor())
        assert not np.array_equal(
            c1.get_paths()[0].vertices, c2.get_paths()[0].vertices,
        )

    def test_pairplot_diag_kde(self):

        vars = ["x", "y", "z"]
        g = ag.pairplot(self.df, diag_kind="kde")

        for ax in g.diag_axes:
            assert len(ax.collections) == 1

        for i, j in zip(*np.triu_indices_from(g.axes, 1)):
            ax = g.axes[i, j]
            x_in = self.df[vars[j]]
            y_in = self.df[vars[i]]
            x_out, y_out = ax.collections[0].get_offsets().T
            npt.assert_array_equal(x_in, x_out)
            npt.assert_array_equal(y_in, y_out)

        for i, j in zip(*np.tril_indices_from(g.axes, -1)):
            ax = g.axes[i, j]
            x_in = self.df[vars[j]]
            y_in = self.df[vars[i]]
            x_out, y_out = ax.collections[0].get_offsets().T
            npt.assert_array_equal(x_in, x_out)
            npt.assert_array_equal(y_in, y_out)

        for i, j in zip(*np.diag_indices_from(g.axes)):
            ax = g.axes[i, j]
            assert len(ax.collections) == 0

    def test_pairplot_kde(self):

        f, ax1 = plt.subplots()
        kdeplot(data=self.df, x="x", y="y", ax=ax1)

        g = ag.pairplot(self.df, kind="kde")
        ax2 = g.axes[1, 0]

        assert_plots_equal(ax1, ax2, labels=False)

    def test_pairplot_hist(self):

        f, ax1 = plt.subplots()
        histplot(data=self.df, x="x", y="y", ax=ax1)

        g = ag.pairplot(self.df, kind="hist")
        ax2 = g.axes[1, 0]

        assert_plots_equal(ax1, ax2, labels=False)

    def test_pairplot_markers(self):

        vars = ["x", "y", "z"]
        markers = ["o", "X", "s"]
        g = ag.pairplot(self.df, hue="a", vars=vars, markers=markers)
        m1 = g._legend.legendHandles[0].get_paths()[0]
        m2 = g._legend.legendHandles[1].get_paths()[0]
        assert m1 != m2

        with pytest.raises(ValueError):
            g = ag.pairplot(self.df, hue="a", vars=vars, markers=markers[:-2])

    def test_corner_despine(self):

        g = ag.PairGrid(self.df, corner=True, despine=False)
        g.map_diag(histplot)
        assert g.axes[0, 0].spines["top"].get_visible()

    def test_corner_set(self):

        g = ag.PairGrid(self.df, corner=True, despine=False)
        g.set(xlim=(0, 10))
        assert g.axes[-1, 0].get_xlim() == (0, 10)

    def test_legend(self):

        g1 = ag.pairplot(self.df, hue="a")
        assert isinstance(g1.legend, mpl.legend.Legend)

        g2 = ag.pairplot(self.df)
        assert g2.legend is None


class TestJointGrid:

    rs = np.random.RandomState(sum(map(ord, "JointGrid")))
    x = rs.randn(100)
    y = rs.randn(100)
    x_na = x.copy()
    x_na[10] = np.nan
    x_na[20] = np.nan
    data = pd.DataFrame(dict(x=x, y=y, x_na=x_na))

    def test_margin_grid_from_lists(self):

        g = ag.JointGrid(x=self.x.tolist(), y=self.y.tolist())
        npt.assert_array_equal(g.x, self.x)
        npt.assert_array_equal(g.y, self.y)

    def test_margin_grid_from_arrays(self):

        g = ag.JointGrid(x=self.x, y=self.y)
        npt.assert_array_equal(g.x, self.x)
        npt.assert_array_equal(g.y, self.y)

    def test_margin_grid_from_series(self):

        g = ag.JointGrid(x=self.data.x, y=self.data.y)
        npt.assert_array_equal(g.x, self.x)
        npt.assert_array_equal(g.y, self.y)

    def test_margin_grid_from_dataframe(self):

        g = ag.JointGrid(x="x", y="y", data=self.data)
        npt.assert_array_equal(g.x, self.x)
        npt.assert_array_equal(g.y, self.y)

    def test_margin_grid_from_dataframe_bad_variable(self):

        with pytest.raises(ValueError):
            ag.JointGrid(x="x", y="bad_column", data=self.data)

    def test_margin_grid_axis_labels(self):

        g = ag.JointGrid(x="x", y="y", data=self.data)

        xlabel, ylabel = g.ax_joint.get_xlabel(), g.ax_joint.get_ylabel()
        assert xlabel == "x"
        assert ylabel == "y"

        g.set_axis_labels("x variable", "y variable")
        xlabel, ylabel = g.ax_joint.get_xlabel(), g.ax_joint.get_ylabel()
        assert xlabel == "x variable"
        assert ylabel == "y variable"

    def test_dropna(self):

        g = ag.JointGrid(x="x_na", y="y", data=self.data, dropna=False)
        assert len(g.x) == len(self.x_na)

        g = ag.JointGrid(x="x_na", y="y", data=self.data, dropna=True)
        assert len(g.x) == pd.notnull(self.x_na).sum()

    def test_axlims(self):

        lim = (-3, 3)
        g = ag.JointGrid(x="x", y="y", data=self.data, xlim=lim, ylim=lim)

        assert g.ax_joint.get_xlim() == lim
        assert g.ax_joint.get_ylim() == lim

        assert g.ax_marg_x.get_xlim() == lim
        assert g.ax_marg_y.get_ylim() == lim

    def test_marginal_ticks(self):

        g = ag.JointGrid(marginal_ticks=False)
        assert not sum(t.get_visible() for t in g.ax_marg_x.get_yticklabels())
        assert not sum(t.get_visible() for t in g.ax_marg_y.get_xticklabels())

        g = ag.JointGrid(marginal_ticks=True)
        assert sum(t.get_visible() for t in g.ax_marg_x.get_yticklabels())
        assert sum(t.get_visible() for t in g.ax_marg_y.get_xticklabels())

    def test_bivariate_plot(self):

        g = ag.JointGrid(x="x", y="y", data=self.data)
        g.plot_joint(plt.plot)

        x, y = g.ax_joint.lines[0].get_xydata().T
        npt.assert_array_equal(x, self.x)
        npt.assert_array_equal(y, self.y)

    def test_univariate_plot(self):

        g = ag.JointGrid(x="x", y="x", data=self.data)
        g.plot_marginals(kdeplot)

        _, y1 = g.ax_marg_x.lines[0].get_xydata().T
        y2, _ = g.ax_marg_y.lines[0].get_xydata().T
        npt.assert_array_equal(y1, y2)

    def test_univariate_plot_distplot(self):

        bins = 10
        g = ag.JointGrid(x="x", y="x", data=self.data)
        with pytest.warns(FutureWarning):
            g.plot_marginals(distplot, bins=bins)
        assert len(g.ax_marg_x.patches) == bins
        assert len(g.ax_marg_y.patches) == bins
        for x, y in zip(g.ax_marg_x.patches, g.ax_marg_y.patches):
            assert x.get_height() == y.get_width()

    def test_univariate_plot_matplotlib(self):

        bins = 10
        g = ag.JointGrid(x="x", y="x", data=self.data)
        g.plot_marginals(plt.hist, bins=bins)
        assert len(g.ax_marg_x.patches) == bins
        assert len(g.ax_marg_y.patches) == bins

    def test_plot(self):

        g = ag.JointGrid(x="x", y="x", data=self.data)
        g.plot(plt.plot, kdeplot)

        x, y = g.ax_joint.lines[0].get_xydata().T
        npt.assert_array_equal(x, self.x)
        npt.assert_array_equal(y, self.x)

        _, y1 = g.ax_marg_x.lines[0].get_xydata().T
        y2, _ = g.ax_marg_y.lines[0].get_xydata().T
        npt.assert_array_equal(y1, y2)

    def test_space(self):

        g = ag.JointGrid(x="x", y="y", data=self.data, space=0)

        joint_bounds = g.ax_joint.bbox.bounds
        marg_x_bounds = g.ax_marg_x.bbox.bounds
        marg_y_bounds = g.ax_marg_y.bbox.bounds

        assert joint_bounds[2] == marg_x_bounds[2]
        assert joint_bounds[3] == marg_y_bounds[3]

    @pytest.mark.parametrize(
        "as_vector", [True, False],
    )
    def test_hue(self, long_df, as_vector):

        if as_vector:
            data = None
            x, y, hue = long_df["x"], long_df["y"], long_df["a"]
        else:
            data = long_df
            x, y, hue = "x", "y", "a"

        g = ag.JointGrid(data=data, x=x, y=y, hue=hue)
        g.plot_joint(scatterplot)
        g.plot_marginals(histplot)

        g2 = ag.JointGrid()
        scatterplot(data=long_df, x=x, y=y, hue=hue, ax=g2.ax_joint)
        histplot(data=long_df, x=x, hue=hue, ax=g2.ax_marg_x)
        histplot(data=long_df, y=y, hue=hue, ax=g2.ax_marg_y)

        assert_plots_equal(g.ax_joint, g2.ax_joint)
        assert_plots_equal(g.ax_marg_x, g2.ax_marg_x, labels=False)
        assert_plots_equal(g.ax_marg_y, g2.ax_marg_y, labels=False)

    def test_refline(self):

        g = ag.JointGrid(x="x", y="y", data=self.data)
        g.plot(scatterplot, histplot)
        g.refline()
        assert not g.ax_joint.lines and not g.ax_marg_x.lines and not g.ax_marg_y.lines

        refx = refy = 0.5
        hline = np.array([[0, refy], [1, refy]])
        vline = np.array([[refx, 0], [refx, 1]])
        g.refline(x=refx, y=refy, joint=False, marginal=False)
        assert not g.ax_joint.lines and not g.ax_marg_x.lines and not g.ax_marg_y.lines

        g.refline(x=refx, y=refy)
        assert g.ax_joint.lines[0].get_color() == '.5'
        assert g.ax_joint.lines[0].get_linestyle() == '--'
        assert len(g.ax_joint.lines) == 2
        assert len(g.ax_marg_x.lines) == 1
        assert len(g.ax_marg_y.lines) == 1
        npt.assert_array_equal(g.ax_joint.lines[0].get_xydata(), vline)
        npt.assert_array_equal(g.ax_joint.lines[1].get_xydata(), hline)
        npt.assert_array_equal(g.ax_marg_x.lines[0].get_xydata(), vline)
        npt.assert_array_equal(g.ax_marg_y.lines[0].get_xydata(), hline)

        color, linestyle = 'red', '-'
        g.refline(x=refx, marginal=False, color=color, linestyle=linestyle)
        npt.assert_array_equal(g.ax_joint.lines[-1].get_xydata(), vline)
        assert g.ax_joint.lines[-1].get_color() == color
        assert g.ax_joint.lines[-1].get_linestyle() == linestyle
        assert len(g.ax_marg_x.lines) == len(g.ax_marg_y.lines)

        g.refline(x=refx, joint=False)
        npt.assert_array_equal(g.ax_marg_x.lines[-1].get_xydata(), vline)
        assert len(g.ax_marg_x.lines) == len(g.ax_marg_y.lines) + 1

        g.refline(y=refy, joint=False)
        npt.assert_array_equal(g.ax_marg_y.lines[-1].get_xydata(), hline)
        assert len(g.ax_marg_x.lines) == len(g.ax_marg_y.lines)

        g.refline(y=refy, marginal=False)
        npt.assert_array_equal(g.ax_joint.lines[-1].get_xydata(), hline)
        assert len(g.ax_marg_x.lines) == len(g.ax_marg_y.lines)


class TestJointPlot:

    rs = np.random.RandomState(sum(map(ord, "jointplot")))
    x = rs.randn(100)
    y = rs.randn(100)
    data = pd.DataFrame(dict(x=x, y=y))

    def test_scatter(self):

        g = ag.jointplot(x="x", y="y", data=self.data)
        assert len(g.ax_joint.collections) == 1

        x, y = g.ax_joint.collections[0].get_offsets().T
        assert_array_equal(self.x, x)
        assert_array_equal(self.y, y)

        assert_array_equal(
            [b.get_x() for b in g.ax_marg_x.patches],
            np.histogram_bin_edges(self.x, "auto")[:-1],
        )

        assert_array_equal(
            [b.get_y() for b in g.ax_marg_y.patches],
            np.histogram_bin_edges(self.y, "auto")[:-1],
        )

    def test_scatter_hue(self, long_df):

        g1 = ag.jointplot(data=long_df, x="x", y="y", hue="a")

        g2 = ag.JointGrid()
        scatterplot(data=long_df, x="x", y="y", hue="a", ax=g2.ax_joint)
        kdeplot(data=long_df, x="x", hue="a", ax=g2.ax_marg_x, fill=True)
        kdeplot(data=long_df, y="y", hue="a", ax=g2.ax_marg_y, fill=True)

        assert_plots_equal(g1.ax_joint, g2.ax_joint)
        assert_plots_equal(g1.ax_marg_x, g2.ax_marg_x, labels=False)
        assert_plots_equal(g1.ax_marg_y, g2.ax_marg_y, labels=False)

    def test_reg(self):

        g = ag.jointplot(x="x", y="y", data=self.data, kind="reg")
        assert len(g.ax_joint.collections) == 2

        x, y = g.ax_joint.collections[0].get_offsets().T
        assert_array_equal(self.x, x)
        assert_array_equal(self.y, y)

        assert g.ax_marg_x.patches
        assert g.ax_marg_y.patches

        assert g.ax_marg_x.lines
        assert g.ax_marg_y.lines

    def test_resid(self):

        g = ag.jointplot(x="x", y="y", data=self.data, kind="resid")
        assert g.ax_joint.collections
        assert g.ax_joint.lines
        assert not g.ax_marg_x.lines
        assert not g.ax_marg_y.lines

    def test_hist(self, long_df):

        bins = 3, 6
        g1 = ag.jointplot(data=long_df, x="x", y="y", kind="hist", bins=bins)

        g2 = ag.JointGrid()
        histplot(data=long_df, x="x", y="y", ax=g2.ax_joint, bins=bins)
        histplot(data=long_df, x="x", ax=g2.ax_marg_x, bins=bins[0])
        histplot(data=long_df, y="y", ax=g2.ax_marg_y, bins=bins[1])

        assert_plots_equal(g1.ax_joint, g2.ax_joint)
        assert_plots_equal(g1.ax_marg_x, g2.ax_marg_x, labels=False)
        assert_plots_equal(g1.ax_marg_y, g2.ax_marg_y, labels=False)

    def test_hex(self):

        g = ag.jointplot(x="x", y="y", data=self.data, kind="hex")
        assert g.ax_joint.collections
        assert g.ax_marg_x.patches
        assert g.ax_marg_y.patches

    def test_kde(self, long_df):

        g1 = ag.jointplot(data=long_df, x="x", y="y", kind="kde")

        g2 = ag.JointGrid()
        kdeplot(data=long_df, x="x", y="y", ax=g2.ax_joint)
        kdeplot(data=long_df, x="x", ax=g2.ax_marg_x)
        kdeplot(data=long_df, y="y", ax=g2.ax_marg_y)

        assert_plots_equal(g1.ax_joint, g2.ax_joint)
        assert_plots_equal(g1.ax_marg_x, g2.ax_marg_x, labels=False)
        assert_plots_equal(g1.ax_marg_y, g2.ax_marg_y, labels=False)

    def test_kde_hue(self, long_df):

        g1 = ag.jointplot(data=long_df, x="x", y="y", hue="a", kind="kde")

        g2 = ag.JointGrid()
        kdeplot(data=long_df, x="x", y="y", hue="a", ax=g2.ax_joint)
        kdeplot(data=long_df, x="x", hue="a", ax=g2.ax_marg_x)
        kdeplot(data=long_df, y="y", hue="a", ax=g2.ax_marg_y)

        assert_plots_equal(g1.ax_joint, g2.ax_joint)
        assert_plots_equal(g1.ax_marg_x, g2.ax_marg_x, labels=False)
        assert_plots_equal(g1.ax_marg_y, g2.ax_marg_y, labels=False)

    def test_color(self):

        g = ag.jointplot(x="x", y="y", data=self.data, color="purple")

        purple = mpl.colors.colorConverter.to_rgb("purple")
        scatter_color = g.ax_joint.collections[0].get_facecolor()[0, :3]
        assert tuple(scatter_color) == purple

        hist_color = g.ax_marg_x.patches[0].get_facecolor()[:3]
        assert hist_color == purple

    def test_palette(self, long_df):

        kws = dict(data=long_df, hue="a", palette="Set2")

        g1 = ag.jointplot(x="x", y="y", **kws)

        g2 = ag.JointGrid()
        scatterplot(x="x", y="y", ax=g2.ax_joint, **kws)
        kdeplot(x="x", ax=g2.ax_marg_x, fill=True, **kws)
        kdeplot(y="y", ax=g2.ax_marg_y, fill=True, **kws)

        assert_plots_equal(g1.ax_joint, g2.ax_joint)
        assert_plots_equal(g1.ax_marg_x, g2.ax_marg_x, labels=False)
        assert_plots_equal(g1.ax_marg_y, g2.ax_marg_y, labels=False)

    def test_hex_customise(self):

        # test that default gridsize can be overridden
        g = ag.jointplot(x="x", y="y", data=self.data, kind="hex",
                         joint_kws=dict(gridsize=5))
        assert len(g.ax_joint.collections) == 1
        a = g.ax_joint.collections[0].get_array()
        assert a.shape[0] == 28  # 28 hexagons expected for gridsize 5

    def test_bad_kind(self):

        with pytest.raises(ValueError):
            ag.jointplot(x="x", y="y", data=self.data, kind="not_a_kind")

    def test_unsupported_hue_kind(self):

        for kind in ["reg", "resid", "hex"]:
            with pytest.raises(ValueError):
                ag.jointplot(x="x", y="y", hue="a", data=self.data, kind=kind)

    def test_leaky_dict(self):
        # Validate input dicts are unchanged by jointplot plotting function

        for kwarg in ("joint_kws", "marginal_kws"):
            for kind in ("hex", "kde", "resid", "reg", "scatter"):
                empty_dict = {}
                ag.jointplot(x="x", y="y", data=self.data, kind=kind,
                             **{kwarg: empty_dict})
                assert empty_dict == {}

    def test_distplot_kwarg_warning(self, long_df):

        with pytest.warns(UserWarning):
            g = ag.jointplot(data=long_df, x="x", y="y", marginal_kws=dict(rug=True))
            assert g.ax_marg_x.patches
