import warnings

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib as mpl
import matplotlib.pyplot as plt
from distutils.version import LooseVersion

import nose.tools as nt
import numpy.testing as npt
from numpy.testing.decorators import skipif
try:
    import pandas.testing as tm
except ImportError:
    import pandas.util.testing as tm

from . import PlotTestCase
from .. import axisgrid as ag
from .. import rcmod
from ..palettes import color_palette
from ..distributions import kdeplot, _freedman_diaconis_bins
from ..categorical import pointplot
from ..utils import categorical_order

rs = np.random.RandomState(0)

old_matplotlib = LooseVersion(mpl.__version__) < "1.4"
pandas_has_categoricals = LooseVersion(pd.__version__) >= "0.15"


class TestFacetGrid(PlotTestCase):

    df = pd.DataFrame(dict(x=rs.normal(size=60),
                           y=rs.gamma(4, size=60),
                           a=np.repeat(list("abc"), 20),
                           b=np.tile(list("mn"), 30),
                           c=np.tile(list("tuv"), 20),
                           d=np.tile(list("abcdefghij"), 6)))

    def test_self_data(self):

        g = ag.FacetGrid(self.df)
        nt.assert_is(g.data, self.df)

    def test_self_fig(self):

        g = ag.FacetGrid(self.df)
        nt.assert_is_instance(g.fig, plt.Figure)

    def test_self_axes(self):

        g = ag.FacetGrid(self.df, row="a", col="b", hue="c")
        for ax in g.axes.flat:
            nt.assert_is_instance(ax, plt.Axes)

    def test_axes_array_size(self):

        g1 = ag.FacetGrid(self.df)
        nt.assert_equal(g1.axes.shape, (1, 1))

        g2 = ag.FacetGrid(self.df, row="a")
        nt.assert_equal(g2.axes.shape, (3, 1))

        g3 = ag.FacetGrid(self.df, col="b")
        nt.assert_equal(g3.axes.shape, (1, 2))

        g4 = ag.FacetGrid(self.df, hue="c")
        nt.assert_equal(g4.axes.shape, (1, 1))

        g5 = ag.FacetGrid(self.df, row="a", col="b", hue="c")
        nt.assert_equal(g5.axes.shape, (3, 2))

        for ax in g5.axes.flat:
            nt.assert_is_instance(ax, plt.Axes)

    def test_single_axes(self):

        g1 = ag.FacetGrid(self.df)
        nt.assert_is_instance(g1.ax, plt.Axes)

        g2 = ag.FacetGrid(self.df, row="a")
        with nt.assert_raises(AttributeError):
            g2.ax

        g3 = ag.FacetGrid(self.df, col="a")
        with nt.assert_raises(AttributeError):
            g3.ax

        g4 = ag.FacetGrid(self.df, col="a", row="b")
        with nt.assert_raises(AttributeError):
            g4.ax

    def test_col_wrap(self):

        g = ag.FacetGrid(self.df, col="d")
        nt.assert_equal(g.axes.shape, (1, 10))
        nt.assert_is(g.facet_axis(0, 8), g.axes[0, 8])

        g_wrap = ag.FacetGrid(self.df, col="d", col_wrap=4)
        nt.assert_equal(g_wrap.axes.shape, (10,))
        nt.assert_is(g_wrap.facet_axis(0, 8), g_wrap.axes[8])
        nt.assert_equal(g_wrap._ncol, 4)
        nt.assert_equal(g_wrap._nrow, 3)

        with nt.assert_raises(ValueError):
            g = ag.FacetGrid(self.df, row="b", col="d", col_wrap=4)

        df = self.df.copy()
        df.loc[df.d == "j"] = np.nan
        g_missing = ag.FacetGrid(df, col="d")
        nt.assert_equal(g_missing.axes.shape, (1, 9))

        g_missing_wrap = ag.FacetGrid(df, col="d", col_wrap=4)
        nt.assert_equal(g_missing_wrap.axes.shape, (9,))

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

    def test_figure_size(self):

        g = ag.FacetGrid(self.df, row="a", col="b")
        npt.assert_array_equal(g.fig.get_size_inches(), (6, 9))

        g = ag.FacetGrid(self.df, row="a", col="b", size=6)
        npt.assert_array_equal(g.fig.get_size_inches(), (12, 18))

        g = ag.FacetGrid(self.df, col="c", size=4, aspect=.5)
        npt.assert_array_equal(g.fig.get_size_inches(), (6, 4))

    def test_figure_size_with_legend(self):

        g1 = ag.FacetGrid(self.df, col="a", hue="c", size=4, aspect=.5)
        npt.assert_array_equal(g1.fig.get_size_inches(), (6, 4))
        g1.add_legend()
        nt.assert_greater(g1.fig.get_size_inches()[0], 6)

        g2 = ag.FacetGrid(self.df, col="a", hue="c", size=4, aspect=.5,
                          legend_out=False)
        npt.assert_array_equal(g2.fig.get_size_inches(), (6, 4))
        g2.add_legend()
        npt.assert_array_equal(g2.fig.get_size_inches(), (6, 4))

    def test_legend_data(self):

        g1 = ag.FacetGrid(self.df, hue="a")
        g1.map(plt.plot, "x", "y")
        g1.add_legend()
        palette = color_palette(n_colors=3)

        nt.assert_equal(g1._legend.get_title().get_text(), "a")

        a_levels = sorted(self.df.a.unique())

        lines = g1._legend.get_lines()
        nt.assert_equal(len(lines), len(a_levels))

        for line, hue in zip(lines, palette):
            nt.assert_equal(line.get_color(), hue)

        labels = g1._legend.get_texts()
        nt.assert_equal(len(labels), len(a_levels))

        for label, level in zip(labels, a_levels):
            nt.assert_equal(label.get_text(), level)

    def test_legend_data_missing_level(self):

        g1 = ag.FacetGrid(self.df, hue="a", hue_order=list("azbc"))
        g1.map(plt.plot, "x", "y")
        g1.add_legend()

        b, g, r, p = color_palette(n_colors=4)
        palette = [b, r, p]

        nt.assert_equal(g1._legend.get_title().get_text(), "a")

        a_levels = sorted(self.df.a.unique())

        lines = g1._legend.get_lines()
        nt.assert_equal(len(lines), len(a_levels))

        for line, hue in zip(lines, palette):
            nt.assert_equal(line.get_color(), hue)

        labels = g1._legend.get_texts()
        nt.assert_equal(len(labels), 4)

        for label, level in zip(labels, list("azbc")):
            nt.assert_equal(label.get_text(), level)

    def test_get_boolean_legend_data(self):

        self.df["b_bool"] = self.df.b == "m"
        g1 = ag.FacetGrid(self.df, hue="b_bool")
        g1.map(plt.plot, "x", "y")
        g1.add_legend()
        palette = color_palette(n_colors=2)

        nt.assert_equal(g1._legend.get_title().get_text(), "b_bool")

        b_levels = list(map(str, categorical_order(self.df.b_bool)))

        lines = g1._legend.get_lines()
        nt.assert_equal(len(lines), len(b_levels))

        for line, hue in zip(lines, palette):
            nt.assert_equal(line.get_color(), hue)

        labels = g1._legend.get_texts()
        nt.assert_equal(len(labels), len(b_levels))

        for label, level in zip(labels, b_levels):
            nt.assert_equal(label.get_text(), level)

    def test_legend_options(self):

        g1 = ag.FacetGrid(self.df, hue="b")
        g1.map(plt.plot, "x", "y")
        g1.add_legend()

    def test_legendout_with_colwrap(self):

        g = ag.FacetGrid(self.df, col="d", hue='b',
                         col_wrap=4, legend_out=False)
        g.map(plt.plot, "x", "y", linewidth=3)
        g.add_legend()

    def test_subplot_kws(self):

        g = ag.FacetGrid(self.df, despine=False,
                         subplot_kws=dict(projection="polar"))
        for ax in g.axes.flat:
            nt.assert_true("PolarAxesSubplot" in str(type(ax)))

    @skipif(old_matplotlib)
    def test_gridspec_kws(self):
        ratios = [3, 1, 2]
        sizes = [0.46, 0.15, 0.31]

        gskws = dict(width_ratios=ratios, height_ratios=ratios)
        g = ag.FacetGrid(self.df, col='c', row='a', gridspec_kws=gskws)

        # clear out all ticks
        for ax in g.axes.flat:
            ax.set_xticks([])
            ax.set_yticks([])

        g.fig.tight_layout()
        widths, heights = np.meshgrid(sizes, sizes)
        for n, ax in enumerate(g.axes.flat):
            npt.assert_almost_equal(
                ax.get_position().width,
                widths.flatten()[n],
                decimal=2
            )
            npt.assert_almost_equal(
                ax.get_position().height,
                heights.flatten()[n],
                decimal=2
            )

    @skipif(old_matplotlib)
    def test_gridspec_kws_col_wrap(self):
        ratios = [3, 1, 2, 1, 1]
        sizes = [0.46, 0.15, 0.31]

        gskws = dict(width_ratios=ratios)
        with warnings.catch_warnings():
            warnings.resetwarnings()
            warnings.simplefilter("always")
            npt.assert_warns(UserWarning, ag.FacetGrid, self.df, col='d',
                             col_wrap=5, gridspec_kws=gskws)

    @skipif(not old_matplotlib)
    def test_gridsic_kws_old_mpl(self):
        ratios = [3, 1, 2]
        sizes = [0.46, 0.15, 0.31]

        gskws = dict(width_ratios=ratios, height_ratios=ratios)
        with warnings.catch_warnings():
            warnings.resetwarnings()
            warnings.simplefilter("always")
            npt.assert_warns(UserWarning, ag.FacetGrid, self.df, col='c',
                             row='a', gridspec_kws=gskws)

    def test_data_generator(self):

        g = ag.FacetGrid(self.df, row="a")
        d = list(g.facet_data())
        nt.assert_equal(len(d), 3)

        tup, data = d[0]
        nt.assert_equal(tup, (0, 0, 0))
        nt.assert_true((data["a"] == "a").all())

        tup, data = d[1]
        nt.assert_equal(tup, (1, 0, 0))
        nt.assert_true((data["a"] == "b").all())

        g = ag.FacetGrid(self.df, row="a", col="b")
        d = list(g.facet_data())
        nt.assert_equal(len(d), 6)

        tup, data = d[0]
        nt.assert_equal(tup, (0, 0, 0))
        nt.assert_true((data["a"] == "a").all())
        nt.assert_true((data["b"] == "m").all())

        tup, data = d[1]
        nt.assert_equal(tup, (0, 1, 0))
        nt.assert_true((data["a"] == "a").all())
        nt.assert_true((data["b"] == "n").all())

        tup, data = d[2]
        nt.assert_equal(tup, (1, 0, 0))
        nt.assert_true((data["a"] == "b").all())
        nt.assert_true((data["b"] == "m").all())

        g = ag.FacetGrid(self.df, hue="c")
        d = list(g.facet_data())
        nt.assert_equal(len(d), 3)
        tup, data = d[1]
        nt.assert_equal(tup, (0, 0, 1))
        nt.assert_true((data["c"] == "u").all())

    def test_map(self):

        g = ag.FacetGrid(self.df, row="a", col="b", hue="c")
        g.map(plt.plot, "x", "y", linewidth=3)

        lines = g.axes[0, 0].lines
        nt.assert_equal(len(lines), 3)

        line1, _, _ = lines
        nt.assert_equal(line1.get_linewidth(), 3)
        x, y = line1.get_data()
        mask = (self.df.a == "a") & (self.df.b == "m") & (self.df.c == "t")
        npt.assert_array_equal(x, self.df.x[mask])
        npt.assert_array_equal(y, self.df.y[mask])

    def test_map_dataframe(self):

        g = ag.FacetGrid(self.df, row="a", col="b", hue="c")
        plot = lambda x, y, data=None, **kws: plt.plot(data[x], data[y], **kws)
        g.map_dataframe(plot, "x", "y", linestyle="--")

        lines = g.axes[0, 0].lines
        nt.assert_equal(len(lines), 3)

        line1, _, _ = lines
        nt.assert_equal(line1.get_linestyle(), "--")
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
        nt.assert_equal(g.axes[0, 0].get_title(), "a = a | b = m")
        nt.assert_equal(g.axes[0, 1].get_title(), "a = a | b = n")
        nt.assert_equal(g.axes[1, 0].get_title(), "a = b | b = m")

        # Test a provided title
        g.set_titles("{row_var} == {row_name} \/ {col_var} == {col_name}")
        nt.assert_equal(g.axes[0, 0].get_title(), "a == a \/ b == m")
        nt.assert_equal(g.axes[0, 1].get_title(), "a == a \/ b == n")
        nt.assert_equal(g.axes[1, 0].get_title(), "a == b \/ b == m")

        # Test a single row
        g = ag.FacetGrid(self.df,  col="b")
        g.map(plt.plot, "x", "y")

        # Test the default titles
        nt.assert_equal(g.axes[0, 0].get_title(), "b = m")
        nt.assert_equal(g.axes[0, 1].get_title(), "b = n")

        # test with dropna=False
        g = ag.FacetGrid(self.df, col="b", hue="b", dropna=False)
        g.map(plt.plot, 'x', 'y')

    def test_set_titles_margin_titles(self):

        g = ag.FacetGrid(self.df, row="a", col="b", margin_titles=True)
        g.map(plt.plot, "x", "y")

        # Test the default titles
        nt.assert_equal(g.axes[0, 0].get_title(), "b = m")
        nt.assert_equal(g.axes[0, 1].get_title(), "b = n")
        nt.assert_equal(g.axes[1, 0].get_title(), "")

        # Test the row "titles"
        nt.assert_equal(g.axes[0, 1].texts[0].get_text(), "a = a")
        nt.assert_equal(g.axes[1, 1].texts[0].get_text(), "a = b")

        # Test a provided title
        g.set_titles(col_template="{col_var} == {col_name}")
        nt.assert_equal(g.axes[0, 0].get_title(), "b == m")
        nt.assert_equal(g.axes[0, 1].get_title(), "b == n")
        nt.assert_equal(g.axes[1, 0].get_title(), "")

    def test_set_ticklabels(self):

        g = ag.FacetGrid(self.df, row="a", col="b")
        g.map(plt.plot, "x", "y")
        xlab = [l.get_text() + "h" for l in g.axes[1, 0].get_xticklabels()]
        ylab = [l.get_text() for l in g.axes[1, 0].get_yticklabels()]

        g.set_xticklabels(xlab)
        g.set_yticklabels(rotation=90)

        got_x = [l.get_text() + "h" for l in g.axes[1, 1].get_xticklabels()]
        got_y = [l.get_text() for l in g.axes[0, 0].get_yticklabels()]
        npt.assert_array_equal(got_x, xlab)
        npt.assert_array_equal(got_y, ylab)

        x, y = np.arange(10), np.arange(10)
        df = pd.DataFrame(np.c_[x, y], columns=["x", "y"])
        g = ag.FacetGrid(df).map(pointplot, "x", "y", order=x)
        g.set_xticklabels(step=2)
        got_x = [int(l.get_text()) for l in g.axes[0, 0].get_xticklabels()]
        npt.assert_array_equal(x[::2], got_x)

        g = ag.FacetGrid(self.df, col="d", col_wrap=5)
        g.map(plt.plot, "x", "y")
        g.set_xticklabels(rotation=45)
        g.set_yticklabels(rotation=75)
        for ax in g._bottom_axes:
            for l in ax.get_xticklabels():
                nt.assert_equal(l.get_rotation(), 45)
        for ax in g._left_axes:
            for l in ax.get_yticklabels():
                nt.assert_equal(l.get_rotation(), 75)

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

    def test_axis_lims(self):

        g = ag.FacetGrid(self.df, row="a", col="b", xlim=(0, 4), ylim=(-2, 3))
        nt.assert_equal(g.axes[0, 0].get_xlim(), (0, 4))
        nt.assert_equal(g.axes[0, 0].get_ylim(), (-2, 3))

    def test_data_orders(self):

        g = ag.FacetGrid(self.df, row="a", col="b", hue="c")

        nt.assert_equal(g.row_names, list("abc"))
        nt.assert_equal(g.col_names, list("mn"))
        nt.assert_equal(g.hue_names, list("tuv"))
        nt.assert_equal(g.axes.shape, (3, 2))

        g = ag.FacetGrid(self.df, row="a", col="b", hue="c",
                         row_order=list("bca"),
                         col_order=list("nm"),
                         hue_order=list("vtu"))

        nt.assert_equal(g.row_names, list("bca"))
        nt.assert_equal(g.col_names, list("nm"))
        nt.assert_equal(g.hue_names, list("vtu"))
        nt.assert_equal(g.axes.shape, (3, 2))

        g = ag.FacetGrid(self.df, row="a", col="b", hue="c",
                         row_order=list("bcda"),
                         col_order=list("nom"),
                         hue_order=list("qvtu"))

        nt.assert_equal(g.row_names, list("bcda"))
        nt.assert_equal(g.col_names, list("nom"))
        nt.assert_equal(g.hue_names, list("qvtu"))
        nt.assert_equal(g.axes.shape, (4, 3))

    def test_palette(self):

        rcmod.set()

        g = ag.FacetGrid(self.df, hue="c")
        nt.assert_equal(g._colors, color_palette(n_colors=3))

        g = ag.FacetGrid(self.df, hue="d")
        nt.assert_equal(g._colors, color_palette("husl", 10))

        g = ag.FacetGrid(self.df, hue="c", palette="Set2")
        nt.assert_equal(g._colors, color_palette("Set2", 3))

        dict_pal = dict(t="red", u="green", v="blue")
        list_pal = color_palette(["red", "green", "blue"], 3)
        g = ag.FacetGrid(self.df, hue="c", palette=dict_pal)
        nt.assert_equal(g._colors, list_pal)

        list_pal = color_palette(["green", "blue", "red"], 3)
        g = ag.FacetGrid(self.df, hue="c", hue_order=list("uvt"),
                         palette=dict_pal)
        nt.assert_equal(g._colors, list_pal)

    def test_hue_kws(self):

        kws = dict(marker=["o", "s", "D"])
        g = ag.FacetGrid(self.df, hue="c", hue_kws=kws)
        g.map(plt.plot, "x", "y")

        for line, marker in zip(g.axes[0, 0].lines, kws["marker"]):
            nt.assert_equal(line.get_marker(), marker)

    def test_dropna(self):

        df = self.df.copy()
        hasna = pd.Series(np.tile(np.arange(6), 10), dtype=np.float)
        hasna[hasna == 5] = np.nan
        df["hasna"] = hasna
        g = ag.FacetGrid(df, dropna=False, row="hasna")
        nt.assert_equal(g._not_na.sum(), 60)

        g = ag.FacetGrid(df, dropna=True, row="hasna")
        nt.assert_equal(g._not_na.sum(), 50)

    def test_unicode_column_label_with_rows(self):

        # use a smaller copy of the default testing data frame:
        df = self.df.copy()
        df = df[["a", "b", "x"]]

        # rename column 'a' (which will be used for the columns in the grid)
        # by using a Unicode string:
        unicode_column_label = u"\u01ff\u02ff\u03ff"
        df = df.rename(columns={"a": unicode_column_label})

        # ensure that the data frame columns have the expected names:
        nt.assert_equal(list(df.columns), [unicode_column_label, "b", "x"])

        # plot the grid -- if successful, no UnicodeEncodingError should
        # occur:
        g = ag.FacetGrid(df, col=unicode_column_label, row="b")
        g = g.map(plt.plot, "x")

    def test_unicode_column_label_no_rows(self):

        # use a smaller copy of the default testing data frame:
        df = self.df.copy()
        df = df[["a", "x"]]

        # rename column 'a' (which will be used for the columns in the grid)
        # by using a Unicode string:
        unicode_column_label = u"\u01ff\u02ff\u03ff"
        df = df.rename(columns={"a": unicode_column_label})

        # ensure that the data frame columns have the expected names:
        nt.assert_equal(list(df.columns), [unicode_column_label, "x"])

        # plot the grid -- if successful, no UnicodeEncodingError should
        # occur:
        g = ag.FacetGrid(df, col=unicode_column_label)
        g = g.map(plt.plot, "x")

    def test_unicode_row_label_with_columns(self):

        # use a smaller copy of the default testing data frame:
        df = self.df.copy()
        df = df[["a", "b", "x"]]

        # rename column 'b' (which will be used for the rows in the grid)
        # by using a Unicode string:
        unicode_row_label = u"\u01ff\u02ff\u03ff"
        df = df.rename(columns={"b": unicode_row_label})

        # ensure that the data frame columns have the expected names:
        nt.assert_equal(list(df.columns), ["a", unicode_row_label, "x"])

        # plot the grid -- if successful, no UnicodeEncodingError should
        # occur:
        g = ag.FacetGrid(df, col="a", row=unicode_row_label)
        g = g.map(plt.plot, "x")

    def test_unicode_row_label_no_columns(self):

        # use a smaller copy of the default testing data frame:
        df = self.df.copy()
        df = df[["b", "x"]]

        # rename column 'b' (which will be used for the rows in the grid)
        # by using a Unicode string:
        unicode_row_label = u"\u01ff\u02ff\u03ff"
        df = df.rename(columns={"b": unicode_row_label})

        # ensure that the data frame columns have the expected names:
        nt.assert_equal(list(df.columns), [unicode_row_label, "x"])

        # plot the grid -- if successful, no UnicodeEncodingError should
        # occur:
        g = ag.FacetGrid(df, row=unicode_row_label)
        g = g.map(plt.plot, "x")

    def test_unicode_content_with_row_and_column(self):

        df = self.df.copy()

        # replace content of column 'a' (which will form the columns in the
        # grid) by Unicode characters:
        unicode_column_val = np.repeat((u'\u01ff', u'\u02ff', u'\u03ff'), 20)
        df["a"] = unicode_column_val

        # make sure that the replacement worked as expected:
        nt.assert_equal(
            list(df["a"]),
            [u'\u01ff'] * 20 + [u'\u02ff'] * 20 + [u'\u03ff'] * 20)

        # plot the grid -- if successful, no UnicodeEncodingError should
        # occur:
        g = ag.FacetGrid(df, col="a", row="b")
        g = g.map(plt.plot, "x")

    def test_unicode_content_no_rows(self):

        df = self.df.copy()

        # replace content of column 'a' (which will form the columns in the
        # grid) by Unicode characters:
        unicode_column_val = np.repeat((u'\u01ff', u'\u02ff', u'\u03ff'), 20)
        df["a"] = unicode_column_val

        # make sure that the replacement worked as expected:
        nt.assert_equal(
            list(df["a"]),
            [u'\u01ff'] * 20 + [u'\u02ff'] * 20 + [u'\u03ff'] * 20)

        # plot the grid -- if successful, no UnicodeEncodingError should
        # occur:
        g = ag.FacetGrid(df, col="a")
        g = g.map(plt.plot, "x")

    def test_unicode_content_no_columns(self):

        df = self.df.copy()

        # replace content of column 'a' (which will form the rows in the
        # grid) by Unicode characters:
        unicode_column_val = np.repeat((u'\u01ff', u'\u02ff', u'\u03ff'), 20)
        df["b"] = unicode_column_val

        # make sure that the replacement worked as expected:
        nt.assert_equal(
            list(df["b"]),
            [u'\u01ff'] * 20 + [u'\u02ff'] * 20 + [u'\u03ff'] * 20)

        # plot the grid -- if successful, no UnicodeEncodingError should
        # occur:
        g = ag.FacetGrid(df, row="b")
        g = g.map(plt.plot, "x")

    @skipif(not pandas_has_categoricals)
    def test_categorical_column_missing_categories(self):

        df = self.df.copy()
        df['a'] = df['a'].astype('category')

        g = ag.FacetGrid(df[df['a'] == 'a'], col="a", col_wrap=1)

        nt.assert_equal(g.axes.shape, (len(df['a'].cat.categories),))

    def test_categorical_warning(self):

        g = ag.FacetGrid(self.df, col="b")
        with warnings.catch_warnings():
            warnings.resetwarnings()
            warnings.simplefilter("always")
            npt.assert_warns(UserWarning, g.map, pointplot, "b", "x")


class TestPairGrid(PlotTestCase):

    rs = np.random.RandomState(sum(map(ord, "PairGrid")))
    df = pd.DataFrame(dict(x=rs.normal(size=80),
                           y=rs.randint(0, 4, size=(80)),
                           z=rs.gamma(3, size=80),
                           a=np.repeat(list("abcd"), 20),
                           b=np.repeat(list("abcdefgh"), 10)))

    def test_self_data(self):

        g = ag.PairGrid(self.df)
        nt.assert_is(g.data, self.df)

    def test_ignore_datelike_data(self):

        df = self.df.copy()
        df['date'] = pd.date_range('2010-01-01', periods=len(df), freq='d')
        result = ag.PairGrid(self.df).data
        expected = df.drop('date', axis=1)
        tm.assert_frame_equal(result, expected)

    def test_self_fig(self):

        g = ag.PairGrid(self.df)
        nt.assert_is_instance(g.fig, plt.Figure)

    def test_self_axes(self):

        g = ag.PairGrid(self.df)
        for ax in g.axes.flat:
            nt.assert_is_instance(ax, plt.Axes)

    def test_default_axes(self):

        g = ag.PairGrid(self.df)
        nt.assert_equal(g.axes.shape, (3, 3))
        nt.assert_equal(g.x_vars, ["x", "y", "z"])
        nt.assert_equal(g.y_vars, ["x", "y", "z"])
        nt.assert_true(g.square_grid)

    def test_specific_square_axes(self):

        vars = ["z", "x"]
        g = ag.PairGrid(self.df, vars=vars)
        nt.assert_equal(g.axes.shape, (len(vars), len(vars)))
        nt.assert_equal(g.x_vars, vars)
        nt.assert_equal(g.y_vars, vars)
        nt.assert_true(g.square_grid)

    def test_specific_nonsquare_axes(self):

        x_vars = ["x", "y"]
        y_vars = ["z", "y", "x"]
        g = ag.PairGrid(self.df, x_vars=x_vars, y_vars=y_vars)
        nt.assert_equal(g.axes.shape, (len(y_vars), len(x_vars)))
        nt.assert_equal(g.x_vars, x_vars)
        nt.assert_equal(g.y_vars, y_vars)
        nt.assert_true(not g.square_grid)

        x_vars = ["x", "y"]
        y_vars = "z"
        g = ag.PairGrid(self.df, x_vars=x_vars, y_vars=y_vars)
        nt.assert_equal(g.axes.shape, (len(y_vars), len(x_vars)))
        nt.assert_equal(g.x_vars, list(x_vars))
        nt.assert_equal(g.y_vars, list(y_vars))
        nt.assert_true(not g.square_grid)

    def test_specific_square_axes_with_array(self):

        vars = np.array(["z", "x"])
        g = ag.PairGrid(self.df, vars=vars)
        nt.assert_equal(g.axes.shape, (len(vars), len(vars)))
        nt.assert_equal(g.x_vars, list(vars))
        nt.assert_equal(g.y_vars, list(vars))
        nt.assert_true(g.square_grid)

    def test_specific_nonsquare_axes_with_array(self):

        x_vars = np.array(["x", "y"])
        y_vars = np.array(["z", "y", "x"])
        g = ag.PairGrid(self.df, x_vars=x_vars, y_vars=y_vars)
        nt.assert_equal(g.axes.shape, (len(y_vars), len(x_vars)))
        nt.assert_equal(g.x_vars, list(x_vars))
        nt.assert_equal(g.y_vars, list(y_vars))
        nt.assert_true(not g.square_grid)

    def test_size(self):

        g1 = ag.PairGrid(self.df, size=3)
        npt.assert_array_equal(g1.fig.get_size_inches(), (9, 9))

        g2 = ag.PairGrid(self.df, size=4, aspect=.5)
        npt.assert_array_equal(g2.fig.get_size_inches(), (6, 12))

        g3 = ag.PairGrid(self.df, y_vars=["z"], x_vars=["x", "y"],
                         size=2, aspect=2)
        npt.assert_array_equal(g3.fig.get_size_inches(), (8, 2))

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

        g2 = ag.PairGrid(self.df, "a")
        g2.map(plt.scatter)

        for i, axes_i in enumerate(g2.axes):
            for j, ax in enumerate(axes_i):
                x_in = self.df[vars[j]]
                y_in = self.df[vars[i]]
                for k, k_level in enumerate("abcd"):
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
            nt.assert_equal(len(ax.collections), 0)

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
            nt.assert_equal(len(ax.collections), 0)

    @skipif(old_matplotlib)
    def test_map_diag(self):

        g1 = ag.PairGrid(self.df)
        g1.map_diag(plt.hist)

        for ax in g1.diag_axes:
            nt.assert_equal(len(ax.patches), 10)

        g2 = ag.PairGrid(self.df)
        g2.map_diag(plt.hist, bins=15)

        for ax in g2.diag_axes:
            nt.assert_equal(len(ax.patches), 15)

        g3 = ag.PairGrid(self.df, hue="a")
        g3.map_diag(plt.hist)

        for ax in g3.diag_axes:
            nt.assert_equal(len(ax.patches), 40)

        g4 = ag.PairGrid(self.df, hue="a")
        g4.map_diag(plt.hist, histtype='step')

        for ax in g4.diag_axes:
            for ptch in ax.patches:
                nt.assert_equal(ptch.fill, False)

    @skipif(old_matplotlib)
    def test_map_diag_color(self):

        color = "red"
        rgb_color = mpl.colors.colorConverter.to_rgba(color)

        g1 = ag.PairGrid(self.df)
        g1.map_diag(plt.hist, color=color)

        for ax in g1.diag_axes:
            for patch in ax.patches:
                nt.assert_equals(patch.get_facecolor(), rgb_color)

        g2 = ag.PairGrid(self.df)
        g2.map_diag(kdeplot, color='red')

        for ax in g2.diag_axes:
            for line in ax.lines:
                nt.assert_equals(line.get_color(), color)

    @skipif(old_matplotlib)
    def test_map_diag_palette(self):

        pal = color_palette(n_colors=len(self.df.a.unique()))
        g = ag.PairGrid(self.df, hue="a")
        g.map_diag(kdeplot)

        for ax in g.diag_axes:
            for line, color in zip(ax.lines, pal):
                nt.assert_equals(line.get_color(), color)

    @skipif(old_matplotlib)
    def test_map_diag_and_offdiag(self):

        vars = ["x", "y", "z"]
        g = ag.PairGrid(self.df)
        g.map_offdiag(plt.scatter)
        g.map_diag(plt.hist)

        for ax in g.diag_axes:
            nt.assert_equal(len(ax.patches), 10)

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
            nt.assert_equal(len(ax.collections), 0)

    def test_palette(self):

        rcmod.set()

        g = ag.PairGrid(self.df, hue="a")
        nt.assert_equal(g.palette, color_palette(n_colors=4))

        g = ag.PairGrid(self.df, hue="b")
        nt.assert_equal(g.palette, color_palette("husl", 8))

        g = ag.PairGrid(self.df, hue="a", palette="Set2")
        nt.assert_equal(g.palette, color_palette("Set2", 4))

        dict_pal = dict(a="red", b="green", c="blue", d="purple")
        list_pal = color_palette(["red", "green", "blue", "purple"], 4)
        g = ag.PairGrid(self.df, hue="a", palette=dict_pal)
        nt.assert_equal(g.palette, list_pal)

        list_pal = color_palette(["purple", "blue", "red", "green"], 4)
        g = ag.PairGrid(self.df, hue="a", hue_order=list("dcab"),
                        palette=dict_pal)
        nt.assert_equal(g.palette, list_pal)

    def test_hue_kws(self):

        kws = dict(marker=["o", "s", "d", "+"])
        g = ag.PairGrid(self.df, hue="a", hue_kws=kws)
        g.map(plt.plot)

        for line, marker in zip(g.axes[0, 0].lines, kws["marker"]):
            nt.assert_equal(line.get_marker(), marker)

        g = ag.PairGrid(self.df, hue="a", hue_kws=kws,
                        hue_order=list("dcab"))
        g.map(plt.plot)

        for line, marker in zip(g.axes[0, 0].lines, kws["marker"]):
            nt.assert_equal(line.get_marker(), marker)

    @skipif(old_matplotlib)
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

    @skipif(old_matplotlib)
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

        vars = ["x", "y", "z"]
        g1 = ag.PairGrid(df)
        g1.map(plt.scatter)

        for i, axes_i in enumerate(g1.axes):
            for j, ax in enumerate(axes_i):
                x_in = self.df[vars[j]]
                y_in = self.df[vars[i]]
                x_out, y_out = ax.collections[0].get_offsets().T
                npt.assert_array_equal(x_in, x_out)
                npt.assert_array_equal(y_in, y_out)

        g2 = ag.PairGrid(df, "a")
        g2.map(plt.scatter)

        for i, axes_i in enumerate(g2.axes):
            for j, ax in enumerate(axes_i):
                x_in = self.df[vars[j]]
                y_in = self.df[vars[i]]
                for k, k_level in enumerate("abcd"):
                    x_in_k = x_in[self.df.a == k_level]
                    y_in_k = y_in[self.df.a == k_level]
                    x_out, y_out = ax.collections[k].get_offsets().T
                npt.assert_array_equal(x_in_k, x_out)
                npt.assert_array_equal(y_in_k, y_out)

    @skipif(old_matplotlib)
    def test_pairplot(self):

        vars = ["x", "y", "z"]
        g = ag.pairplot(self.df)

        for ax in g.diag_axes:
            nt.assert_equal(len(ax.patches), 10)

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
            nt.assert_equal(len(ax.collections), 0)

    @skipif(old_matplotlib)
    def test_pairplot_reg(self):

        vars = ["x", "y", "z"]
        g = ag.pairplot(self.df, kind="reg")

        for ax in g.diag_axes:
            nt.assert_equal(len(ax.patches), 10)

        for i, j in zip(*np.triu_indices_from(g.axes, 1)):
            ax = g.axes[i, j]
            x_in = self.df[vars[j]]
            y_in = self.df[vars[i]]
            x_out, y_out = ax.collections[0].get_offsets().T
            npt.assert_array_equal(x_in, x_out)
            npt.assert_array_equal(y_in, y_out)

            nt.assert_equal(len(ax.lines), 1)
            nt.assert_equal(len(ax.collections), 2)

        for i, j in zip(*np.tril_indices_from(g.axes, -1)):
            ax = g.axes[i, j]
            x_in = self.df[vars[j]]
            y_in = self.df[vars[i]]
            x_out, y_out = ax.collections[0].get_offsets().T
            npt.assert_array_equal(x_in, x_out)
            npt.assert_array_equal(y_in, y_out)

            nt.assert_equal(len(ax.lines), 1)
            nt.assert_equal(len(ax.collections), 2)

        for i, j in zip(*np.diag_indices_from(g.axes)):
            ax = g.axes[i, j]
            nt.assert_equal(len(ax.collections), 0)

    @skipif(old_matplotlib)
    def test_pairplot_kde(self):

        vars = ["x", "y", "z"]
        g = ag.pairplot(self.df, diag_kind="kde")

        for ax in g.diag_axes:
            nt.assert_equal(len(ax.lines), 1)

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
            nt.assert_equal(len(ax.collections), 0)

    @skipif(old_matplotlib)
    def test_pairplot_markers(self):

        vars = ["x", "y", "z"]
        markers = ["o", "x", "s", "d"]
        g = ag.pairplot(self.df, hue="a", vars=vars, markers=markers)
        nt.assert_equal(g.hue_kws["marker"], markers)
        plt.close("all")

        with nt.assert_raises(ValueError):
            g = ag.pairplot(self.df, hue="a", vars=vars, markers=markers[:-2])


class TestJointGrid(PlotTestCase):

    rs = np.random.RandomState(sum(map(ord, "JointGrid")))
    x = rs.randn(100)
    y = rs.randn(100)
    x_na = x.copy()
    x_na[10] = np.nan
    x_na[20] = np.nan
    data = pd.DataFrame(dict(x=x, y=y, x_na=x_na))

    def test_margin_grid_from_arrays(self):

        g = ag.JointGrid(self.x, self.y)
        npt.assert_array_equal(g.x, self.x)
        npt.assert_array_equal(g.y, self.y)

    def test_margin_grid_from_series(self):

        g = ag.JointGrid(self.data.x, self.data.y)
        npt.assert_array_equal(g.x, self.x)
        npt.assert_array_equal(g.y, self.y)

    def test_margin_grid_from_dataframe(self):

        g = ag.JointGrid("x", "y", self.data)
        npt.assert_array_equal(g.x, self.x)
        npt.assert_array_equal(g.y, self.y)

    def test_margin_grid_from_dataframe_bad_variable(self):

        with nt.assert_raises(ValueError):
            g = ag.JointGrid("x", "bad_column", self.data)

    def test_margin_grid_axis_labels(self):

        g = ag.JointGrid("x", "y", self.data)

        xlabel, ylabel = g.ax_joint.get_xlabel(), g.ax_joint.get_ylabel()
        nt.assert_equal(xlabel, "x")
        nt.assert_equal(ylabel, "y")

        g.set_axis_labels("x variable", "y variable")
        xlabel, ylabel = g.ax_joint.get_xlabel(), g.ax_joint.get_ylabel()
        nt.assert_equal(xlabel, "x variable")
        nt.assert_equal(ylabel, "y variable")

    def test_dropna(self):

        g = ag.JointGrid("x_na", "y", self.data, dropna=False)
        nt.assert_equal(len(g.x), len(self.x_na))

        g = ag.JointGrid("x_na", "y", self.data, dropna=True)
        nt.assert_equal(len(g.x), pd.notnull(self.x_na).sum())

    def test_axlims(self):

        lim = (-3, 3)
        g = ag.JointGrid("x", "y", self.data, xlim=lim, ylim=lim)

        nt.assert_equal(g.ax_joint.get_xlim(), lim)
        nt.assert_equal(g.ax_joint.get_ylim(), lim)

        nt.assert_equal(g.ax_marg_x.get_xlim(), lim)
        nt.assert_equal(g.ax_marg_y.get_ylim(), lim)

    def test_marginal_ticks(self):

        g = ag.JointGrid("x", "y", self.data)
        nt.assert_true(~len(g.ax_marg_x.get_xticks()))
        nt.assert_true(~len(g.ax_marg_y.get_yticks()))

    def test_bivariate_plot(self):

        g = ag.JointGrid("x", "y", self.data)
        g.plot_joint(plt.plot)

        x, y = g.ax_joint.lines[0].get_xydata().T
        npt.assert_array_equal(x, self.x)
        npt.assert_array_equal(y, self.y)

    def test_univariate_plot(self):

        g = ag.JointGrid("x", "x", self.data)
        g.plot_marginals(kdeplot)

        _, y1 = g.ax_marg_x.lines[0].get_xydata().T
        y2, _ = g.ax_marg_y.lines[0].get_xydata().T
        npt.assert_array_equal(y1, y2)

    def test_plot(self):

        g = ag.JointGrid("x", "x", self.data)
        g.plot(plt.plot, kdeplot)

        x, y = g.ax_joint.lines[0].get_xydata().T
        npt.assert_array_equal(x, self.x)
        npt.assert_array_equal(y, self.x)

        _, y1 = g.ax_marg_x.lines[0].get_xydata().T
        y2, _ = g.ax_marg_y.lines[0].get_xydata().T
        npt.assert_array_equal(y1, y2)

    def test_annotate(self):

        g = ag.JointGrid("x", "y", self.data)
        rp = stats.pearsonr(self.x, self.y)

        g.annotate(stats.pearsonr)
        annotation = g.ax_joint.legend_.texts[0].get_text()
        nt.assert_equal(annotation, "pearsonr = %.2g; p = %.2g" % rp)

        g.annotate(stats.pearsonr, stat="correlation")
        annotation = g.ax_joint.legend_.texts[0].get_text()
        nt.assert_equal(annotation, "correlation = %.2g; p = %.2g" % rp)

        def rsquared(x, y):
            return stats.pearsonr(x, y)[0] ** 2

        r2 = rsquared(self.x, self.y)
        g.annotate(rsquared)
        annotation = g.ax_joint.legend_.texts[0].get_text()
        nt.assert_equal(annotation, "rsquared = %.2g" % r2)

        template = "{stat} = {val:.3g} (p = {p:.3g})"
        g.annotate(stats.pearsonr, template=template)
        annotation = g.ax_joint.legend_.texts[0].get_text()
        nt.assert_equal(annotation, template.format(stat="pearsonr",
                                                    val=rp[0], p=rp[1]))

    def test_space(self):

        g = ag.JointGrid("x", "y", self.data, space=0)

        joint_bounds = g.ax_joint.bbox.bounds
        marg_x_bounds = g.ax_marg_x.bbox.bounds
        marg_y_bounds = g.ax_marg_y.bbox.bounds

        nt.assert_equal(joint_bounds[2], marg_x_bounds[2])
        nt.assert_equal(joint_bounds[3], marg_y_bounds[3])


class TestJointPlot(PlotTestCase):

    rs = np.random.RandomState(sum(map(ord, "jointplot")))
    x = rs.randn(100)
    y = rs.randn(100)
    data = pd.DataFrame(dict(x=x, y=y))

    def test_scatter(self):

        g = ag.jointplot("x", "y", self.data)
        nt.assert_equal(len(g.ax_joint.collections), 1)

        x, y = g.ax_joint.collections[0].get_offsets().T
        npt.assert_array_equal(self.x, x)
        npt.assert_array_equal(self.y, y)

        x_bins = _freedman_diaconis_bins(self.x)
        nt.assert_equal(len(g.ax_marg_x.patches), x_bins)

        y_bins = _freedman_diaconis_bins(self.y)
        nt.assert_equal(len(g.ax_marg_y.patches), y_bins)

    def test_reg(self):

        g = ag.jointplot("x", "y", self.data, kind="reg")
        nt.assert_equal(len(g.ax_joint.collections), 2)

        x, y = g.ax_joint.collections[0].get_offsets().T
        npt.assert_array_equal(self.x, x)
        npt.assert_array_equal(self.y, y)

        x_bins = _freedman_diaconis_bins(self.x)
        nt.assert_equal(len(g.ax_marg_x.patches), x_bins)

        y_bins = _freedman_diaconis_bins(self.y)
        nt.assert_equal(len(g.ax_marg_y.patches), y_bins)

        nt.assert_equal(len(g.ax_joint.lines), 1)
        nt.assert_equal(len(g.ax_marg_x.lines), 1)
        nt.assert_equal(len(g.ax_marg_y.lines), 1)

    def test_resid(self):

        g = ag.jointplot("x", "y", self.data, kind="resid")
        nt.assert_equal(len(g.ax_joint.collections), 1)
        nt.assert_equal(len(g.ax_joint.lines), 1)
        nt.assert_equal(len(g.ax_marg_x.lines), 0)
        nt.assert_equal(len(g.ax_marg_y.lines), 1)

    def test_hex(self):

        g = ag.jointplot("x", "y", self.data, kind="hex")
        nt.assert_equal(len(g.ax_joint.collections), 1)

        x_bins = _freedman_diaconis_bins(self.x)
        nt.assert_equal(len(g.ax_marg_x.patches), x_bins)

        y_bins = _freedman_diaconis_bins(self.y)
        nt.assert_equal(len(g.ax_marg_y.patches), y_bins)

    def test_kde(self):

        g = ag.jointplot("x", "y", self.data, kind="kde")

        nt.assert_true(len(g.ax_joint.collections) > 0)
        nt.assert_equal(len(g.ax_marg_x.collections), 1)
        nt.assert_equal(len(g.ax_marg_y.collections), 1)

        nt.assert_equal(len(g.ax_marg_x.lines), 1)
        nt.assert_equal(len(g.ax_marg_y.lines), 1)

    def test_color(self):

        g = ag.jointplot("x", "y", self.data, color="purple")

        purple = mpl.colors.colorConverter.to_rgb("purple")
        scatter_color = g.ax_joint.collections[0].get_facecolor()[0, :3]
        nt.assert_equal(tuple(scatter_color), purple)

        hist_color = g.ax_marg_x.patches[0].get_facecolor()[:3]
        nt.assert_equal(hist_color, purple)

    def test_annotation(self):

        g = ag.jointplot("x", "y", self.data)
        nt.assert_equal(len(g.ax_joint.legend_.get_texts()), 1)

        g = ag.jointplot("x", "y", self.data, stat_func=None)
        nt.assert_is(g.ax_joint.legend_, None)

    def test_hex_customise(self):

        # test that default gridsize can be overridden
        g = ag.jointplot("x", "y", self.data, kind="hex",
                         joint_kws=dict(gridsize=5))
        nt.assert_equal(len(g.ax_joint.collections), 1)
        a = g.ax_joint.collections[0].get_array()
        nt.assert_equal(28, a.shape[0])  # 28 hexagons expected for gridsize 5

    def test_bad_kind(self):

        with nt.assert_raises(ValueError):
            ag.jointplot("x", "y", self.data, kind="not_a_kind")
