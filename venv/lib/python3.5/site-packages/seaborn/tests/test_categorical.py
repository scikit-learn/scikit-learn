import numpy as np
import pandas as pd
import scipy
from scipy import stats, spatial
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import rgb2hex

from distutils.version import LooseVersion

import nose.tools as nt
import numpy.testing as npt
from numpy.testing.decorators import skipif

from . import PlotTestCase
from .. import categorical as cat
from .. import palettes


pandas_has_categoricals = LooseVersion(pd.__version__) >= "0.15"
mpl_barplot_change = LooseVersion("2.0.1")


class CategoricalFixture(PlotTestCase):
    """Test boxplot (also base class for things like violinplots)."""
    rs = np.random.RandomState(30)
    n_total = 60
    x = rs.randn(int(n_total / 3), 3)
    x_df = pd.DataFrame(x, columns=pd.Series(list("XYZ"), name="big"))
    y = pd.Series(rs.randn(n_total), name="y_data")
    g = pd.Series(np.repeat(list("abc"), int(n_total / 3)), name="small")
    h = pd.Series(np.tile(list("mn"), int(n_total / 2)), name="medium")
    u = pd.Series(np.tile(list("jkh"), int(n_total / 3)))
    df = pd.DataFrame(dict(y=y, g=g, h=h, u=u))
    x_df["W"] = g


class TestCategoricalPlotter(CategoricalFixture):

    def test_wide_df_data(self):

        p = cat._CategoricalPlotter()

        # Test basic wide DataFrame
        p.establish_variables(data=self.x_df)

        # Check data attribute
        for x, y, in zip(p.plot_data, self.x_df[["X", "Y", "Z"]].values.T):
            npt.assert_array_equal(x, y)

        # Check semantic attributes
        nt.assert_equal(p.orient, "v")
        nt.assert_is(p.plot_hues, None)
        nt.assert_is(p.group_label, "big")
        nt.assert_is(p.value_label, None)

        # Test wide dataframe with forced horizontal orientation
        p.establish_variables(data=self.x_df, orient="horiz")
        nt.assert_equal(p.orient, "h")

        # Text exception by trying to hue-group with a wide dataframe
        with nt.assert_raises(ValueError):
            p.establish_variables(hue="d", data=self.x_df)

    def test_1d_input_data(self):

        p = cat._CategoricalPlotter()

        # Test basic vector data
        x_1d_array = self.x.ravel()
        p.establish_variables(data=x_1d_array)
        nt.assert_equal(len(p.plot_data), 1)
        nt.assert_equal(len(p.plot_data[0]), self.n_total)
        nt.assert_is(p.group_label, None)
        nt.assert_is(p.value_label, None)

        # Test basic vector data in list form
        x_1d_list = x_1d_array.tolist()
        p.establish_variables(data=x_1d_list)
        nt.assert_equal(len(p.plot_data), 1)
        nt.assert_equal(len(p.plot_data[0]), self.n_total)
        nt.assert_is(p.group_label, None)
        nt.assert_is(p.value_label, None)

        # Test an object array that looks 1D but isn't
        x_notreally_1d = np.array([self.x.ravel(),
                                   self.x.ravel()[:int(self.n_total / 2)]])
        p.establish_variables(data=x_notreally_1d)
        nt.assert_equal(len(p.plot_data), 2)
        nt.assert_equal(len(p.plot_data[0]), self.n_total)
        nt.assert_equal(len(p.plot_data[1]), self.n_total / 2)
        nt.assert_is(p.group_label, None)
        nt.assert_is(p.value_label, None)

    def test_2d_input_data(self):

        p = cat._CategoricalPlotter()

        x = self.x[:, 0]

        # Test vector data that looks 2D but doesn't really have columns
        p.establish_variables(data=x[:, np.newaxis])
        nt.assert_equal(len(p.plot_data), 1)
        nt.assert_equal(len(p.plot_data[0]), self.x.shape[0])
        nt.assert_is(p.group_label, None)
        nt.assert_is(p.value_label, None)

        # Test vector data that looks 2D but doesn't really have rows
        p.establish_variables(data=x[np.newaxis, :])
        nt.assert_equal(len(p.plot_data), 1)
        nt.assert_equal(len(p.plot_data[0]), self.x.shape[0])
        nt.assert_is(p.group_label, None)
        nt.assert_is(p.value_label, None)

    def test_3d_input_data(self):

        p = cat._CategoricalPlotter()

        # Test that passing actually 3D data raises
        x = np.zeros((5, 5, 5))
        with nt.assert_raises(ValueError):
            p.establish_variables(data=x)

    def test_list_of_array_input_data(self):

        p = cat._CategoricalPlotter()

        # Test 2D input in list form
        x_list = self.x.T.tolist()
        p.establish_variables(data=x_list)
        nt.assert_equal(len(p.plot_data), 3)

        lengths = [len(v_i) for v_i in p.plot_data]
        nt.assert_equal(lengths, [self.n_total / 3] * 3)

        nt.assert_is(p.group_label, None)
        nt.assert_is(p.value_label, None)

    def test_wide_array_input_data(self):

        p = cat._CategoricalPlotter()

        # Test 2D input in array form
        p.establish_variables(data=self.x)
        nt.assert_equal(np.shape(p.plot_data), (3, self.n_total / 3))
        npt.assert_array_equal(p.plot_data, self.x.T)

        nt.assert_is(p.group_label, None)
        nt.assert_is(p.value_label, None)

    def test_single_long_direct_inputs(self):

        p = cat._CategoricalPlotter()

        # Test passing a series to the x variable
        p.establish_variables(x=self.y)
        npt.assert_equal(p.plot_data, [self.y])
        nt.assert_equal(p.orient, "h")
        nt.assert_equal(p.value_label, "y_data")
        nt.assert_is(p.group_label, None)

        # Test passing a series to the y variable
        p.establish_variables(y=self.y)
        npt.assert_equal(p.plot_data, [self.y])
        nt.assert_equal(p.orient, "v")
        nt.assert_equal(p.value_label, "y_data")
        nt.assert_is(p.group_label, None)

        # Test passing an array to the y variable
        p.establish_variables(y=self.y.values)
        npt.assert_equal(p.plot_data, [self.y])
        nt.assert_equal(p.orient, "v")
        nt.assert_is(p.value_label, None)
        nt.assert_is(p.group_label, None)

    def test_single_long_indirect_inputs(self):

        p = cat._CategoricalPlotter()

        # Test referencing a DataFrame series in the x variable
        p.establish_variables(x="y", data=self.df)
        npt.assert_equal(p.plot_data, [self.y])
        nt.assert_equal(p.orient, "h")
        nt.assert_equal(p.value_label, "y")
        nt.assert_is(p.group_label, None)

        # Test referencing a DataFrame series in the y variable
        p.establish_variables(y="y", data=self.df)
        npt.assert_equal(p.plot_data, [self.y])
        nt.assert_equal(p.orient, "v")
        nt.assert_equal(p.value_label, "y")
        nt.assert_is(p.group_label, None)

    def test_longform_groupby(self):

        p = cat._CategoricalPlotter()

        # Test a vertically oriented grouped and nested plot
        p.establish_variables("g", "y", "h", data=self.df)
        nt.assert_equal(len(p.plot_data), 3)
        nt.assert_equal(len(p.plot_hues), 3)
        nt.assert_equal(p.orient, "v")
        nt.assert_equal(p.value_label, "y")
        nt.assert_equal(p.group_label, "g")
        nt.assert_equal(p.hue_title, "h")

        for group, vals in zip(["a", "b", "c"], p.plot_data):
            npt.assert_array_equal(vals, self.y[self.g == group])

        for group, hues in zip(["a", "b", "c"], p.plot_hues):
            npt.assert_array_equal(hues, self.h[self.g == group])

        # Test a grouped and nested plot with direct array value data
        p.establish_variables("g", self.y.values, "h", self.df)
        nt.assert_is(p.value_label, None)
        nt.assert_equal(p.group_label, "g")

        for group, vals in zip(["a", "b", "c"], p.plot_data):
            npt.assert_array_equal(vals, self.y[self.g == group])

        # Test a grouped and nested plot with direct array hue data
        p.establish_variables("g", "y", self.h.values, self.df)

        for group, hues in zip(["a", "b", "c"], p.plot_hues):
            npt.assert_array_equal(hues, self.h[self.g == group])

        # Test categorical grouping data
        if pandas_has_categoricals:
            df = self.df.copy()
            df.g = df.g.astype("category")

            # Test that horizontal orientation is automatically detected
            p.establish_variables("y", "g", "h", data=df)
            nt.assert_equal(len(p.plot_data), 3)
            nt.assert_equal(len(p.plot_hues), 3)
            nt.assert_equal(p.orient, "h")
            nt.assert_equal(p.value_label, "y")
            nt.assert_equal(p.group_label, "g")
            nt.assert_equal(p.hue_title, "h")

            for group, vals in zip(["a", "b", "c"], p.plot_data):
                npt.assert_array_equal(vals, self.y[self.g == group])

            for group, hues in zip(["a", "b", "c"], p.plot_hues):
                npt.assert_array_equal(hues, self.h[self.g == group])

    def test_input_validation(self):

        p = cat._CategoricalPlotter()

        kws = dict(x="g", y="y", hue="h", units="u", data=self.df)
        for input in ["x", "y", "hue", "units"]:
            input_kws = kws.copy()
            input_kws[input] = "bad_input"
            with nt.assert_raises(ValueError):
                p.establish_variables(**input_kws)

    def test_order(self):

        p = cat._CategoricalPlotter()

        # Test inferred order from a wide dataframe input
        p.establish_variables(data=self.x_df)
        nt.assert_equal(p.group_names, ["X", "Y", "Z"])

        # Test specified order with a wide dataframe input
        p.establish_variables(data=self.x_df, order=["Y", "Z", "X"])
        nt.assert_equal(p.group_names, ["Y", "Z", "X"])

        for group, vals in zip(["Y", "Z", "X"], p.plot_data):
            npt.assert_array_equal(vals, self.x_df[group])

        with nt.assert_raises(ValueError):
            p.establish_variables(data=self.x, order=[1, 2, 0])

        # Test inferred order from a grouped longform input
        p.establish_variables("g", "y", data=self.df)
        nt.assert_equal(p.group_names, ["a", "b", "c"])

        # Test specified order from a grouped longform input
        p.establish_variables("g", "y", data=self.df, order=["b", "a", "c"])
        nt.assert_equal(p.group_names, ["b", "a", "c"])

        for group, vals in zip(["b", "a", "c"], p.plot_data):
            npt.assert_array_equal(vals, self.y[self.g == group])

        # Test inferred order from a grouped input with categorical groups
        if pandas_has_categoricals:
            df = self.df.copy()
            df.g = df.g.astype("category")
            df.g = df.g.cat.reorder_categories(["c", "b", "a"])
            p.establish_variables("g", "y", data=df)
            nt.assert_equal(p.group_names, ["c", "b", "a"])

            for group, vals in zip(["c", "b", "a"], p.plot_data):
                npt.assert_array_equal(vals, self.y[self.g == group])

            df.g = (df.g.cat.add_categories("d")
                        .cat.reorder_categories(["c", "b", "d", "a"]))
            p.establish_variables("g", "y", data=df)
            nt.assert_equal(p.group_names, ["c", "b", "d", "a"])

    def test_hue_order(self):

        p = cat._CategoricalPlotter()

        # Test inferred hue order
        p.establish_variables("g", "y", "h", data=self.df)
        nt.assert_equal(p.hue_names, ["m", "n"])

        # Test specified hue order
        p.establish_variables("g", "y", "h", data=self.df,
                              hue_order=["n", "m"])
        nt.assert_equal(p.hue_names, ["n", "m"])

        # Test inferred hue order from a categorical hue input
        if pandas_has_categoricals:
            df = self.df.copy()
            df.h = df.h.astype("category")
            df.h = df.h.cat.reorder_categories(["n", "m"])
            p.establish_variables("g", "y", "h", data=df)
            nt.assert_equal(p.hue_names, ["n", "m"])

            df.h = (df.h.cat.add_categories("o")
                        .cat.reorder_categories(["o", "m", "n"]))
            p.establish_variables("g", "y", "h", data=df)
            nt.assert_equal(p.hue_names, ["o", "m", "n"])

    def test_plot_units(self):

        p = cat._CategoricalPlotter()
        p.establish_variables("g", "y", "h", data=self.df)
        nt.assert_is(p.plot_units, None)

        p.establish_variables("g", "y", "h", data=self.df, units="u")
        for group, units in zip(["a", "b", "c"], p.plot_units):
            npt.assert_array_equal(units, self.u[self.g == group])

    def test_infer_orient(self):

        p = cat._CategoricalPlotter()

        cats = pd.Series(["a", "b", "c"] * 10)
        nums = pd.Series(self.rs.randn(30))

        nt.assert_equal(p.infer_orient(cats, nums), "v")
        nt.assert_equal(p.infer_orient(nums, cats), "h")
        nt.assert_equal(p.infer_orient(nums, None), "h")
        nt.assert_equal(p.infer_orient(None, nums), "v")
        nt.assert_equal(p.infer_orient(nums, nums, "vert"), "v")
        nt.assert_equal(p.infer_orient(nums, nums, "hori"), "h")

        with nt.assert_raises(ValueError):
            p.infer_orient(cats, cats)

        if pandas_has_categoricals:
            cats = pd.Series([0, 1, 2] * 10, dtype="category")
            nt.assert_equal(p.infer_orient(cats, nums), "v")
            nt.assert_equal(p.infer_orient(nums, cats), "h")

            with nt.assert_raises(ValueError):
                p.infer_orient(cats, cats)

    def test_default_palettes(self):

        p = cat._CategoricalPlotter()

        # Test palette mapping the x position
        p.establish_variables("g", "y", data=self.df)
        p.establish_colors(None, None, 1)
        nt.assert_equal(p.colors, palettes.color_palette(n_colors=3))

        # Test palette mapping the hue position
        p.establish_variables("g", "y", "h", data=self.df)
        p.establish_colors(None, None, 1)
        nt.assert_equal(p.colors, palettes.color_palette(n_colors=2))

    def test_default_palette_with_many_levels(self):

        with palettes.color_palette(["blue", "red"], 2):
            p = cat._CategoricalPlotter()
            p.establish_variables("g", "y", data=self.df)
            p.establish_colors(None, None, 1)
            npt.assert_array_equal(p.colors, palettes.husl_palette(3, l=.7))

    def test_specific_color(self):

        p = cat._CategoricalPlotter()

        # Test the same color for each x position
        p.establish_variables("g", "y", data=self.df)
        p.establish_colors("blue", None, 1)
        blue_rgb = mpl.colors.colorConverter.to_rgb("blue")
        nt.assert_equal(p.colors, [blue_rgb] * 3)

        # Test a color-based blend for the hue mapping
        p.establish_variables("g", "y", "h", data=self.df)
        p.establish_colors("#ff0022", None, 1)
        rgba_array = np.array(palettes.light_palette("#ff0022", 2))
        npt.assert_array_almost_equal(p.colors,
                                      rgba_array[:, :3])

    def test_specific_palette(self):

        p = cat._CategoricalPlotter()

        # Test palette mapping the x position
        p.establish_variables("g", "y", data=self.df)
        p.establish_colors(None, "dark", 1)
        nt.assert_equal(p.colors, palettes.color_palette("dark", 3))

        # Test that non-None `color` and `hue` raises an error
        p.establish_variables("g", "y", "h", data=self.df)
        p.establish_colors(None, "muted", 1)
        nt.assert_equal(p.colors, palettes.color_palette("muted", 2))

        # Test that specified palette overrides specified color
        p = cat._CategoricalPlotter()
        p.establish_variables("g", "y", data=self.df)
        p.establish_colors("blue", "deep", 1)
        nt.assert_equal(p.colors, palettes.color_palette("deep", 3))

    def test_dict_as_palette(self):

        p = cat._CategoricalPlotter()
        p.establish_variables("g", "y", "h", data=self.df)
        pal = {"m": (0, 0, 1), "n": (1, 0, 0)}
        p.establish_colors(None, pal, 1)
        nt.assert_equal(p.colors, [(0, 0, 1), (1, 0, 0)])

    def test_palette_desaturation(self):

        p = cat._CategoricalPlotter()
        p.establish_variables("g", "y", data=self.df)
        p.establish_colors((0, 0, 1), None, .5)
        nt.assert_equal(p.colors, [(.25, .25, .75)] * 3)

        p.establish_colors(None, [(0, 0, 1), (1, 0, 0), "w"], .5)
        nt.assert_equal(p.colors, [(.25, .25, .75),
                                   (.75, .25, .25),
                                   (1, 1, 1)])


class TestCategoricalStatPlotter(CategoricalFixture):

    def test_no_bootstrappig(self):

        p = cat._CategoricalStatPlotter()
        p.establish_variables("g", "y", data=self.df)
        p.estimate_statistic(np.mean, None, 100)
        npt.assert_array_equal(p.confint, np.array([]))

        p.establish_variables("g", "y", "h", data=self.df)
        p.estimate_statistic(np.mean, None, 100)
        npt.assert_array_equal(p.confint, np.array([[], [], []]))

    def test_single_layer_stats(self):

        p = cat._CategoricalStatPlotter()

        g = pd.Series(np.repeat(list("abc"), 100))
        y = pd.Series(np.random.RandomState(0).randn(300))

        p.establish_variables(g, y)
        p.estimate_statistic(np.mean, 95, 10000)

        nt.assert_equal(p.statistic.shape, (3,))
        nt.assert_equal(p.confint.shape, (3, 2))

        npt.assert_array_almost_equal(p.statistic,
                                      y.groupby(g).mean())

        for ci, (_, grp_y) in zip(p.confint, y.groupby(g)):
            sem = stats.sem(grp_y)
            mean = grp_y.mean()
            stats.norm.ppf(.975)
            half_ci = stats.norm.ppf(.975) * sem
            ci_want = mean - half_ci, mean + half_ci
            npt.assert_array_almost_equal(ci_want, ci, 2)

    def test_single_layer_stats_with_units(self):

        p = cat._CategoricalStatPlotter()

        g = pd.Series(np.repeat(list("abc"), 90))
        y = pd.Series(np.random.RandomState(0).randn(270))
        u = pd.Series(np.repeat(np.tile(list("xyz"), 30), 3))
        y[u == "x"] -= 3
        y[u == "y"] += 3

        p.establish_variables(g, y)
        p.estimate_statistic(np.mean, 95, 10000)
        stat1, ci1 = p.statistic, p.confint

        p.establish_variables(g, y, units=u)
        p.estimate_statistic(np.mean, 95, 10000)
        stat2, ci2 = p.statistic, p.confint

        npt.assert_array_equal(stat1, stat2)
        ci1_size = ci1[:, 1] - ci1[:, 0]
        ci2_size = ci2[:, 1] - ci2[:, 0]
        npt.assert_array_less(ci1_size, ci2_size)

    def test_single_layer_stats_with_missing_data(self):

        p = cat._CategoricalStatPlotter()

        g = pd.Series(np.repeat(list("abc"), 100))
        y = pd.Series(np.random.RandomState(0).randn(300))

        p.establish_variables(g, y, order=list("abdc"))
        p.estimate_statistic(np.mean, 95, 10000)

        nt.assert_equal(p.statistic.shape, (4,))
        nt.assert_equal(p.confint.shape, (4, 2))

        mean = y[g == "b"].mean()
        sem = stats.sem(y[g == "b"])
        half_ci = stats.norm.ppf(.975) * sem
        ci = mean - half_ci, mean + half_ci
        npt.assert_almost_equal(p.statistic[1], mean)
        npt.assert_array_almost_equal(p.confint[1], ci, 2)

        npt.assert_equal(p.statistic[2], np.nan)
        npt.assert_array_equal(p.confint[2], (np.nan, np.nan))

    def test_nested_stats(self):

        p = cat._CategoricalStatPlotter()

        g = pd.Series(np.repeat(list("abc"), 100))
        h = pd.Series(np.tile(list("xy"), 150))
        y = pd.Series(np.random.RandomState(0).randn(300))

        p.establish_variables(g, y, h)
        p.estimate_statistic(np.mean, 95, 50000)

        nt.assert_equal(p.statistic.shape, (3, 2))
        nt.assert_equal(p.confint.shape, (3, 2, 2))

        npt.assert_array_almost_equal(p.statistic,
                                      y.groupby([g, h]).mean().unstack())

        for ci_g, (_, grp_y) in zip(p.confint, y.groupby(g)):
            for ci, hue_y in zip(ci_g, [grp_y[::2], grp_y[1::2]]):
                sem = stats.sem(hue_y)
                mean = hue_y.mean()
                half_ci = stats.norm.ppf(.975) * sem
                ci_want = mean - half_ci, mean + half_ci
                npt.assert_array_almost_equal(ci_want, ci, 2)

    def test_nested_stats_with_units(self):

        p = cat._CategoricalStatPlotter()

        g = pd.Series(np.repeat(list("abc"), 90))
        h = pd.Series(np.tile(list("xy"), 135))
        u = pd.Series(np.repeat(list("ijkijk"), 45))
        y = pd.Series(np.random.RandomState(0).randn(270))
        y[u == "i"] -= 3
        y[u == "k"] += 3

        p.establish_variables(g, y, h)
        p.estimate_statistic(np.mean, 95, 10000)
        stat1, ci1 = p.statistic, p.confint

        p.establish_variables(g, y, h, units=u)
        p.estimate_statistic(np.mean, 95, 10000)
        stat2, ci2 = p.statistic, p.confint

        npt.assert_array_equal(stat1, stat2)
        ci1_size = ci1[:, 0, 1] - ci1[:, 0, 0]
        ci2_size = ci2[:, 0, 1] - ci2[:, 0, 0]
        npt.assert_array_less(ci1_size, ci2_size)

    def test_nested_stats_with_missing_data(self):

        p = cat._CategoricalStatPlotter()

        g = pd.Series(np.repeat(list("abc"), 100))
        y = pd.Series(np.random.RandomState(0).randn(300))
        h = pd.Series(np.tile(list("xy"), 150))

        p.establish_variables(g, y, h,
                              order=list("abdc"),
                              hue_order=list("zyx"))
        p.estimate_statistic(np.mean, 95, 50000)

        nt.assert_equal(p.statistic.shape, (4, 3))
        nt.assert_equal(p.confint.shape, (4, 3, 2))

        mean = y[(g == "b") & (h == "x")].mean()
        sem = stats.sem(y[(g == "b") & (h == "x")])
        half_ci = stats.norm.ppf(.975) * sem
        ci = mean - half_ci, mean + half_ci
        npt.assert_almost_equal(p.statistic[1, 2], mean)
        npt.assert_array_almost_equal(p.confint[1, 2], ci, 2)

        npt.assert_array_equal(p.statistic[:, 0], [np.nan] * 4)
        npt.assert_array_equal(p.statistic[2], [np.nan] * 3)
        npt.assert_array_equal(p.confint[:, 0],
                               np.zeros((4, 2)) * np.nan)
        npt.assert_array_equal(p.confint[2],
                               np.zeros((3, 2)) * np.nan)

    def test_sd_error_bars(self):

        p = cat._CategoricalStatPlotter()

        g = pd.Series(np.repeat(list("abc"), 100))
        y = pd.Series(np.random.RandomState(0).randn(300))

        p.establish_variables(g, y)
        p.estimate_statistic(np.mean, "sd", None)

        nt.assert_equal(p.statistic.shape, (3,))
        nt.assert_equal(p.confint.shape, (3, 2))

        npt.assert_array_almost_equal(p.statistic,
                                      y.groupby(g).mean())

        for ci, (_, grp_y) in zip(p.confint, y.groupby(g)):
            mean = grp_y.mean()
            half_ci = np.std(grp_y)
            ci_want = mean - half_ci, mean + half_ci
            npt.assert_array_almost_equal(ci_want, ci, 2)

    def test_nested_sd_error_bars(self):

        p = cat._CategoricalStatPlotter()

        g = pd.Series(np.repeat(list("abc"), 100))
        h = pd.Series(np.tile(list("xy"), 150))
        y = pd.Series(np.random.RandomState(0).randn(300))

        p.establish_variables(g, y, h)
        p.estimate_statistic(np.mean, "sd", None)

        nt.assert_equal(p.statistic.shape, (3, 2))
        nt.assert_equal(p.confint.shape, (3, 2, 2))

        npt.assert_array_almost_equal(p.statistic,
                                      y.groupby([g, h]).mean().unstack())

        for ci_g, (_, grp_y) in zip(p.confint, y.groupby(g)):
            for ci, hue_y in zip(ci_g, [grp_y[::2], grp_y[1::2]]):
                mean = hue_y.mean()
                half_ci = np.std(hue_y)
                ci_want = mean - half_ci, mean + half_ci
                npt.assert_array_almost_equal(ci_want, ci, 2)

    def test_draw_cis(self):

        p = cat._CategoricalStatPlotter()

        # Test vertical CIs
        p.orient = "v"

        f, ax = plt.subplots()
        at_group = [0, 1]
        confints = [(.5, 1.5), (.25, .8)]
        colors = [".2", ".3"]
        p.draw_confints(ax, at_group, confints, colors)

        lines = ax.lines
        for line, at, ci, c in zip(lines, at_group, confints, colors):
            x, y = line.get_xydata().T
            npt.assert_array_equal(x, [at, at])
            npt.assert_array_equal(y, ci)
            nt.assert_equal(line.get_color(), c)

        plt.close("all")

        # Test horizontal CIs
        p.orient = "h"

        f, ax = plt.subplots()
        p.draw_confints(ax, at_group, confints, colors)

        lines = ax.lines
        for line, at, ci, c in zip(lines, at_group, confints, colors):
            x, y = line.get_xydata().T
            npt.assert_array_equal(x, ci)
            npt.assert_array_equal(y, [at, at])
            nt.assert_equal(line.get_color(), c)

        plt.close("all")

        # Test vertical CIs with endcaps
        p.orient = "v"

        f, ax = plt.subplots()
        p.draw_confints(ax, at_group, confints, colors, capsize=0.3)
        capline = ax.lines[len(ax.lines) - 1]
        caplinestart = capline.get_xdata()[0]
        caplineend = capline.get_xdata()[1]
        caplinelength = abs(caplineend - caplinestart)
        nt.assert_almost_equal(caplinelength, 0.3)
        nt.assert_equal(len(ax.lines), 6)

        plt.close("all")

        # Test horizontal CIs with endcaps
        p.orient = "h"

        f, ax = plt.subplots()
        p.draw_confints(ax, at_group, confints, colors, capsize=0.3)
        capline = ax.lines[len(ax.lines) - 1]
        caplinestart = capline.get_ydata()[0]
        caplineend = capline.get_ydata()[1]
        caplinelength = abs(caplineend - caplinestart)
        nt.assert_almost_equal(caplinelength, 0.3)
        nt.assert_equal(len(ax.lines), 6)

        # Test extra keyword arguments
        f, ax = plt.subplots()
        p.draw_confints(ax, at_group, confints, colors, lw=4)
        line = ax.lines[0]
        nt.assert_equal(line.get_linewidth(), 4)

        plt.close("all")

        # Test errwidth is set appropriately
        f, ax = plt.subplots()
        p.draw_confints(ax, at_group, confints, colors, errwidth=2)
        capline = ax.lines[len(ax.lines)-1]
        nt.assert_equal(capline._linewidth, 2)
        nt.assert_equal(len(ax.lines), 2)

        plt.close("all")


class TestBoxPlotter(CategoricalFixture):

    default_kws = dict(x=None, y=None, hue=None, data=None,
                       order=None, hue_order=None,
                       orient=None, color=None, palette=None,
                       saturation=.75, width=.8, dodge=True,
                       fliersize=5, linewidth=None)

    def test_nested_width(self):

        kws = self.default_kws.copy()
        p = cat._BoxPlotter(**kws)
        p.establish_variables("g", "y", "h", data=self.df)
        nt.assert_equal(p.nested_width, .4 * .98)

        kws = self.default_kws.copy()
        kws["width"] = .6
        p = cat._BoxPlotter(**kws)
        p.establish_variables("g", "y", "h", data=self.df)
        nt.assert_equal(p.nested_width, .3 * .98)

        kws = self.default_kws.copy()
        kws["dodge"] = False
        p = cat._BoxPlotter(**kws)
        p.establish_variables("g", "y", "h", data=self.df)
        nt.assert_equal(p.nested_width, .8)

    def test_hue_offsets(self):

        p = cat._BoxPlotter(**self.default_kws)
        p.establish_variables("g", "y", "h", data=self.df)
        npt.assert_array_equal(p.hue_offsets, [-.2, .2])

        kws = self.default_kws.copy()
        kws["width"] = .6
        p = cat._BoxPlotter(**kws)
        p.establish_variables("g", "y", "h", data=self.df)
        npt.assert_array_equal(p.hue_offsets, [-.15, .15])

        p = cat._BoxPlotter(**kws)
        p.establish_variables("h", "y", "g", data=self.df)
        npt.assert_array_almost_equal(p.hue_offsets, [-.2, 0, .2])

    def test_axes_data(self):

        ax = cat.boxplot("g", "y", data=self.df)
        nt.assert_equal(len(ax.artists), 3)

        plt.close("all")

        ax = cat.boxplot("g", "y", "h", data=self.df)
        nt.assert_equal(len(ax.artists), 6)

        plt.close("all")

    def test_box_colors(self):

        ax = cat.boxplot("g", "y", data=self.df, saturation=1)
        pal = palettes.color_palette(n_colors=3)
        for patch, color in zip(ax.artists, pal):
            nt.assert_equal(patch.get_facecolor()[:3], color)

        plt.close("all")

        ax = cat.boxplot("g", "y", "h", data=self.df, saturation=1)
        pal = palettes.color_palette(n_colors=2)
        for patch, color in zip(ax.artists, pal * 2):
            nt.assert_equal(patch.get_facecolor()[:3], color)

        plt.close("all")

    def test_draw_missing_boxes(self):

        ax = cat.boxplot("g", "y", data=self.df,
                         order=["a", "b", "c", "d"])
        nt.assert_equal(len(ax.artists), 3)

    def test_missing_data(self):

        x = ["a", "a", "b", "b", "c", "c", "d", "d"]
        h = ["x", "y", "x", "y", "x", "y", "x", "y"]
        y = self.rs.randn(8)
        y[-2:] = np.nan

        ax = cat.boxplot(x, y)
        nt.assert_equal(len(ax.artists), 3)

        plt.close("all")

        y[-1] = 0
        ax = cat.boxplot(x, y, h)
        nt.assert_equal(len(ax.artists), 7)

        plt.close("all")

    def test_boxplots(self):

        # Smoke test the high level boxplot options

        cat.boxplot("y", data=self.df)
        plt.close("all")

        cat.boxplot(y="y", data=self.df)
        plt.close("all")

        cat.boxplot("g", "y", data=self.df)
        plt.close("all")

        cat.boxplot("y", "g", data=self.df, orient="h")
        plt.close("all")

        cat.boxplot("g", "y", "h", data=self.df)
        plt.close("all")

        cat.boxplot("g", "y", "h", order=list("nabc"), data=self.df)
        plt.close("all")

        cat.boxplot("g", "y", "h", hue_order=list("omn"), data=self.df)
        plt.close("all")

        cat.boxplot("y", "g", "h", data=self.df, orient="h")
        plt.close("all")

    def test_axes_annotation(self):

        ax = cat.boxplot("g", "y", data=self.df)
        nt.assert_equal(ax.get_xlabel(), "g")
        nt.assert_equal(ax.get_ylabel(), "y")
        nt.assert_equal(ax.get_xlim(), (-.5, 2.5))
        npt.assert_array_equal(ax.get_xticks(), [0, 1, 2])
        npt.assert_array_equal([l.get_text() for l in ax.get_xticklabels()],
                               ["a", "b", "c"])

        plt.close("all")

        ax = cat.boxplot("g", "y", "h", data=self.df)
        nt.assert_equal(ax.get_xlabel(), "g")
        nt.assert_equal(ax.get_ylabel(), "y")
        npt.assert_array_equal(ax.get_xticks(), [0, 1, 2])
        npt.assert_array_equal([l.get_text() for l in ax.get_xticklabels()],
                               ["a", "b", "c"])
        npt.assert_array_equal([l.get_text() for l in ax.legend_.get_texts()],
                               ["m", "n"])

        plt.close("all")

        ax = cat.boxplot("y", "g", data=self.df, orient="h")
        nt.assert_equal(ax.get_xlabel(), "y")
        nt.assert_equal(ax.get_ylabel(), "g")
        nt.assert_equal(ax.get_ylim(), (2.5, -.5))
        npt.assert_array_equal(ax.get_yticks(), [0, 1, 2])
        npt.assert_array_equal([l.get_text() for l in ax.get_yticklabels()],
                               ["a", "b", "c"])

        plt.close("all")


class TestViolinPlotter(CategoricalFixture):

    default_kws = dict(x=None, y=None, hue=None, data=None,
                       order=None, hue_order=None,
                       bw="scott", cut=2, scale="area", scale_hue=True,
                       gridsize=100, width=.8, inner="box", split=False,
                       dodge=True, orient=None, linewidth=None,
                       color=None, palette=None, saturation=.75)

    def test_split_error(self):

        kws = self.default_kws.copy()
        kws.update(dict(x="h", y="y", hue="g", data=self.df, split=True))

        with nt.assert_raises(ValueError):
            cat._ViolinPlotter(**kws)

    def test_no_observations(self):

        p = cat._ViolinPlotter(**self.default_kws)

        x = ["a", "a", "b"]
        y = self.rs.randn(3)
        y[-1] = np.nan
        p.establish_variables(x, y)
        p.estimate_densities("scott", 2, "area", True, 20)

        nt.assert_equal(len(p.support[0]), 20)
        nt.assert_equal(len(p.support[1]), 0)

        nt.assert_equal(len(p.density[0]), 20)
        nt.assert_equal(len(p.density[1]), 1)

        nt.assert_equal(p.density[1].item(), 1)

        p.estimate_densities("scott", 2, "count", True, 20)
        nt.assert_equal(p.density[1].item(), 0)

        x = ["a"] * 4 + ["b"] * 2
        y = self.rs.randn(6)
        h = ["m", "n"] * 2 + ["m"] * 2

        p.establish_variables(x, y, h)
        p.estimate_densities("scott", 2, "area", True, 20)

        nt.assert_equal(len(p.support[1][0]), 20)
        nt.assert_equal(len(p.support[1][1]), 0)

        nt.assert_equal(len(p.density[1][0]), 20)
        nt.assert_equal(len(p.density[1][1]), 1)

        nt.assert_equal(p.density[1][1].item(), 1)

        p.estimate_densities("scott", 2, "count", False, 20)
        nt.assert_equal(p.density[1][1].item(), 0)

    def test_single_observation(self):

        p = cat._ViolinPlotter(**self.default_kws)

        x = ["a", "a", "b"]
        y = self.rs.randn(3)
        p.establish_variables(x, y)
        p.estimate_densities("scott", 2, "area", True, 20)

        nt.assert_equal(len(p.support[0]), 20)
        nt.assert_equal(len(p.support[1]), 1)

        nt.assert_equal(len(p.density[0]), 20)
        nt.assert_equal(len(p.density[1]), 1)

        nt.assert_equal(p.density[1].item(), 1)

        p.estimate_densities("scott", 2, "count", True, 20)
        nt.assert_equal(p.density[1].item(), .5)

        x = ["b"] * 4 + ["a"] * 3
        y = self.rs.randn(7)
        h = (["m", "n"] * 4)[:-1]

        p.establish_variables(x, y, h)
        p.estimate_densities("scott", 2, "area", True, 20)

        nt.assert_equal(len(p.support[1][0]), 20)
        nt.assert_equal(len(p.support[1][1]), 1)

        nt.assert_equal(len(p.density[1][0]), 20)
        nt.assert_equal(len(p.density[1][1]), 1)

        nt.assert_equal(p.density[1][1].item(), 1)

        p.estimate_densities("scott", 2, "count", False, 20)
        nt.assert_equal(p.density[1][1].item(), .5)

    def test_dwidth(self):

        kws = self.default_kws.copy()
        kws.update(dict(x="g", y="y", data=self.df))

        p = cat._ViolinPlotter(**kws)
        nt.assert_equal(p.dwidth, .4)

        kws.update(dict(width=.4))
        p = cat._ViolinPlotter(**kws)
        nt.assert_equal(p.dwidth, .2)

        kws.update(dict(hue="h", width=.8))
        p = cat._ViolinPlotter(**kws)
        nt.assert_equal(p.dwidth, .2)

        kws.update(dict(split=True))
        p = cat._ViolinPlotter(**kws)
        nt.assert_equal(p.dwidth, .4)

    def test_scale_area(self):

        kws = self.default_kws.copy()
        kws["scale"] = "area"
        p = cat._ViolinPlotter(**kws)

        # Test single layer of grouping
        p.hue_names = None
        density = [self.rs.uniform(0, .8, 50), self.rs.uniform(0, .2, 50)]
        max_before = np.array([d.max() for d in density])
        p.scale_area(density, max_before, False)
        max_after = np.array([d.max() for d in density])
        nt.assert_equal(max_after[0], 1)

        before_ratio = max_before[1] / max_before[0]
        after_ratio = max_after[1] / max_after[0]
        nt.assert_equal(before_ratio, after_ratio)

        # Test nested grouping scaling across all densities
        p.hue_names = ["foo", "bar"]
        density = [[self.rs.uniform(0, .8, 50), self.rs.uniform(0, .2, 50)],
                   [self.rs.uniform(0, .1, 50), self.rs.uniform(0, .02, 50)]]

        max_before = np.array([[r.max() for r in row] for row in density])
        p.scale_area(density, max_before, False)
        max_after = np.array([[r.max() for r in row] for row in density])
        nt.assert_equal(max_after[0, 0], 1)

        before_ratio = max_before[1, 1] / max_before[0, 0]
        after_ratio = max_after[1, 1] / max_after[0, 0]
        nt.assert_equal(before_ratio, after_ratio)

        # Test nested grouping scaling within hue
        p.hue_names = ["foo", "bar"]
        density = [[self.rs.uniform(0, .8, 50), self.rs.uniform(0, .2, 50)],
                   [self.rs.uniform(0, .1, 50), self.rs.uniform(0, .02, 50)]]

        max_before = np.array([[r.max() for r in row] for row in density])
        p.scale_area(density, max_before, True)
        max_after = np.array([[r.max() for r in row] for row in density])
        nt.assert_equal(max_after[0, 0], 1)
        nt.assert_equal(max_after[1, 0], 1)

        before_ratio = max_before[1, 1] / max_before[1, 0]
        after_ratio = max_after[1, 1] / max_after[1, 0]
        nt.assert_equal(before_ratio, after_ratio)

    def test_scale_width(self):

        kws = self.default_kws.copy()
        kws["scale"] = "width"
        p = cat._ViolinPlotter(**kws)

        # Test single layer of grouping
        p.hue_names = None
        density = [self.rs.uniform(0, .8, 50), self.rs.uniform(0, .2, 50)]
        p.scale_width(density)
        max_after = np.array([d.max() for d in density])
        npt.assert_array_equal(max_after, [1, 1])

        # Test nested grouping
        p.hue_names = ["foo", "bar"]
        density = [[self.rs.uniform(0, .8, 50), self.rs.uniform(0, .2, 50)],
                   [self.rs.uniform(0, .1, 50), self.rs.uniform(0, .02, 50)]]

        p.scale_width(density)
        max_after = np.array([[r.max() for r in row] for row in density])
        npt.assert_array_equal(max_after, [[1, 1], [1, 1]])

    def test_scale_count(self):

        kws = self.default_kws.copy()
        kws["scale"] = "count"
        p = cat._ViolinPlotter(**kws)

        # Test single layer of grouping
        p.hue_names = None
        density = [self.rs.uniform(0, .8, 20), self.rs.uniform(0, .2, 40)]
        counts = np.array([20, 40])
        p.scale_count(density, counts, False)
        max_after = np.array([d.max() for d in density])
        npt.assert_array_equal(max_after, [.5, 1])

        # Test nested grouping scaling across all densities
        p.hue_names = ["foo", "bar"]
        density = [[self.rs.uniform(0, .8, 5), self.rs.uniform(0, .2, 40)],
                   [self.rs.uniform(0, .1, 100), self.rs.uniform(0, .02, 50)]]

        counts = np.array([[5, 40], [100, 50]])
        p.scale_count(density, counts, False)
        max_after = np.array([[r.max() for r in row] for row in density])
        npt.assert_array_equal(max_after, [[.05, .4], [1, .5]])

        # Test nested grouping scaling within hue
        p.hue_names = ["foo", "bar"]
        density = [[self.rs.uniform(0, .8, 5), self.rs.uniform(0, .2, 40)],
                   [self.rs.uniform(0, .1, 100), self.rs.uniform(0, .02, 50)]]

        counts = np.array([[5, 40], [100, 50]])
        p.scale_count(density, counts, True)
        max_after = np.array([[r.max() for r in row] for row in density])
        npt.assert_array_equal(max_after, [[.125, 1], [1, .5]])

    def test_bad_scale(self):

        kws = self.default_kws.copy()
        kws["scale"] = "not_a_scale_type"
        with nt.assert_raises(ValueError):
            cat._ViolinPlotter(**kws)

    def test_kde_fit(self):

        p = cat._ViolinPlotter(**self.default_kws)
        data = self.y
        data_std = data.std(ddof=1)

        # Bandwidth behavior depends on scipy version
        if LooseVersion(scipy.__version__) < "0.11":
            # Test ignoring custom bandwidth on old scipy
            kde, bw = p.fit_kde(self.y, .2)
            nt.assert_is_instance(kde, stats.gaussian_kde)
            nt.assert_equal(kde.factor, kde.scotts_factor())

        else:
            # Test reference rule bandwidth
            kde, bw = p.fit_kde(data, "scott")
            nt.assert_is_instance(kde, stats.gaussian_kde)
            nt.assert_equal(kde.factor, kde.scotts_factor())
            nt.assert_equal(bw, kde.scotts_factor() * data_std)

            # Test numeric scale factor
            kde, bw = p.fit_kde(self.y, .2)
            nt.assert_is_instance(kde, stats.gaussian_kde)
            nt.assert_equal(kde.factor, .2)
            nt.assert_equal(bw, .2 * data_std)

    def test_draw_to_density(self):

        p = cat._ViolinPlotter(**self.default_kws)
        # p.dwidth will be 1 for easier testing
        p.width = 2

        # Test verical plots
        support = np.array([.2, .6])
        density = np.array([.1, .4])

        # Test full vertical plot
        _, ax = plt.subplots()
        p.draw_to_density(ax, 0, .5, support, density, False)
        x, y = ax.lines[0].get_xydata().T
        npt.assert_array_equal(x, [.99 * -.4, .99 * .4])
        npt.assert_array_equal(y, [.5, .5])
        plt.close("all")

        # Test left vertical plot
        _, ax = plt.subplots()
        p.draw_to_density(ax, 0, .5, support, density, "left")
        x, y = ax.lines[0].get_xydata().T
        npt.assert_array_equal(x, [.99 * -.4, 0])
        npt.assert_array_equal(y, [.5, .5])
        plt.close("all")

        # Test right vertical plot
        _, ax = plt.subplots()
        p.draw_to_density(ax, 0, .5, support, density, "right")
        x, y = ax.lines[0].get_xydata().T
        npt.assert_array_equal(x, [0, .99 * .4])
        npt.assert_array_equal(y, [.5, .5])
        plt.close("all")

        # Switch orientation to test horizontal plots
        p.orient = "h"
        support = np.array([.2, .5])
        density = np.array([.3, .7])

        # Test full horizontal plot
        _, ax = plt.subplots()
        p.draw_to_density(ax, 0, .6, support, density, False)
        x, y = ax.lines[0].get_xydata().T
        npt.assert_array_equal(x, [.6, .6])
        npt.assert_array_equal(y, [.99 * -.7, .99 * .7])
        plt.close("all")

        # Test left horizontal plot
        _, ax = plt.subplots()
        p.draw_to_density(ax, 0, .6, support, density, "left")
        x, y = ax.lines[0].get_xydata().T
        npt.assert_array_equal(x, [.6, .6])
        npt.assert_array_equal(y, [.99 * -.7, 0])
        plt.close("all")

        # Test right horizontal plot
        _, ax = plt.subplots()
        p.draw_to_density(ax, 0, .6, support, density, "right")
        x, y = ax.lines[0].get_xydata().T
        npt.assert_array_equal(x, [.6, .6])
        npt.assert_array_equal(y, [0, .99 * .7])
        plt.close("all")

    def test_draw_single_observations(self):

        p = cat._ViolinPlotter(**self.default_kws)
        p.width = 2

        # Test vertical plot
        _, ax = plt.subplots()
        p.draw_single_observation(ax, 1, 1.5, 1)
        x, y = ax.lines[0].get_xydata().T
        npt.assert_array_equal(x, [0, 2])
        npt.assert_array_equal(y, [1.5, 1.5])
        plt.close("all")

        # Test horizontal plot
        p.orient = "h"
        _, ax = plt.subplots()
        p.draw_single_observation(ax, 2, 2.2, .5)
        x, y = ax.lines[0].get_xydata().T
        npt.assert_array_equal(x, [2.2, 2.2])
        npt.assert_array_equal(y, [1.5, 2.5])
        plt.close("all")

    def test_draw_box_lines(self):

        # Test vertical plot
        kws = self.default_kws.copy()
        kws.update(dict(y="y", data=self.df, inner=None))
        p = cat._ViolinPlotter(**kws)

        _, ax = plt.subplots()
        p.draw_box_lines(ax, self.y, p.support[0], p.density[0], 0)
        nt.assert_equal(len(ax.lines), 2)

        q25, q50, q75 = np.percentile(self.y, [25, 50, 75])
        _, y = ax.lines[1].get_xydata().T
        npt.assert_array_equal(y, [q25, q75])

        _, y = ax.collections[0].get_offsets().T
        nt.assert_equal(y, q50)

        plt.close("all")

        # Test horizontal plot
        kws = self.default_kws.copy()
        kws.update(dict(x="y", data=self.df, inner=None))
        p = cat._ViolinPlotter(**kws)

        _, ax = plt.subplots()
        p.draw_box_lines(ax, self.y, p.support[0], p.density[0], 0)
        nt.assert_equal(len(ax.lines), 2)

        q25, q50, q75 = np.percentile(self.y, [25, 50, 75])
        x, _ = ax.lines[1].get_xydata().T
        npt.assert_array_equal(x, [q25, q75])

        x, _ = ax.collections[0].get_offsets().T
        nt.assert_equal(x, q50)

        plt.close("all")

    def test_draw_quartiles(self):

        kws = self.default_kws.copy()
        kws.update(dict(y="y", data=self.df, inner=None))
        p = cat._ViolinPlotter(**kws)

        _, ax = plt.subplots()
        p.draw_quartiles(ax, self.y, p.support[0], p.density[0], 0)
        for val, line in zip(np.percentile(self.y, [25, 50, 75]), ax.lines):
            _, y = line.get_xydata().T
            npt.assert_array_equal(y, [val, val])

    def test_draw_points(self):

        p = cat._ViolinPlotter(**self.default_kws)

        # Test vertical plot
        _, ax = plt.subplots()
        p.draw_points(ax, self.y, 0)
        x, y = ax.collections[0].get_offsets().T
        npt.assert_array_equal(x, np.zeros_like(self.y))
        npt.assert_array_equal(y, self.y)
        plt.close("all")

        # Test horizontal plot
        p.orient = "h"
        _, ax = plt.subplots()
        p.draw_points(ax, self.y, 0)
        x, y = ax.collections[0].get_offsets().T
        npt.assert_array_equal(x, self.y)
        npt.assert_array_equal(y, np.zeros_like(self.y))
        plt.close("all")

    def test_draw_sticks(self):

        kws = self.default_kws.copy()
        kws.update(dict(y="y", data=self.df, inner=None))
        p = cat._ViolinPlotter(**kws)

        # Test vertical plot
        _, ax = plt.subplots()
        p.draw_stick_lines(ax, self.y, p.support[0], p.density[0], 0)
        for val, line in zip(self.y, ax.lines):
            _, y = line.get_xydata().T
            npt.assert_array_equal(y, [val, val])
        plt.close("all")

        # Test horizontal plot
        p.orient = "h"
        _, ax = plt.subplots()
        p.draw_stick_lines(ax, self.y, p.support[0], p.density[0], 0)
        for val, line in zip(self.y, ax.lines):
            x, _ = line.get_xydata().T
            npt.assert_array_equal(x, [val, val])
        plt.close("all")

    def test_validate_inner(self):

        kws = self.default_kws.copy()
        kws.update(dict(inner="bad_inner"))
        with nt.assert_raises(ValueError):
            cat._ViolinPlotter(**kws)

    def test_draw_violinplots(self):

        kws = self.default_kws.copy()

        # Test single vertical violin
        kws.update(dict(y="y", data=self.df, inner=None,
                        saturation=1, color=(1, 0, 0, 1)))
        p = cat._ViolinPlotter(**kws)

        _, ax = plt.subplots()
        p.draw_violins(ax)
        nt.assert_equal(len(ax.collections), 1)
        npt.assert_array_equal(ax.collections[0].get_facecolors(),
                               [(1, 0, 0, 1)])
        plt.close("all")

        # Test single horizontal violin
        kws.update(dict(x="y", y=None, color=(0, 1, 0, 1)))
        p = cat._ViolinPlotter(**kws)

        _, ax = plt.subplots()
        p.draw_violins(ax)
        nt.assert_equal(len(ax.collections), 1)
        npt.assert_array_equal(ax.collections[0].get_facecolors(),
                               [(0, 1, 0, 1)])
        plt.close("all")

        # Test multiple vertical violins
        kws.update(dict(x="g", y="y", color=None,))
        p = cat._ViolinPlotter(**kws)

        _, ax = plt.subplots()
        p.draw_violins(ax)
        nt.assert_equal(len(ax.collections), 3)
        for violin, color in zip(ax.collections, palettes.color_palette()):
            npt.assert_array_equal(violin.get_facecolors()[0, :-1], color)
        plt.close("all")

        # Test multiple violins with hue nesting
        kws.update(dict(hue="h"))
        p = cat._ViolinPlotter(**kws)

        _, ax = plt.subplots()
        p.draw_violins(ax)
        nt.assert_equal(len(ax.collections), 6)
        for violin, color in zip(ax.collections,
                                 palettes.color_palette(n_colors=2) * 3):
            npt.assert_array_equal(violin.get_facecolors()[0, :-1], color)
        plt.close("all")

        # Test multiple split violins
        kws.update(dict(split=True, palette="muted"))
        p = cat._ViolinPlotter(**kws)

        _, ax = plt.subplots()
        p.draw_violins(ax)
        nt.assert_equal(len(ax.collections), 6)
        for violin, color in zip(ax.collections,
                                 palettes.color_palette("muted",
                                                        n_colors=2) * 3):
            npt.assert_array_equal(violin.get_facecolors()[0, :-1], color)
        plt.close("all")

    def test_draw_violinplots_no_observations(self):

        kws = self.default_kws.copy()
        kws["inner"] = None

        # Test single layer of grouping
        x = ["a", "a", "b"]
        y = self.rs.randn(3)
        y[-1] = np.nan
        kws.update(x=x, y=y)
        p = cat._ViolinPlotter(**kws)

        _, ax = plt.subplots()
        p.draw_violins(ax)
        nt.assert_equal(len(ax.collections), 1)
        nt.assert_equal(len(ax.lines), 0)
        plt.close("all")

        # Test nested hue grouping
        x = ["a"] * 4 + ["b"] * 2
        y = self.rs.randn(6)
        h = ["m", "n"] * 2 + ["m"] * 2
        kws.update(x=x, y=y, hue=h)
        p = cat._ViolinPlotter(**kws)

        _, ax = plt.subplots()
        p.draw_violins(ax)
        nt.assert_equal(len(ax.collections), 3)
        nt.assert_equal(len(ax.lines), 0)
        plt.close("all")

    def test_draw_violinplots_single_observations(self):

        kws = self.default_kws.copy()
        kws["inner"] = None

        # Test single layer of grouping
        x = ["a", "a", "b"]
        y = self.rs.randn(3)
        kws.update(x=x, y=y)
        p = cat._ViolinPlotter(**kws)

        _, ax = plt.subplots()
        p.draw_violins(ax)
        nt.assert_equal(len(ax.collections), 1)
        nt.assert_equal(len(ax.lines), 1)
        plt.close("all")

        # Test nested hue grouping
        x = ["b"] * 4 + ["a"] * 3
        y = self.rs.randn(7)
        h = (["m", "n"] * 4)[:-1]
        kws.update(x=x, y=y, hue=h)
        p = cat._ViolinPlotter(**kws)

        _, ax = plt.subplots()
        p.draw_violins(ax)
        nt.assert_equal(len(ax.collections), 3)
        nt.assert_equal(len(ax.lines), 1)
        plt.close("all")

        # Test nested hue grouping with split
        kws["split"] = True
        p = cat._ViolinPlotter(**kws)

        _, ax = plt.subplots()
        p.draw_violins(ax)
        nt.assert_equal(len(ax.collections), 3)
        nt.assert_equal(len(ax.lines), 1)
        plt.close("all")

    def test_violinplots(self):

        # Smoke test the high level violinplot options

        cat.violinplot("y", data=self.df)
        plt.close("all")

        cat.violinplot(y="y", data=self.df)
        plt.close("all")

        cat.violinplot("g", "y", data=self.df)
        plt.close("all")

        cat.violinplot("y", "g", data=self.df, orient="h")
        plt.close("all")

        cat.violinplot("g", "y", "h", data=self.df)
        plt.close("all")

        cat.violinplot("g", "y", "h", order=list("nabc"), data=self.df)
        plt.close("all")

        cat.violinplot("g", "y", "h", hue_order=list("omn"), data=self.df)
        plt.close("all")

        cat.violinplot("y", "g", "h", data=self.df, orient="h")
        plt.close("all")

        for inner in ["box", "quart", "point", "stick", None]:
            cat.violinplot("g", "y", data=self.df, inner=inner)
            plt.close("all")

            cat.violinplot("g", "y", "h", data=self.df, inner=inner)
            plt.close("all")

            cat.violinplot("g", "y", "h", data=self.df,
                           inner=inner, split=True)
            plt.close("all")


class TestCategoricalScatterPlotter(CategoricalFixture):

    def test_group_point_colors(self):

        p = cat._CategoricalScatterPlotter()

        p.establish_variables(x="g", y="y", data=self.df)
        p.establish_colors(None, "deep", 1)

        point_colors = p.point_colors
        nt.assert_equal(len(point_colors), self.g.unique().size)
        deep_colors = palettes.color_palette("deep", self.g.unique().size)

        for i, group_colors in enumerate(point_colors):
            nt.assert_equal(tuple(deep_colors[i]), tuple(group_colors[0]))
            for channel in group_colors.T:
                nt.assert_equals(np.unique(channel).size, 1)

    def test_hue_point_colors(self):

        p = cat._CategoricalScatterPlotter()

        hue_order = self.h.unique().tolist()
        p.establish_variables(x="g", y="y", hue="h",
                              hue_order=hue_order, data=self.df)
        p.establish_colors(None, "deep", 1)

        point_colors = p.point_colors
        nt.assert_equal(len(point_colors), self.g.unique().size)
        deep_colors = palettes.color_palette("deep", len(hue_order))

        for i, group_colors in enumerate(point_colors):
            for j, point_color in enumerate(group_colors):
                hue_level = p.plot_hues[i][j]
                nt.assert_equal(tuple(point_color),
                                deep_colors[hue_order.index(hue_level)])

    def test_scatterplot_legend(self):

        p = cat._CategoricalScatterPlotter()

        hue_order = ["m", "n"]
        p.establish_variables(x="g", y="y", hue="h",
                              hue_order=hue_order, data=self.df)
        p.establish_colors(None, "deep", 1)
        deep_colors = palettes.color_palette("deep", self.h.unique().size)

        f, ax = plt.subplots()
        p.add_legend_data(ax)
        leg = ax.legend()

        for i, t in enumerate(leg.get_texts()):
            nt.assert_equal(t.get_text(), hue_order[i])

        for i, h in enumerate(leg.legendHandles):
            rgb = h.get_facecolor()[0, :3]
            nt.assert_equal(tuple(rgb), tuple(deep_colors[i]))


class TestStripPlotter(CategoricalFixture):

    def test_stripplot_vertical(self):

        pal = palettes.color_palette()

        ax = cat.stripplot("g", "y", data=self.df)
        for i, (_, vals) in enumerate(self.y.groupby(self.g)):

            x, y = ax.collections[i].get_offsets().T

            npt.assert_array_equal(x, np.ones(len(x)) * i)
            npt.assert_array_equal(y, vals)

            npt.assert_equal(ax.collections[i].get_facecolors()[0, :3], pal[i])

    @skipif(not pandas_has_categoricals)
    def test_stripplot_horiztonal(self):

        df = self.df.copy()
        df.g = df.g.astype("category")

        ax = cat.stripplot("y", "g", data=df)
        for i, (_, vals) in enumerate(self.y.groupby(self.g)):

            x, y = ax.collections[i].get_offsets().T

            npt.assert_array_equal(x, vals)
            npt.assert_array_equal(y, np.ones(len(x)) * i)

    def test_stripplot_jitter(self):

        pal = palettes.color_palette()

        ax = cat.stripplot("g", "y", data=self.df, jitter=True)
        for i, (_, vals) in enumerate(self.y.groupby(self.g)):

            x, y = ax.collections[i].get_offsets().T

            npt.assert_array_less(np.ones(len(x)) * i - .1, x)
            npt.assert_array_less(x, np.ones(len(x)) * i + .1)
            npt.assert_array_equal(y, vals)

            npt.assert_equal(ax.collections[i].get_facecolors()[0, :3], pal[i])

    def test_dodge_nested_stripplot_vertical(self):

        pal = palettes.color_palette()

        ax = cat.stripplot("g", "y", "h", data=self.df, dodge=True)
        for i, (_, group_vals) in enumerate(self.y.groupby(self.g)):
            for j, (_, vals) in enumerate(group_vals.groupby(self.h)):

                x, y = ax.collections[i * 2 + j].get_offsets().T

                npt.assert_array_equal(x, np.ones(len(x)) * i + [-.2, .2][j])
                npt.assert_array_equal(y, vals)

                fc = ax.collections[i * 2 + j].get_facecolors()[0, :3]
                npt.assert_equal(fc, pal[j])

    @skipif(not pandas_has_categoricals)
    def test_dodge_nested_stripplot_horizontal(self):

        df = self.df.copy()
        df.g = df.g.astype("category")

        ax = cat.stripplot("y", "g", "h", data=df, dodge=True)
        for i, (_, group_vals) in enumerate(self.y.groupby(self.g)):
            for j, (_, vals) in enumerate(group_vals.groupby(self.h)):

                x, y = ax.collections[i * 2 + j].get_offsets().T

                npt.assert_array_equal(x, vals)
                npt.assert_array_equal(y, np.ones(len(x)) * i + [-.2, .2][j])

    def test_nested_stripplot_vertical(self):

        # Test a simple vertical strip plot
        ax = cat.stripplot("g", "y", "h", data=self.df, dodge=False)
        for i, (_, group_vals) in enumerate(self.y.groupby(self.g)):

            x, y = ax.collections[i].get_offsets().T

            npt.assert_array_equal(x, np.ones(len(x)) * i)
            npt.assert_array_equal(y, group_vals)

    @skipif(not pandas_has_categoricals)
    def test_nested_stripplot_horizontal(self):

        df = self.df.copy()
        df.g = df.g.astype("category")

        ax = cat.stripplot("y", "g", "h", data=df, dodge=False)
        for i, (_, group_vals) in enumerate(self.y.groupby(self.g)):

            x, y = ax.collections[i].get_offsets().T

            npt.assert_array_equal(x, group_vals)
            npt.assert_array_equal(y, np.ones(len(x)) * i)

    def test_three_strip_points(self):

        x = np.arange(3)
        ax = cat.stripplot(x=x)
        facecolors = ax.collections[0].get_facecolor()
        nt.assert_equal(facecolors.shape, (3, 4))
        npt.assert_array_equal(facecolors[0], facecolors[1])


class TestSwarmPlotter(CategoricalFixture):

    default_kws = dict(x=None, y=None, hue=None, data=None,
                       order=None, hue_order=None, dodge=False,
                       orient=None, color=None, palette=None)

    def test_could_overlap(self):

        p = cat._SwarmPlotter(**self.default_kws)
        neighbors = p.could_overlap((1, 1), [(0, 0), (1, .5), (.5, .5)], 1)
        npt.assert_array_equal(neighbors, [(1, .5), (.5, .5)])

    def test_position_candidates(self):

        p = cat._SwarmPlotter(**self.default_kws)
        xy_i = (0, 1)
        neighbors = [(0, 1), (0, 1.5)]
        candidates = p.position_candidates(xy_i, neighbors, 1)
        dx1 = 1.05
        dx2 = np.sqrt(1 - .5 ** 2) * 1.05
        npt.assert_array_equal(candidates,
                               [(0, 1), (-dx1, 1), (dx1, 1),
                                (dx2, 1), (-dx2, 1)])

    def test_find_first_non_overlapping_candidate(self):

        p = cat._SwarmPlotter(**self.default_kws)
        candidates = [(.5, 1), (1, 1), (1.5, 1)]
        neighbors = np.array([(0, 1)])

        first = p.first_non_overlapping_candidate(candidates, neighbors, 1)
        npt.assert_array_equal(first, (1, 1))

    def test_beeswarm(self):

        p = cat._SwarmPlotter(**self.default_kws)
        d = self.y.diff().mean() * 1.5
        x = np.zeros(self.y.size)
        y = np.sort(self.y)
        orig_xy = np.c_[x, y]
        swarm = p.beeswarm(orig_xy, d)
        dmat = spatial.distance.cdist(swarm, swarm)
        triu = dmat[np.triu_indices_from(dmat, 1)]
        npt.assert_array_less(d, triu)
        npt.assert_array_equal(y, swarm[:, 1])

    def test_add_gutters(self):

        p = cat._SwarmPlotter(**self.default_kws)
        points = np.array([0, -1, .4, .8])
        points = p.add_gutters(points, 0, 1)
        npt.assert_array_equal(points,
                               np.array([0, -.5, .4, .5]))

    def test_swarmplot_vertical(self):

        pal = palettes.color_palette()

        ax = cat.swarmplot("g", "y", data=self.df)
        for i, (_, vals) in enumerate(self.y.groupby(self.g)):

            x, y = ax.collections[i].get_offsets().T
            npt.assert_array_almost_equal(y, np.sort(vals))

            fc = ax.collections[i].get_facecolors()[0, :3]
            npt.assert_equal(fc, pal[i])

    def test_swarmplot_horizontal(self):

        pal = palettes.color_palette()

        ax = cat.swarmplot("y", "g", data=self.df, orient="h")
        for i, (_, vals) in enumerate(self.y.groupby(self.g)):

            x, y = ax.collections[i].get_offsets().T
            npt.assert_array_almost_equal(x, np.sort(vals))

            fc = ax.collections[i].get_facecolors()[0, :3]
            npt.assert_equal(fc, pal[i])

    def test_dodge_nested_swarmplot_vetical(self):

        pal = palettes.color_palette()

        ax = cat.swarmplot("g", "y", "h", data=self.df, dodge=True)
        for i, (_, group_vals) in enumerate(self.y.groupby(self.g)):
            for j, (_, vals) in enumerate(group_vals.groupby(self.h)):

                x, y = ax.collections[i * 2 + j].get_offsets().T
                npt.assert_array_almost_equal(y, np.sort(vals))

                fc = ax.collections[i * 2 + j].get_facecolors()[0, :3]
                npt.assert_equal(fc, pal[j])

    def test_dodge_nested_swarmplot_horizontal(self):

        pal = palettes.color_palette()

        ax = cat.swarmplot("y", "g", "h", data=self.df, orient="h", dodge=True)
        for i, (_, group_vals) in enumerate(self.y.groupby(self.g)):
            for j, (_, vals) in enumerate(group_vals.groupby(self.h)):

                x, y = ax.collections[i * 2 + j].get_offsets().T
                npt.assert_array_almost_equal(x, np.sort(vals))

                fc = ax.collections[i * 2 + j].get_facecolors()[0, :3]
                npt.assert_equal(fc, pal[j])

    def test_nested_swarmplot_vertical(self):

        ax = cat.swarmplot("g", "y", "h", data=self.df)

        pal = palettes.color_palette()
        hue_names = self.h.unique().tolist()
        grouped_hues = list(self.h.groupby(self.g))

        for i, (_, vals) in enumerate(self.y.groupby(self.g)):

            points = ax.collections[i]
            x, y = points.get_offsets().T
            sorter = np.argsort(vals)
            npt.assert_array_almost_equal(y, vals.iloc[sorter])

            _, hue_vals = grouped_hues[i]
            for hue, fc in zip(hue_vals.values[sorter.values],
                               points.get_facecolors()):

                npt.assert_equal(fc[:3], pal[hue_names.index(hue)])

    def test_nested_swarmplot_horizontal(self):

        ax = cat.swarmplot("y", "g", "h", data=self.df, orient="h")

        pal = palettes.color_palette()
        hue_names = self.h.unique().tolist()
        grouped_hues = list(self.h.groupby(self.g))

        for i, (_, vals) in enumerate(self.y.groupby(self.g)):

            points = ax.collections[i]
            x, y = points.get_offsets().T
            sorter = np.argsort(vals)
            npt.assert_array_almost_equal(x, vals.iloc[sorter])

            _, hue_vals = grouped_hues[i]
            for hue, fc in zip(hue_vals.values[sorter.values],
                               points.get_facecolors()):

                npt.assert_equal(fc[:3], pal[hue_names.index(hue)])


class TestBarPlotter(CategoricalFixture):

    default_kws = dict(x=None, y=None, hue=None, data=None,
                       estimator=np.mean, ci=95, n_boot=100, units=None,
                       order=None, hue_order=None,
                       orient=None, color=None, palette=None,
                       saturation=.75, errcolor=".26", errwidth=None,
                       capsize=None, dodge=True)

    def test_nested_width(self):

        kws = self.default_kws.copy()

        p = cat._BarPlotter(**kws)
        p.establish_variables("g", "y", "h", data=self.df)
        nt.assert_equal(p.nested_width, .8 / 2)

        p = cat._BarPlotter(**kws)
        p.establish_variables("h", "y", "g", data=self.df)
        nt.assert_equal(p.nested_width, .8 / 3)

        kws["dodge"] = False
        p = cat._BarPlotter(**kws)
        p.establish_variables("h", "y", "g", data=self.df)
        nt.assert_equal(p.nested_width, .8)

    def test_draw_vertical_bars(self):

        kws = self.default_kws.copy()
        kws.update(x="g", y="y", data=self.df)
        p = cat._BarPlotter(**kws)

        f, ax = plt.subplots()
        p.draw_bars(ax, {})

        nt.assert_equal(len(ax.patches), len(p.plot_data))
        nt.assert_equal(len(ax.lines), len(p.plot_data))

        for bar, color in zip(ax.patches, p.colors):
            nt.assert_equal(bar.get_facecolor()[:-1], color)

        positions = np.arange(len(p.plot_data)) - p.width / 2
        for bar, pos, stat in zip(ax.patches, positions, p.statistic):
            nt.assert_equal(bar.get_x(), pos)
            nt.assert_equal(bar.get_width(), p.width)
            if mpl.__version__ >= mpl_barplot_change:
                nt.assert_equal(bar.get_y(), 0)
                nt.assert_equal(bar.get_height(), stat)
            else:
                nt.assert_equal(bar.get_y(), min(0, stat))
                nt.assert_equal(bar.get_height(), abs(stat))

    def test_draw_horizontal_bars(self):

        kws = self.default_kws.copy()
        kws.update(x="y", y="g", orient="h", data=self.df)
        p = cat._BarPlotter(**kws)

        f, ax = plt.subplots()
        p.draw_bars(ax, {})

        nt.assert_equal(len(ax.patches), len(p.plot_data))
        nt.assert_equal(len(ax.lines), len(p.plot_data))

        for bar, color in zip(ax.patches, p.colors):
            nt.assert_equal(bar.get_facecolor()[:-1], color)

        positions = np.arange(len(p.plot_data)) - p.width / 2
        for bar, pos, stat in zip(ax.patches, positions, p.statistic):
            nt.assert_equal(bar.get_y(), pos)
            nt.assert_equal(bar.get_height(), p.width)
            if mpl.__version__ >= mpl_barplot_change:
                nt.assert_equal(bar.get_x(), 0)
                nt.assert_equal(bar.get_width(), stat)
            else:
                nt.assert_equal(bar.get_x(), min(0, stat))
                nt.assert_equal(bar.get_width(), abs(stat))

    def test_draw_nested_vertical_bars(self):

        kws = self.default_kws.copy()
        kws.update(x="g", y="y", hue="h", data=self.df)
        p = cat._BarPlotter(**kws)

        f, ax = plt.subplots()
        p.draw_bars(ax, {})

        n_groups, n_hues = len(p.plot_data), len(p.hue_names)
        nt.assert_equal(len(ax.patches), n_groups * n_hues)
        nt.assert_equal(len(ax.lines), n_groups * n_hues)

        for bar in ax.patches[:n_groups]:
            nt.assert_equal(bar.get_facecolor()[:-1], p.colors[0])
        for bar in ax.patches[n_groups:]:
            nt.assert_equal(bar.get_facecolor()[:-1], p.colors[1])

        positions = np.arange(len(p.plot_data))
        for bar, pos in zip(ax.patches[:n_groups], positions):
            nt.assert_almost_equal(bar.get_x(), pos - p.width / 2)
            nt.assert_almost_equal(bar.get_width(), p.nested_width)

        for bar, stat in zip(ax.patches, p.statistic.T.flat):
            if LooseVersion(mpl.__version__) >= mpl_barplot_change:
                nt.assert_almost_equal(bar.get_y(), 0)
                nt.assert_almost_equal(bar.get_height(), stat)
            else:
                nt.assert_almost_equal(bar.get_y(), min(0, stat))
                nt.assert_almost_equal(bar.get_height(), abs(stat))

    def test_draw_nested_horizontal_bars(self):

        kws = self.default_kws.copy()
        kws.update(x="y", y="g", hue="h", orient="h", data=self.df)
        p = cat._BarPlotter(**kws)

        f, ax = plt.subplots()
        p.draw_bars(ax, {})

        n_groups, n_hues = len(p.plot_data), len(p.hue_names)
        nt.assert_equal(len(ax.patches), n_groups * n_hues)
        nt.assert_equal(len(ax.lines), n_groups * n_hues)

        for bar in ax.patches[:n_groups]:
            nt.assert_equal(bar.get_facecolor()[:-1], p.colors[0])
        for bar in ax.patches[n_groups:]:
            nt.assert_equal(bar.get_facecolor()[:-1], p.colors[1])

        positions = np.arange(len(p.plot_data))
        for bar, pos in zip(ax.patches[:n_groups], positions):
            nt.assert_almost_equal(bar.get_y(), pos - p.width / 2)
            nt.assert_almost_equal(bar.get_height(), p.nested_width)

        for bar, stat in zip(ax.patches, p.statistic.T.flat):
            if LooseVersion(mpl.__version__) >= mpl_barplot_change:
                nt.assert_almost_equal(bar.get_x(), 0)
                nt.assert_almost_equal(bar.get_width(), stat)
            else:
                nt.assert_almost_equal(bar.get_x(), min(0, stat))
                nt.assert_almost_equal(bar.get_width(), abs(stat))

    def test_draw_missing_bars(self):

        kws = self.default_kws.copy()

        order = list("abcd")
        kws.update(x="g", y="y", order=order, data=self.df)
        p = cat._BarPlotter(**kws)

        f, ax = plt.subplots()
        p.draw_bars(ax, {})

        nt.assert_equal(len(ax.patches), len(order))
        nt.assert_equal(len(ax.lines), len(order))

        plt.close("all")

        hue_order = list("mno")
        kws.update(x="g", y="y", hue="h", hue_order=hue_order, data=self.df)
        p = cat._BarPlotter(**kws)

        f, ax = plt.subplots()
        p.draw_bars(ax, {})

        nt.assert_equal(len(ax.patches), len(p.plot_data) * len(hue_order))
        nt.assert_equal(len(ax.lines),  len(p.plot_data) * len(hue_order))

        plt.close("all")

    def test_barplot_colors(self):

        # Test unnested palette colors
        kws = self.default_kws.copy()
        kws.update(x="g", y="y", data=self.df,
                   saturation=1, palette="muted")
        p = cat._BarPlotter(**kws)

        f, ax = plt.subplots()
        p.draw_bars(ax, {})

        palette = palettes.color_palette("muted", len(self.g.unique()))
        for patch, pal_color in zip(ax.patches, palette):
            nt.assert_equal(patch.get_facecolor()[:-1], pal_color)

        plt.close("all")

        # Test single color
        color = (.2, .2, .3, 1)
        kws = self.default_kws.copy()
        kws.update(x="g", y="y", data=self.df,
                   saturation=1, color=color)
        p = cat._BarPlotter(**kws)

        f, ax = plt.subplots()
        p.draw_bars(ax, {})

        for patch in ax.patches:
            nt.assert_equal(patch.get_facecolor(), color)

        plt.close("all")

        # Test nested palette colors
        kws = self.default_kws.copy()
        kws.update(x="g", y="y", hue="h", data=self.df,
                   saturation=1, palette="Set2")
        p = cat._BarPlotter(**kws)

        f, ax = plt.subplots()
        p.draw_bars(ax, {})

        palette = palettes.color_palette("Set2", len(self.h.unique()))
        for patch in ax.patches[:len(self.g.unique())]:
            nt.assert_equal(patch.get_facecolor()[:-1], palette[0])
        for patch in ax.patches[len(self.g.unique()):]:
            nt.assert_equal(patch.get_facecolor()[:-1], palette[1])

        plt.close("all")

    def test_simple_barplots(self):

        ax = cat.barplot("g", "y", data=self.df)
        nt.assert_equal(len(ax.patches), len(self.g.unique()))
        nt.assert_equal(ax.get_xlabel(), "g")
        nt.assert_equal(ax.get_ylabel(), "y")
        plt.close("all")

        ax = cat.barplot("y", "g", orient="h", data=self.df)
        nt.assert_equal(len(ax.patches), len(self.g.unique()))
        nt.assert_equal(ax.get_xlabel(), "y")
        nt.assert_equal(ax.get_ylabel(), "g")
        plt.close("all")

        ax = cat.barplot("g", "y", "h", data=self.df)
        nt.assert_equal(len(ax.patches),
                        len(self.g.unique()) * len(self.h.unique()))
        nt.assert_equal(ax.get_xlabel(), "g")
        nt.assert_equal(ax.get_ylabel(), "y")
        plt.close("all")

        ax = cat.barplot("y", "g", "h", orient="h", data=self.df)
        nt.assert_equal(len(ax.patches),
                        len(self.g.unique()) * len(self.h.unique()))
        nt.assert_equal(ax.get_xlabel(), "y")
        nt.assert_equal(ax.get_ylabel(), "g")
        plt.close("all")


class TestPointPlotter(CategoricalFixture):

    default_kws = dict(x=None, y=None, hue=None, data=None,
                       estimator=np.mean, ci=95, n_boot=100, units=None,
                       order=None, hue_order=None,
                       markers="o", linestyles="-", dodge=0,
                       join=True, scale=1,
                       orient=None, color=None, palette=None)

    def test_different_defualt_colors(self):

        kws = self.default_kws.copy()
        kws.update(dict(x="g", y="y", data=self.df))
        p = cat._PointPlotter(**kws)
        color = palettes.color_palette()[0]
        npt.assert_array_equal(p.colors, [color, color, color])

    def test_hue_offsets(self):

        kws = self.default_kws.copy()
        kws.update(dict(x="g", y="y", hue="h", data=self.df))

        p = cat._PointPlotter(**kws)
        npt.assert_array_equal(p.hue_offsets, [0, 0])

        kws.update(dict(dodge=.5))

        p = cat._PointPlotter(**kws)
        npt.assert_array_equal(p.hue_offsets, [-.25, .25])

        kws.update(dict(x="h", hue="g", dodge=0))

        p = cat._PointPlotter(**kws)
        npt.assert_array_equal(p.hue_offsets, [0, 0, 0])

        kws.update(dict(dodge=.3))

        p = cat._PointPlotter(**kws)
        npt.assert_array_equal(p.hue_offsets, [-.15, 0, .15])

    def test_draw_vertical_points(self):

        kws = self.default_kws.copy()
        kws.update(x="g", y="y", data=self.df)
        p = cat._PointPlotter(**kws)

        f, ax = plt.subplots()
        p.draw_points(ax)

        nt.assert_equal(len(ax.collections), 1)
        nt.assert_equal(len(ax.lines), len(p.plot_data) + 1)
        points = ax.collections[0]
        nt.assert_equal(len(points.get_offsets()), len(p.plot_data))

        x, y = points.get_offsets().T
        npt.assert_array_equal(x, np.arange(len(p.plot_data)))
        npt.assert_array_equal(y, p.statistic)

        for got_color, want_color in zip(points.get_facecolors(),
                                         p.colors):
            npt.assert_array_equal(got_color[:-1], want_color)

    def test_draw_horizontal_points(self):

        kws = self.default_kws.copy()
        kws.update(x="y", y="g", orient="h", data=self.df)
        p = cat._PointPlotter(**kws)

        f, ax = plt.subplots()
        p.draw_points(ax)

        nt.assert_equal(len(ax.collections), 1)
        nt.assert_equal(len(ax.lines), len(p.plot_data) + 1)
        points = ax.collections[0]
        nt.assert_equal(len(points.get_offsets()), len(p.plot_data))

        x, y = points.get_offsets().T
        npt.assert_array_equal(x, p.statistic)
        npt.assert_array_equal(y, np.arange(len(p.plot_data)))

        for got_color, want_color in zip(points.get_facecolors(),
                                         p.colors):
            npt.assert_array_equal(got_color[:-1], want_color)

    def test_draw_vertical_nested_points(self):

        kws = self.default_kws.copy()
        kws.update(x="g", y="y", hue="h", data=self.df)
        p = cat._PointPlotter(**kws)

        f, ax = plt.subplots()
        p.draw_points(ax)

        nt.assert_equal(len(ax.collections), 2)
        nt.assert_equal(len(ax.lines),
                        len(p.plot_data) * len(p.hue_names) + len(p.hue_names))

        for points, stats, color in zip(ax.collections,
                                        p.statistic.T,
                                        p.colors):

            nt.assert_equal(len(points.get_offsets()), len(p.plot_data))

            x, y = points.get_offsets().T
            npt.assert_array_equal(x, np.arange(len(p.plot_data)))
            npt.assert_array_equal(y, stats)

            for got_color in points.get_facecolors():
                npt.assert_array_equal(got_color[:-1], color)

    def test_draw_horizontal_nested_points(self):

        kws = self.default_kws.copy()
        kws.update(x="y", y="g", hue="h", orient="h", data=self.df)
        p = cat._PointPlotter(**kws)

        f, ax = plt.subplots()
        p.draw_points(ax)

        nt.assert_equal(len(ax.collections), 2)
        nt.assert_equal(len(ax.lines),
                        len(p.plot_data) * len(p.hue_names) + len(p.hue_names))

        for points, stats, color in zip(ax.collections,
                                        p.statistic.T,
                                        p.colors):

            nt.assert_equal(len(points.get_offsets()), len(p.plot_data))

            x, y = points.get_offsets().T
            npt.assert_array_equal(x, stats)
            npt.assert_array_equal(y, np.arange(len(p.plot_data)))

            for got_color in points.get_facecolors():
                npt.assert_array_equal(got_color[:-1], color)

    def test_pointplot_colors(self):

        # Test a single-color unnested plot
        color = (.2, .2, .3, 1)
        kws = self.default_kws.copy()
        kws.update(x="g", y="y", data=self.df, color=color)
        p = cat._PointPlotter(**kws)

        f, ax = plt.subplots()
        p.draw_points(ax)

        for line in ax.lines:
            nt.assert_equal(line.get_color(), color[:-1])

        for got_color in ax.collections[0].get_facecolors():
            npt.assert_array_equal(rgb2hex(got_color), rgb2hex(color))

        plt.close("all")

        # Test a multi-color unnested plot
        palette = palettes.color_palette("Set1", 3)
        kws.update(x="g", y="y", data=self.df, palette="Set1")
        p = cat._PointPlotter(**kws)

        nt.assert_true(not p.join)

        f, ax = plt.subplots()
        p.draw_points(ax)

        for line, pal_color in zip(ax.lines, palette):
            npt.assert_array_equal(line.get_color(), pal_color)

        for point_color, pal_color in zip(ax.collections[0].get_facecolors(),
                                          palette):
            npt.assert_array_equal(rgb2hex(point_color), rgb2hex(pal_color))

        plt.close("all")

        # Test a multi-colored nested plot
        palette = palettes.color_palette("dark", 2)
        kws.update(x="g", y="y", hue="h", data=self.df, palette="dark")
        p = cat._PointPlotter(**kws)

        f, ax = plt.subplots()
        p.draw_points(ax)

        for line in ax.lines[:(len(p.plot_data) + 1)]:
            nt.assert_equal(line.get_color(), palette[0])
        for line in ax.lines[(len(p.plot_data) + 1):]:
            nt.assert_equal(line.get_color(), palette[1])

        for i, pal_color in enumerate(palette):
            for point_color in ax.collections[i].get_facecolors():
                npt.assert_array_equal(point_color[:-1], pal_color)

        plt.close("all")

    def test_simple_pointplots(self):

        ax = cat.pointplot("g", "y", data=self.df)
        nt.assert_equal(len(ax.collections), 1)
        nt.assert_equal(len(ax.lines), len(self.g.unique()) + 1)
        nt.assert_equal(ax.get_xlabel(), "g")
        nt.assert_equal(ax.get_ylabel(), "y")
        plt.close("all")

        ax = cat.pointplot("y", "g", orient="h", data=self.df)
        nt.assert_equal(len(ax.collections), 1)
        nt.assert_equal(len(ax.lines), len(self.g.unique()) + 1)
        nt.assert_equal(ax.get_xlabel(), "y")
        nt.assert_equal(ax.get_ylabel(), "g")
        plt.close("all")

        ax = cat.pointplot("g", "y", "h", data=self.df)
        nt.assert_equal(len(ax.collections), len(self.h.unique()))
        nt.assert_equal(len(ax.lines),
                        (len(self.g.unique()) *
                         len(self.h.unique()) +
                         len(self.h.unique())))
        nt.assert_equal(ax.get_xlabel(), "g")
        nt.assert_equal(ax.get_ylabel(), "y")
        plt.close("all")

        ax = cat.pointplot("y", "g", "h", orient="h", data=self.df)
        nt.assert_equal(len(ax.collections), len(self.h.unique()))
        nt.assert_equal(len(ax.lines),
                        (len(self.g.unique()) *
                         len(self.h.unique()) +
                         len(self.h.unique())))
        nt.assert_equal(ax.get_xlabel(), "y")
        nt.assert_equal(ax.get_ylabel(), "g")
        plt.close("all")


class TestCountPlot(CategoricalFixture):

    def test_plot_elements(self):

        ax = cat.countplot("g", data=self.df)
        nt.assert_equal(len(ax.patches), self.g.unique().size)
        for p in ax.patches:
            nt.assert_equal(p.get_y(), 0)
            nt.assert_equal(p.get_height(),
                            self.g.size / self.g.unique().size)
        plt.close("all")

        ax = cat.countplot(y="g", data=self.df)
        nt.assert_equal(len(ax.patches), self.g.unique().size)
        for p in ax.patches:
            nt.assert_equal(p.get_x(), 0)
            nt.assert_equal(p.get_width(),
                            self.g.size / self.g.unique().size)
        plt.close("all")

        ax = cat.countplot("g", hue="h", data=self.df)
        nt.assert_equal(len(ax.patches),
                        self.g.unique().size * self.h.unique().size)
        plt.close("all")

        ax = cat.countplot(y="g", hue="h", data=self.df)
        nt.assert_equal(len(ax.patches),
                        self.g.unique().size * self.h.unique().size)
        plt.close("all")

    def test_input_error(self):

        with nt.assert_raises(TypeError):
            cat.countplot()

        with nt.assert_raises(TypeError):
            cat.countplot(x="g", y="h", data=self.df)


class TestFactorPlot(CategoricalFixture):

    def test_facet_organization(self):

        g = cat.factorplot("g", "y", data=self.df)
        nt.assert_equal(g.axes.shape, (1, 1))

        g = cat.factorplot("g", "y", col="h", data=self.df)
        nt.assert_equal(g.axes.shape, (1, 2))

        g = cat.factorplot("g", "y", row="h", data=self.df)
        nt.assert_equal(g.axes.shape, (2, 1))

        g = cat.factorplot("g", "y", col="u", row="h", data=self.df)
        nt.assert_equal(g.axes.shape, (2, 3))

    def test_plot_elements(self):

        g = cat.factorplot("g", "y", data=self.df)
        nt.assert_equal(len(g.ax.collections), 1)
        want_lines = self.g.unique().size + 1
        nt.assert_equal(len(g.ax.lines), want_lines)

        g = cat.factorplot("g", "y", "h", data=self.df)
        want_collections = self.h.unique().size
        nt.assert_equal(len(g.ax.collections), want_collections)
        want_lines = (self.g.unique().size + 1) * self.h.unique().size
        nt.assert_equal(len(g.ax.lines), want_lines)

        g = cat.factorplot("g", "y", data=self.df, kind="bar")
        want_elements = self.g.unique().size
        nt.assert_equal(len(g.ax.patches), want_elements)
        nt.assert_equal(len(g.ax.lines), want_elements)

        g = cat.factorplot("g", "y", "h", data=self.df, kind="bar")
        want_elements = self.g.unique().size * self.h.unique().size
        nt.assert_equal(len(g.ax.patches), want_elements)
        nt.assert_equal(len(g.ax.lines), want_elements)

        g = cat.factorplot("g", data=self.df, kind="count")
        want_elements = self.g.unique().size
        nt.assert_equal(len(g.ax.patches), want_elements)
        nt.assert_equal(len(g.ax.lines), 0)

        g = cat.factorplot("g", hue="h", data=self.df, kind="count")
        want_elements = self.g.unique().size * self.h.unique().size
        nt.assert_equal(len(g.ax.patches), want_elements)
        nt.assert_equal(len(g.ax.lines), 0)

        g = cat.factorplot("g", "y", data=self.df, kind="box")
        want_artists = self.g.unique().size
        nt.assert_equal(len(g.ax.artists), want_artists)

        g = cat.factorplot("g", "y", "h", data=self.df, kind="box")
        want_artists = self.g.unique().size * self.h.unique().size
        nt.assert_equal(len(g.ax.artists), want_artists)

        g = cat.factorplot("g", "y", data=self.df,
                           kind="violin", inner=None)
        want_elements = self.g.unique().size
        nt.assert_equal(len(g.ax.collections), want_elements)

        g = cat.factorplot("g", "y", "h", data=self.df,
                           kind="violin", inner=None)
        want_elements = self.g.unique().size * self.h.unique().size
        nt.assert_equal(len(g.ax.collections), want_elements)

        g = cat.factorplot("g", "y", data=self.df, kind="strip")
        want_elements = self.g.unique().size
        nt.assert_equal(len(g.ax.collections), want_elements)

        g = cat.factorplot("g", "y", "h", data=self.df, kind="strip")
        want_elements = self.g.unique().size + self.h.unique().size
        nt.assert_equal(len(g.ax.collections), want_elements)

    def test_bad_plot_kind_error(self):

        with nt.assert_raises(ValueError):
            cat.factorplot("g", "y", data=self.df, kind="not_a_kind")

    def test_count_x_and_y(self):

        with nt.assert_raises(ValueError):
            cat.factorplot("g", "y", data=self.df, kind="count")

    def test_plot_colors(self):

        ax = cat.barplot("g", "y", data=self.df)
        g = cat.factorplot("g", "y", data=self.df, kind="bar")
        for p1, p2 in zip(ax.patches, g.ax.patches):
            nt.assert_equal(p1.get_facecolor(), p2.get_facecolor())
        plt.close("all")

        ax = cat.barplot("g", "y", data=self.df, color="purple")
        g = cat.factorplot("g", "y", data=self.df,
                           kind="bar", color="purple")
        for p1, p2 in zip(ax.patches, g.ax.patches):
            nt.assert_equal(p1.get_facecolor(), p2.get_facecolor())
        plt.close("all")

        ax = cat.barplot("g", "y", data=self.df, palette="Set2")
        g = cat.factorplot("g", "y", data=self.df,
                           kind="bar", palette="Set2")
        for p1, p2 in zip(ax.patches, g.ax.patches):
            nt.assert_equal(p1.get_facecolor(), p2.get_facecolor())
        plt.close("all")

        ax = cat.pointplot("g", "y", data=self.df)
        g = cat.factorplot("g", "y", data=self.df)
        for l1, l2 in zip(ax.lines, g.ax.lines):
            nt.assert_equal(l1.get_color(), l2.get_color())
        plt.close("all")

        ax = cat.pointplot("g", "y", data=self.df, color="purple")
        g = cat.factorplot("g", "y", data=self.df, color="purple")
        for l1, l2 in zip(ax.lines, g.ax.lines):
            nt.assert_equal(l1.get_color(), l2.get_color())
        plt.close("all")

        ax = cat.pointplot("g", "y", data=self.df, palette="Set2")
        g = cat.factorplot("g", "y", data=self.df, palette="Set2")
        for l1, l2 in zip(ax.lines, g.ax.lines):
            nt.assert_equal(l1.get_color(), l2.get_color())
        plt.close("all")


class TestLVPlotter(CategoricalFixture):

    def setUp(self):
        def edge_calc(n, data):
            q = np.asanyarray([0.5**n, 1 - 0.5**n]) * 100
            q = list(np.unique(q))
            return np.percentile(data, q)

        self.ispatch = lambda c: isinstance(c, mpl.collections.PatchCollection)

        self.default_kws = dict(x=None, y=None, hue=None, data=None,
                                order=None, hue_order=None,
                                orient=None, color=None, palette=None,
                                saturation=.75, width=.8, dodge=True,
                                k_depth='proportion', linewidth=None,
                                scale='exponential', outlier_prop=None)
        self.linear_data = np.arange(101)
        self.n = len(self.linear_data)
        self.expected_k = int(np.log2(self.n)) - int(np.log2(self.n*0.007)) + 1
        self.expected_edges_l = map(lambda i: edge_calc(i, self.linear_data),
                                    range(self.expected_k + 2, 1, -1))
        self.outlier_data = np.concatenate((np.arange(100), [200]))
        self.expected_edges_o = map(lambda i: edge_calc(i, self.outlier_data),
                                    range(self.expected_k + 2, 1, -1))

    def test_box_ends_finite(self):
        p = cat._LVPlotter(**self.default_kws)
        p.establish_variables("g", "y", data=self.df)
        box_k = np.asarray([[b, k]
                           for b, k in map(p._lv_box_ends, p.plot_data)])
        box_ends = box_k[:, 0]
        k_vals = box_k[:, 1]

        # Check that all the box ends are finite and are within
        # the bounds of the data
        b_e = map(lambda a: np.all(np.isfinite(a)), box_ends)
        npt.assert_equal(np.sum(list(b_e)), len(box_ends))

        def within(t):
            a, d = t
            return ((np.ravel(a) <= d.max()) &
                    (np.ravel(a) >= d.min())).all()

        b_w = map(within, zip(box_ends, p.plot_data))
        npt.assert_equal(np.sum(list(b_w)), len(box_ends))

        k_f = map(lambda k: (k > 0.) & np.isfinite(k), k_vals)
        npt.assert_equal(np.sum(list(k_f)), len(k_vals))

    def test_box_ends_correct(self):
        p = cat._LVPlotter(**self.default_kws)
        calc_edges, calc_k = p._lv_box_ends(self.linear_data)

        npt.assert_equal(list(self.expected_edges_l), calc_edges)

        npt.assert_equal(self.expected_k, calc_k)

    def test_outliers(self):
        p = cat._LVPlotter(**self.default_kws)
        calc_edges, calc_k = p._lv_box_ends(self.outlier_data)

        npt.assert_equal(list(self.expected_edges_o), calc_edges)

        npt.assert_equal(self.expected_k, calc_k)

        out_calc = p._lv_outliers(self.outlier_data, calc_k)
        out_exp = p._lv_outliers(self.outlier_data, self.expected_k)

        npt.assert_equal(out_exp, out_calc)

    def test_hue_offsets(self):

        p = cat._LVPlotter(**self.default_kws)
        p.establish_variables("g", "y", "h", data=self.df)
        npt.assert_array_equal(p.hue_offsets, [-.2, .2])

        kws = self.default_kws.copy()
        kws["width"] = .6
        p = cat._LVPlotter(**kws)
        p.establish_variables("g", "y", "h", data=self.df)
        npt.assert_array_equal(p.hue_offsets, [-.15, .15])

        p = cat._LVPlotter(**kws)
        p.establish_variables("h", "y", "g", data=self.df)
        npt.assert_array_almost_equal(p.hue_offsets, [-.2, 0, .2])

    def test_axes_data(self):

        ax = cat.lvplot("g", "y", data=self.df)
        patches = filter(self.ispatch, ax.collections)
        nt.assert_equal(len(list(patches)), 3)

        plt.close("all")

        ax = cat.lvplot("g", "y", "h", data=self.df)
        patches = filter(self.ispatch, ax.collections)
        nt.assert_equal(len(list(patches)), 6)

        plt.close("all")

    def test_box_colors(self):

        ax = cat.lvplot("g", "y", data=self.df, saturation=1)
        pal = palettes.color_palette(n_colors=3)
        for patch, color in zip(ax.artists, pal):
            nt.assert_equal(patch.get_facecolor()[:3], color)

        plt.close("all")

        ax = cat.lvplot("g", "y", "h", data=self.df, saturation=1)
        pal = palettes.color_palette(n_colors=2)
        for patch, color in zip(ax.artists, pal * 2):
            nt.assert_equal(patch.get_facecolor()[:3], color)

        plt.close("all")

    def test_draw_missing_boxes(self):

        ax = cat.lvplot("g", "y", data=self.df,
                        order=["a", "b", "c", "d"])

        patches = filter(self.ispatch, ax.collections)
        nt.assert_equal(len(list(patches)), 3)
        plt.close("all")

    def test_missing_data(self):

        x = ["a", "a", "b", "b", "c", "c", "d", "d"]
        h = ["x", "y", "x", "y", "x", "y", "x", "y"]
        y = self.rs.randn(8)
        y[-2:] = np.nan

        ax = cat.lvplot(x, y)
        patches = filter(self.ispatch, ax.collections)
        nt.assert_equal(len(ax.lines), 3)

        plt.close("all")

        y[-1] = 0
        ax = cat.lvplot(x, y, h)
        nt.assert_equal(len(ax.lines), 7)

        plt.close("all")

    def test_lvplots(self):

        # Smoke test the high level lvplot options

        cat.lvplot("y", data=self.df)
        plt.close("all")

        cat.lvplot(y="y", data=self.df)
        plt.close("all")

        cat.lvplot("g", "y", data=self.df)
        plt.close("all")

        cat.lvplot("y", "g", data=self.df, orient="h")
        plt.close("all")

        cat.lvplot("g", "y", "h", data=self.df)
        plt.close("all")

        cat.lvplot("g", "y", "h", order=list("nabc"), data=self.df)
        plt.close("all")

        cat.lvplot("g", "y", "h", hue_order=list("omn"), data=self.df)
        plt.close("all")

        cat.lvplot("y", "g", "h", data=self.df, orient="h")
        plt.close("all")

    def test_axes_annotation(self):

        ax = cat.lvplot("g", "y", data=self.df)
        nt.assert_equal(ax.get_xlabel(), "g")
        nt.assert_equal(ax.get_ylabel(), "y")
        nt.assert_equal(ax.get_xlim(), (-.5, 2.5))
        npt.assert_array_equal(ax.get_xticks(), [0, 1, 2])
        npt.assert_array_equal([l.get_text() for l in ax.get_xticklabels()],
                               ["a", "b", "c"])

        plt.close("all")

        ax = cat.lvplot("g", "y", "h", data=self.df)
        nt.assert_equal(ax.get_xlabel(), "g")
        nt.assert_equal(ax.get_ylabel(), "y")
        npt.assert_array_equal(ax.get_xticks(), [0, 1, 2])
        npt.assert_array_equal([l.get_text() for l in ax.get_xticklabels()],
                               ["a", "b", "c"])
        npt.assert_array_equal([l.get_text() for l in ax.legend_.get_texts()],
                               ["m", "n"])

        plt.close("all")

        ax = cat.lvplot("y", "g", data=self.df, orient="h")
        nt.assert_equal(ax.get_xlabel(), "y")
        nt.assert_equal(ax.get_ylabel(), "g")
        nt.assert_equal(ax.get_ylim(), (2.5, -.5))
        npt.assert_array_equal(ax.get_yticks(), [0, 1, 2])
        npt.assert_array_equal([l.get_text() for l in ax.get_yticklabels()],
                               ["a", "b", "c"])

        plt.close("all")
