import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

import nose.tools as nt
import numpy.testing as npt
try:
    import pandas.testing as pdt
except ImportError:
    import pandas.util.testing as pdt
from numpy.testing.decorators import skipif
from nose import SkipTest

try:
    import statsmodels.regression.linear_model as smlm
    _no_statsmodels = False
except ImportError:
    _no_statsmodels = True

from . import PlotTestCase
from .. import regression as lm
from ..palettes import color_palette

rs = np.random.RandomState(0)


class TestLinearPlotter(PlotTestCase):

    rs = np.random.RandomState(77)
    df = pd.DataFrame(dict(x=rs.normal(size=60),
                           d=rs.randint(-2, 3, 60),
                           y=rs.gamma(4, size=60),
                           s=np.tile(list("abcdefghij"), 6)))
    df["z"] = df.y + rs.randn(60)
    df["y_na"] = df.y.copy()
    df.loc[[10, 20, 30], 'y_na'] = np.nan

    def test_establish_variables_from_frame(self):

        p = lm._LinearPlotter()
        p.establish_variables(self.df, x="x", y="y")
        pdt.assert_series_equal(p.x, self.df.x)
        pdt.assert_series_equal(p.y, self.df.y)
        pdt.assert_frame_equal(p.data, self.df)

    def test_establish_variables_from_series(self):

        p = lm._LinearPlotter()
        p.establish_variables(None, x=self.df.x, y=self.df.y)
        pdt.assert_series_equal(p.x, self.df.x)
        pdt.assert_series_equal(p.y, self.df.y)
        nt.assert_is(p.data, None)

    def test_establish_variables_from_array(self):

        p = lm._LinearPlotter()
        p.establish_variables(None,
                              x=self.df.x.values,
                              y=self.df.y.values)
        npt.assert_array_equal(p.x, self.df.x)
        npt.assert_array_equal(p.y, self.df.y)
        nt.assert_is(p.data, None)

    def test_establish_variables_from_mix(self):

        p = lm._LinearPlotter()
        p.establish_variables(self.df, x="x", y=self.df.y)
        pdt.assert_series_equal(p.x, self.df.x)
        pdt.assert_series_equal(p.y, self.df.y)
        pdt.assert_frame_equal(p.data, self.df)

    def test_establish_variables_from_bad(self):

        p = lm._LinearPlotter()
        with nt.assert_raises(ValueError):
            p.establish_variables(None, x="x", y=self.df.y)

    def test_dropna(self):

        p = lm._LinearPlotter()
        p.establish_variables(self.df, x="x", y_na="y_na")
        pdt.assert_series_equal(p.x, self.df.x)
        pdt.assert_series_equal(p.y_na, self.df.y_na)

        p.dropna("x", "y_na")
        mask = self.df.y_na.notnull()
        pdt.assert_series_equal(p.x, self.df.x[mask])
        pdt.assert_series_equal(p.y_na, self.df.y_na[mask])


class TestRegressionPlotter(PlotTestCase):

    rs = np.random.RandomState(49)

    grid = np.linspace(-3, 3, 30)
    n_boot = 100
    bins_numeric = 3
    bins_given = [-1, 0, 1]

    df = pd.DataFrame(dict(x=rs.normal(size=60),
                           d=rs.randint(-2, 3, 60),
                           y=rs.gamma(4, size=60),
                           s=np.tile(list(range(6)), 10)))
    df["z"] = df.y + rs.randn(60)
    df["y_na"] = df.y.copy()

    bw_err = rs.randn(6)[df.s.values] * 2
    df.y += bw_err

    p = 1 / (1 + np.exp(-(df.x * 2 + rs.randn(60))))
    df["c"] = [rs.binomial(1, p_i) for p_i in p]
    df.loc[[10, 20, 30], 'y_na'] = np.nan

    def test_variables_from_frame(self):

        p = lm._RegressionPlotter("x", "y", data=self.df, units="s")

        pdt.assert_series_equal(p.x, self.df.x)
        pdt.assert_series_equal(p.y, self.df.y)
        pdt.assert_series_equal(p.units, self.df.s)
        pdt.assert_frame_equal(p.data, self.df)

    def test_variables_from_series(self):

        p = lm._RegressionPlotter(self.df.x, self.df.y, units=self.df.s)

        npt.assert_array_equal(p.x, self.df.x)
        npt.assert_array_equal(p.y, self.df.y)
        npt.assert_array_equal(p.units, self.df.s)
        nt.assert_is(p.data, None)

    def test_variables_from_mix(self):

        p = lm._RegressionPlotter("x", self.df.y + 1, data=self.df)

        npt.assert_array_equal(p.x, self.df.x)
        npt.assert_array_equal(p.y, self.df.y + 1)
        pdt.assert_frame_equal(p.data, self.df)

    def test_dropna(self):

        p = lm._RegressionPlotter("x", "y_na", data=self.df)
        nt.assert_equal(len(p.x), pd.notnull(self.df.y_na).sum())

        p = lm._RegressionPlotter("x", "y_na", data=self.df, dropna=False)
        nt.assert_equal(len(p.x), len(self.df.y_na))

    def test_ci(self):

        p = lm._RegressionPlotter("x", "y", data=self.df, ci=95)
        nt.assert_equal(p.ci, 95)
        nt.assert_equal(p.x_ci, 95)

        p = lm._RegressionPlotter("x", "y", data=self.df, ci=95, x_ci=68)
        nt.assert_equal(p.ci, 95)
        nt.assert_equal(p.x_ci, 68)

        p = lm._RegressionPlotter("x", "y", data=self.df, ci=95, x_ci="sd")
        nt.assert_equal(p.ci, 95)
        nt.assert_equal(p.x_ci, "sd")

    @skipif(_no_statsmodels)
    def test_fast_regression(self):

        p = lm._RegressionPlotter("x", "y", data=self.df, n_boot=self.n_boot)

        # Fit with the "fast" function, which just does linear algebra
        yhat_fast, _ = p.fit_fast(self.grid)

        # Fit using the statsmodels function with an OLS model
        yhat_smod, _ = p.fit_statsmodels(self.grid, smlm.OLS)

        # Compare the vector of y_hat values
        npt.assert_array_almost_equal(yhat_fast, yhat_smod)

    @skipif(_no_statsmodels)
    def test_regress_poly(self):

        p = lm._RegressionPlotter("x", "y", data=self.df, n_boot=self.n_boot)

        # Fit an first-order polynomial
        yhat_poly, _ = p.fit_poly(self.grid, 1)

        # Fit using the statsmodels function with an OLS model
        yhat_smod, _ = p.fit_statsmodels(self.grid, smlm.OLS)

        # Compare the vector of y_hat values
        npt.assert_array_almost_equal(yhat_poly, yhat_smod)

    def test_regress_logx(self):

        x = np.arange(1, 10)
        y = np.arange(1, 10)
        grid = np.linspace(1, 10, 100)
        p = lm._RegressionPlotter(x, y, n_boot=self.n_boot)

        yhat_lin, _ = p.fit_fast(grid)
        yhat_log, _ = p.fit_logx(grid)

        nt.assert_greater(yhat_lin[0], yhat_log[0])
        nt.assert_greater(yhat_log[20], yhat_lin[20])
        nt.assert_greater(yhat_lin[90], yhat_log[90])

    @skipif(_no_statsmodels)
    def test_regress_n_boot(self):

        p = lm._RegressionPlotter("x", "y", data=self.df, n_boot=self.n_boot)

        # Fast (linear algebra) version
        _, boots_fast = p.fit_fast(self.grid)
        npt.assert_equal(boots_fast.shape, (self.n_boot, self.grid.size))

        # Slower (np.polyfit) version
        _, boots_poly = p.fit_poly(self.grid, 1)
        npt.assert_equal(boots_poly.shape, (self.n_boot, self.grid.size))

        # Slowest (statsmodels) version
        _, boots_smod = p.fit_statsmodels(self.grid, smlm.OLS)
        npt.assert_equal(boots_smod.shape, (self.n_boot, self.grid.size))

    @skipif(_no_statsmodels)
    def test_regress_without_bootstrap(self):

        p = lm._RegressionPlotter("x", "y", data=self.df,
                                  n_boot=self.n_boot, ci=None)

        # Fast (linear algebra) version
        _, boots_fast = p.fit_fast(self.grid)
        nt.assert_is(boots_fast, None)

        # Slower (np.polyfit) version
        _, boots_poly = p.fit_poly(self.grid, 1)
        nt.assert_is(boots_poly, None)

        # Slowest (statsmodels) version
        _, boots_smod = p.fit_statsmodels(self.grid, smlm.OLS)
        nt.assert_is(boots_smod, None)

    def test_numeric_bins(self):

        p = lm._RegressionPlotter(self.df.x, self.df.y)
        x_binned, bins = p.bin_predictor(self.bins_numeric)
        npt.assert_equal(len(bins), self.bins_numeric)
        npt.assert_array_equal(np.unique(x_binned), bins)

    def test_provided_bins(self):

        p = lm._RegressionPlotter(self.df.x, self.df.y)
        x_binned, bins = p.bin_predictor(self.bins_given)
        npt.assert_array_equal(np.unique(x_binned), self.bins_given)

    def test_bin_results(self):

        p = lm._RegressionPlotter(self.df.x, self.df.y)
        x_binned, bins = p.bin_predictor(self.bins_given)
        nt.assert_greater(self.df.x[x_binned == 0].min(),
                          self.df.x[x_binned == -1].max())
        nt.assert_greater(self.df.x[x_binned == 1].min(),
                          self.df.x[x_binned == 0].max())

    def test_scatter_data(self):

        p = lm._RegressionPlotter(self.df.x, self.df.y)
        x, y = p.scatter_data
        npt.assert_array_equal(x, self.df.x)
        npt.assert_array_equal(y, self.df.y)

        p = lm._RegressionPlotter(self.df.d, self.df.y)
        x, y = p.scatter_data
        npt.assert_array_equal(x, self.df.d)
        npt.assert_array_equal(y, self.df.y)

        p = lm._RegressionPlotter(self.df.d, self.df.y, x_jitter=.1)
        x, y = p.scatter_data
        nt.assert_true((x != self.df.d).any())
        npt.assert_array_less(np.abs(self.df.d - x), np.repeat(.1, len(x)))
        npt.assert_array_equal(y, self.df.y)

        p = lm._RegressionPlotter(self.df.d, self.df.y, y_jitter=.05)
        x, y = p.scatter_data
        npt.assert_array_equal(x, self.df.d)
        npt.assert_array_less(np.abs(self.df.y - y), np.repeat(.1, len(y)))

    def test_estimate_data(self):

        p = lm._RegressionPlotter(self.df.d, self.df.y, x_estimator=np.mean)

        x, y, ci = p.estimate_data

        npt.assert_array_equal(x, np.sort(np.unique(self.df.d)))
        npt.assert_array_almost_equal(y, self.df.groupby("d").y.mean())
        npt.assert_array_less(np.array(ci)[:, 0], y)
        npt.assert_array_less(y, np.array(ci)[:, 1])

    def test_estimate_cis(self):

        # set known good seed to avoid the test stochastically failing
        np.random.seed(123)

        p = lm._RegressionPlotter(self.df.d, self.df.y,
                                  x_estimator=np.mean, ci=95)
        _, _, ci_big = p.estimate_data

        p = lm._RegressionPlotter(self.df.d, self.df.y,
                                  x_estimator=np.mean, ci=50)
        _, _, ci_wee = p.estimate_data
        npt.assert_array_less(np.diff(ci_wee), np.diff(ci_big))

        p = lm._RegressionPlotter(self.df.d, self.df.y,
                                  x_estimator=np.mean, ci=None)
        _, _, ci_nil = p.estimate_data
        npt.assert_array_equal(ci_nil, [None] * len(ci_nil))

    def test_estimate_units(self):

        # Seed the RNG locally
        np.random.seed(345)

        p = lm._RegressionPlotter("x", "y", data=self.df,
                                  units="s", x_bins=3)
        _, _, ci_big = p.estimate_data
        ci_big = np.diff(ci_big, axis=1)

        p = lm._RegressionPlotter("x", "y", data=self.df, x_bins=3)
        _, _, ci_wee = p.estimate_data
        ci_wee = np.diff(ci_wee, axis=1)

        npt.assert_array_less(ci_wee, ci_big)

    def test_partial(self):

        x = self.rs.randn(100)
        y = x + self.rs.randn(100)
        z = x + self.rs.randn(100)

        p = lm._RegressionPlotter(y, z)
        _, r_orig = np.corrcoef(p.x, p.y)[0]

        p = lm._RegressionPlotter(y, z, y_partial=x)
        _, r_semipartial = np.corrcoef(p.x, p.y)[0]
        nt.assert_less(r_semipartial, r_orig)

        p = lm._RegressionPlotter(y, z, x_partial=x, y_partial=x)
        _, r_partial = np.corrcoef(p.x, p.y)[0]
        nt.assert_less(r_partial, r_orig)

    @skipif(_no_statsmodels)
    def test_logistic_regression(self):

        p = lm._RegressionPlotter("x", "c", data=self.df,
                                  logistic=True, n_boot=self.n_boot)
        _, yhat, _ = p.fit_regression(x_range=(-3, 3))
        npt.assert_array_less(yhat, 1)
        npt.assert_array_less(0, yhat)

    @skipif(_no_statsmodels)
    def test_logistic_perfect_separation(self):

        y = self.df.x > self.df.x.mean()
        p = lm._RegressionPlotter("x", y, data=self.df,
                                  logistic=True, n_boot=10)
        _, yhat, _ = p.fit_regression(x_range=(-3, 3))
        nt.assert_true(np.isnan(yhat).all())

    @skipif(_no_statsmodels)
    def test_robust_regression(self):

        p_ols = lm._RegressionPlotter("x", "y", data=self.df,
                                      n_boot=self.n_boot)
        _, ols_yhat, _ = p_ols.fit_regression(x_range=(-3, 3))

        p_robust = lm._RegressionPlotter("x", "y", data=self.df,
                                         robust=True, n_boot=self.n_boot)
        _, robust_yhat, _ = p_robust.fit_regression(x_range=(-3, 3))

        nt.assert_equal(len(ols_yhat), len(robust_yhat))

    @skipif(_no_statsmodels)
    def test_lowess_regression(self):

        p = lm._RegressionPlotter("x", "y", data=self.df, lowess=True)
        grid, yhat, err_bands = p.fit_regression(x_range=(-3, 3))

        nt.assert_equal(len(grid), len(yhat))
        nt.assert_is(err_bands, None)

    def test_regression_options(self):

        with nt.assert_raises(ValueError):
            lm._RegressionPlotter("x", "y", data=self.df,
                                  lowess=True, order=2)

        with nt.assert_raises(ValueError):
            lm._RegressionPlotter("x", "y", data=self.df,
                                  lowess=True, logistic=True)

    def test_regression_limits(self):

        f, ax = plt.subplots()
        ax.scatter(self.df.x, self.df.y)
        p = lm._RegressionPlotter("x", "y", data=self.df)
        grid, _, _ = p.fit_regression(ax)
        xlim = ax.get_xlim()
        nt.assert_equal(grid.min(), xlim[0])
        nt.assert_equal(grid.max(), xlim[1])

        p = lm._RegressionPlotter("x", "y", data=self.df, truncate=True)
        grid, _, _ = p.fit_regression()
        nt.assert_equal(grid.min(), self.df.x.min())
        nt.assert_equal(grid.max(), self.df.x.max())


class TestRegressionPlots(PlotTestCase):

    rs = np.random.RandomState(56)
    df = pd.DataFrame(dict(x=rs.randn(90),
                           y=rs.randn(90) + 5,
                           z=rs.randint(0, 1, 90),
                           g=np.repeat(list("abc"), 30),
                           h=np.tile(list("xy"), 45),
                           u=np.tile(np.arange(6), 15)))
    bw_err = rs.randn(6)[df.u.values]
    df.y += bw_err

    def test_regplot_basic(self):

        f, ax = plt.subplots()
        lm.regplot("x", "y", self.df)
        nt.assert_equal(len(ax.lines), 1)
        nt.assert_equal(len(ax.collections), 2)

        x, y = ax.collections[0].get_offsets().T
        npt.assert_array_equal(x, self.df.x)
        npt.assert_array_equal(y, self.df.y)

    def test_regplot_selective(self):

        f, ax = plt.subplots()
        ax = lm.regplot("x", "y", self.df, scatter=False, ax=ax)
        nt.assert_equal(len(ax.lines), 1)
        nt.assert_equal(len(ax.collections), 1)
        ax.clear()

        f, ax = plt.subplots()
        ax = lm.regplot("x", "y", self.df, fit_reg=False)
        nt.assert_equal(len(ax.lines), 0)
        nt.assert_equal(len(ax.collections), 1)
        ax.clear()

        f, ax = plt.subplots()
        ax = lm.regplot("x", "y", self.df, ci=None)
        nt.assert_equal(len(ax.lines), 1)
        nt.assert_equal(len(ax.collections), 1)
        ax.clear()

    def test_regplot_scatter_kws_alpha(self):

        f, ax = plt.subplots()
        color = np.array([[0.3, 0.8, 0.5, 0.5]])
        ax = lm.regplot("x", "y", self.df, scatter_kws={'color': color})
        nt.assert_is(ax.collections[0]._alpha, None)
        nt.assert_equal(ax.collections[0]._facecolors[0, 3], 0.5)

        f, ax = plt.subplots()
        color = np.array([[0.3, 0.8, 0.5]])
        ax = lm.regplot("x", "y", self.df, scatter_kws={'color': color})
        nt.assert_equal(ax.collections[0]._alpha, 0.8)

        f, ax = plt.subplots()
        color = np.array([[0.3, 0.8, 0.5]])
        ax = lm.regplot("x", "y", self.df, scatter_kws={'color': color,
                                                        'alpha': 0.4})
        nt.assert_equal(ax.collections[0]._alpha, 0.4)

        f, ax = plt.subplots()
        color = 'r'
        ax = lm.regplot("x", "y", self.df, scatter_kws={'color': color})
        nt.assert_equal(ax.collections[0]._alpha, 0.8)

    def test_regplot_binned(self):

        ax = lm.regplot("x", "y", self.df, x_bins=5)
        nt.assert_equal(len(ax.lines), 6)
        nt.assert_equal(len(ax.collections), 2)

    def test_lmplot_basic(self):

        g = lm.lmplot("x", "y", self.df)
        ax = g.axes[0, 0]
        nt.assert_equal(len(ax.lines), 1)
        nt.assert_equal(len(ax.collections), 2)

        x, y = ax.collections[0].get_offsets().T
        npt.assert_array_equal(x, self.df.x)
        npt.assert_array_equal(y, self.df.y)

    def test_lmplot_hue(self):

        g = lm.lmplot("x", "y", data=self.df, hue="h")
        ax = g.axes[0, 0]

        nt.assert_equal(len(ax.lines), 2)
        nt.assert_equal(len(ax.collections), 4)

    def test_lmplot_markers(self):

        g1 = lm.lmplot("x", "y", data=self.df, hue="h", markers="s")
        nt.assert_equal(g1.hue_kws, {"marker": ["s", "s"]})

        g2 = lm.lmplot("x", "y", data=self.df, hue="h", markers=["o", "s"])
        nt.assert_equal(g2.hue_kws, {"marker": ["o", "s"]})

        with nt.assert_raises(ValueError):
            lm.lmplot("x", "y", data=self.df, hue="h", markers=["o", "s", "d"])

    def test_lmplot_marker_linewidths(self):

        if mpl.__version__ == "1.4.2":
            raise SkipTest

        g = lm.lmplot("x", "y", data=self.df, hue="h",
                      fit_reg=False, markers=["o", "+"])
        c = g.axes[0, 0].collections
        nt.assert_equal(c[0].get_linewidths()[0], 0)
        rclw = mpl.rcParams["lines.linewidth"]
        nt.assert_equal(c[1].get_linewidths()[0], rclw)

    def test_lmplot_facets(self):

        g = lm.lmplot("x", "y", data=self.df, row="g", col="h")
        nt.assert_equal(g.axes.shape, (3, 2))

        g = lm.lmplot("x", "y", data=self.df, col="u", col_wrap=4)
        nt.assert_equal(g.axes.shape, (6,))

        g = lm.lmplot("x", "y", data=self.df, hue="h", col="u")
        nt.assert_equal(g.axes.shape, (1, 6))

    def test_lmplot_hue_col_nolegend(self):

        g = lm.lmplot("x", "y", data=self.df, col="h", hue="h")
        nt.assert_is(g._legend, None)

    def test_lmplot_scatter_kws(self):

        g = lm.lmplot("x", "y", hue="h", data=self.df, ci=None)
        red_scatter, blue_scatter = g.axes[0, 0].collections

        red, blue = color_palette(n_colors=2)
        npt.assert_array_equal(red, red_scatter.get_facecolors()[0, :3])
        npt.assert_array_equal(blue, blue_scatter.get_facecolors()[0, :3])

    def test_residplot(self):

        x, y = self.df.x, self.df.y
        ax = lm.residplot(x, y)

        resid = y - np.polyval(np.polyfit(x, y, 1), x)
        x_plot, y_plot = ax.collections[0].get_offsets().T

        npt.assert_array_equal(x, x_plot)
        npt.assert_array_almost_equal(resid, y_plot)

    @skipif(_no_statsmodels)
    def test_residplot_lowess(self):

        ax = lm.residplot("x", "y", self.df, lowess=True)
        nt.assert_equal(len(ax.lines), 2)

        x, y = ax.lines[1].get_xydata().T
        npt.assert_array_equal(x, np.sort(self.df.x))

    def test_three_point_colors(self):

        x, y = np.random.randn(2, 3)
        ax = lm.regplot(x, y, color=(1, 0, 0))
        color = ax.collections[0].get_facecolors()
        npt.assert_almost_equal(color[0, :3],
                                (1, 0, 0))
