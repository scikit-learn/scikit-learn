import numpy as np
from scipy import integrate

try:
    import statsmodels.distributions as smdist
except ImportError:
    smdist = None

import pytest
from numpy.testing import assert_array_equal, assert_array_almost_equal

from .._statistics import (
    KDE,
    Histogram,
    ECDF,
)


class DistributionFixtures:

    @pytest.fixture
    def x(self, rng):
        return rng.normal(0, 1, 100)

    @pytest.fixture
    def y(self, rng):
        return rng.normal(0, 5, 100)

    @pytest.fixture
    def weights(self, rng):
        return rng.uniform(0, 5, 100)


class TestKDE:

    def test_gridsize(self, rng):

        x = rng.normal(0, 3, 1000)

        n = 200
        kde = KDE(gridsize=n)
        density, support = kde(x)
        assert density.size == n
        assert support.size == n

    def test_cut(self, rng):

        x = rng.normal(0, 3, 1000)

        kde = KDE(cut=0)
        _, support = kde(x)
        assert support.min() == x.min()
        assert support.max() == x.max()

        cut = 2
        bw_scale = .5
        bw = x.std() * bw_scale
        kde = KDE(cut=cut, bw_method=bw_scale, gridsize=1000)
        _, support = kde(x)
        assert support.min() == pytest.approx(x.min() - bw * cut, abs=1e-2)
        assert support.max() == pytest.approx(x.max() + bw * cut, abs=1e-2)

    def test_clip(self, rng):

        x = rng.normal(0, 3, 100)
        clip = -1, 1
        kde = KDE(clip=clip)
        _, support = kde(x)

        assert support.min() >= clip[0]
        assert support.max() <= clip[1]

    def test_density_normalization(self, rng):

        x = rng.normal(0, 3, 1000)
        kde = KDE()
        density, support = kde(x)
        assert integrate.trapz(density, support) == pytest.approx(1, abs=1e-5)

    def test_cumulative(self, rng):

        x = rng.normal(0, 3, 1000)
        kde = KDE(cumulative=True)
        density, _ = kde(x)
        assert density[0] == pytest.approx(0, abs=1e-5)
        assert density[-1] == pytest.approx(1, abs=1e-5)

    def test_cached_support(self, rng):

        x = rng.normal(0, 3, 100)
        kde = KDE()
        kde.define_support(x)
        _, support = kde(x[(x > -1) & (x < 1)])
        assert_array_equal(support, kde.support)

    def test_bw_method(self, rng):

        x = rng.normal(0, 3, 100)
        kde1 = KDE(bw_method=.2)
        kde2 = KDE(bw_method=2)

        d1, _ = kde1(x)
        d2, _ = kde2(x)

        assert np.abs(np.diff(d1)).mean() > np.abs(np.diff(d2)).mean()

    def test_bw_adjust(self, rng):

        x = rng.normal(0, 3, 100)
        kde1 = KDE(bw_adjust=.2)
        kde2 = KDE(bw_adjust=2)

        d1, _ = kde1(x)
        d2, _ = kde2(x)

        assert np.abs(np.diff(d1)).mean() > np.abs(np.diff(d2)).mean()

    def test_bivariate_grid(self, rng):

        n = 100
        x, y = rng.normal(0, 3, (2, 50))
        kde = KDE(gridsize=n)
        density, (xx, yy) = kde(x, y)

        assert density.shape == (n, n)
        assert xx.size == n
        assert yy.size == n

    def test_bivariate_normalization(self, rng):

        x, y = rng.normal(0, 3, (2, 50))
        kde = KDE(gridsize=100)
        density, (xx, yy) = kde(x, y)

        dx = xx[1] - xx[0]
        dy = yy[1] - yy[0]

        total = density.sum() * (dx * dy)
        assert total == pytest.approx(1, abs=1e-2)

    def test_bivariate_cumulative(self, rng):

        x, y = rng.normal(0, 3, (2, 50))
        kde = KDE(gridsize=100, cumulative=True)
        density, _ = kde(x, y)

        assert density[0, 0] == pytest.approx(0, abs=1e-2)
        assert density[-1, -1] == pytest.approx(1, abs=1e-2)


class TestHistogram(DistributionFixtures):

    def test_string_bins(self, x):

        h = Histogram(bins="sqrt")
        bin_kws = h.define_bin_params(x)
        assert bin_kws["range"] == (x.min(), x.max())
        assert bin_kws["bins"] == int(np.sqrt(len(x)))

    def test_int_bins(self, x):

        n = 24
        h = Histogram(bins=n)
        bin_kws = h.define_bin_params(x)
        assert bin_kws["range"] == (x.min(), x.max())
        assert bin_kws["bins"] == n

    def test_array_bins(self, x):

        bins = [-3, -2, 1, 2, 3]
        h = Histogram(bins=bins)
        bin_kws = h.define_bin_params(x)
        assert_array_equal(bin_kws["bins"], bins)

    def test_bivariate_string_bins(self, x, y):

        s1, s2 = "sqrt", "fd"

        h = Histogram(bins=s1)
        e1, e2 = h.define_bin_params(x, y)["bins"]
        assert_array_equal(e1, np.histogram_bin_edges(x, s1))
        assert_array_equal(e2, np.histogram_bin_edges(y, s1))

        h = Histogram(bins=(s1, s2))
        e1, e2 = h.define_bin_params(x, y)["bins"]
        assert_array_equal(e1, np.histogram_bin_edges(x, s1))
        assert_array_equal(e2, np.histogram_bin_edges(y, s2))

    def test_bivariate_int_bins(self, x, y):

        b1, b2 = 5, 10

        h = Histogram(bins=b1)
        e1, e2 = h.define_bin_params(x, y)["bins"]
        assert len(e1) == b1 + 1
        assert len(e2) == b1 + 1

        h = Histogram(bins=(b1, b2))
        e1, e2 = h.define_bin_params(x, y)["bins"]
        assert len(e1) == b1 + 1
        assert len(e2) == b2 + 1

    def test_bivariate_array_bins(self, x, y):

        b1 = [-3, -2, 1, 2, 3]
        b2 = [-5, -2, 3, 6]

        h = Histogram(bins=b1)
        e1, e2 = h.define_bin_params(x, y)["bins"]
        assert_array_equal(e1, b1)
        assert_array_equal(e2, b1)

        h = Histogram(bins=(b1, b2))
        e1, e2 = h.define_bin_params(x, y)["bins"]
        assert_array_equal(e1, b1)
        assert_array_equal(e2, b2)

    def test_binwidth(self, x):

        binwidth = .5
        h = Histogram(binwidth=binwidth)
        bin_kws = h.define_bin_params(x)
        n_bins = bin_kws["bins"]
        left, right = bin_kws["range"]
        assert (right - left) / n_bins == pytest.approx(binwidth)

    def test_bivariate_binwidth(self, x, y):

        w1, w2 = .5, 1

        h = Histogram(binwidth=w1)
        e1, e2 = h.define_bin_params(x, y)["bins"]
        assert np.all(np.diff(e1) == w1)
        assert np.all(np.diff(e2) == w1)

        h = Histogram(binwidth=(w1, w2))
        e1, e2 = h.define_bin_params(x, y)["bins"]
        assert np.all(np.diff(e1) == w1)
        assert np.all(np.diff(e2) == w2)

    def test_binrange(self, x):

        binrange = (-4, 4)
        h = Histogram(binrange=binrange)
        bin_kws = h.define_bin_params(x)
        assert bin_kws["range"] == binrange

    def test_bivariate_binrange(self, x, y):

        r1, r2 = (-4, 4), (-10, 10)

        h = Histogram(binrange=r1)
        e1, e2 = h.define_bin_params(x, y)["bins"]
        assert e1.min() == r1[0]
        assert e1.max() == r1[1]
        assert e2.min() == r1[0]
        assert e2.max() == r1[1]

        h = Histogram(binrange=(r1, r2))
        e1, e2 = h.define_bin_params(x, y)["bins"]
        assert e1.min() == r1[0]
        assert e1.max() == r1[1]
        assert e2.min() == r2[0]
        assert e2.max() == r2[1]

    def test_discrete_bins(self, rng):

        x = rng.binomial(20, .5, 100)
        h = Histogram(discrete=True)
        bin_kws = h.define_bin_params(x)
        assert bin_kws["range"] == (x.min() - .5, x.max() + .5)
        assert bin_kws["bins"] == (x.max() - x.min() + 1)

    def test_histogram(self, x):

        h = Histogram()
        heights, edges = h(x)
        heights_mpl, edges_mpl = np.histogram(x, bins="auto")

        assert_array_equal(heights, heights_mpl)
        assert_array_equal(edges, edges_mpl)

    def test_count_stat(self, x):

        h = Histogram(stat="count")
        heights, _ = h(x)
        assert heights.sum() == len(x)

    def test_density_stat(self, x):

        h = Histogram(stat="density")
        heights, edges = h(x)
        assert (heights * np.diff(edges)).sum() == 1

    def test_probability_stat(self, x):

        h = Histogram(stat="probability")
        heights, _ = h(x)
        assert heights.sum() == 1

    def test_frequency_stat(self, x):

        h = Histogram(stat="frequency")
        heights, edges = h(x)
        assert (heights * np.diff(edges)).sum() == len(x)

    def test_cumulative_count(self, x):

        h = Histogram(stat="count", cumulative=True)
        heights, _ = h(x)
        assert heights[-1] == len(x)

    def test_cumulative_density(self, x):

        h = Histogram(stat="density", cumulative=True)
        heights, _ = h(x)
        assert heights[-1] == 1

    def test_cumulative_probability(self, x):

        h = Histogram(stat="probability", cumulative=True)
        heights, _ = h(x)
        assert heights[-1] == 1

    def test_cumulative_frequency(self, x):

        h = Histogram(stat="frequency", cumulative=True)
        heights, _ = h(x)
        assert heights[-1] == len(x)

    def test_bivariate_histogram(self, x, y):

        h = Histogram()
        heights, edges = h(x, y)
        bins_mpl = (
            np.histogram_bin_edges(x, "auto"),
            np.histogram_bin_edges(y, "auto"),
        )
        heights_mpl, *edges_mpl = np.histogram2d(x, y, bins_mpl)
        assert_array_equal(heights, heights_mpl)
        assert_array_equal(edges[0], edges_mpl[0])
        assert_array_equal(edges[1], edges_mpl[1])

    def test_bivariate_count_stat(self, x, y):

        h = Histogram(stat="count")
        heights, _ = h(x, y)
        assert heights.sum() == len(x)

    def test_bivariate_density_stat(self, x, y):

        h = Histogram(stat="density")
        heights, (edges_x, edges_y) = h(x, y)
        areas = np.outer(np.diff(edges_x), np.diff(edges_y))
        assert (heights * areas).sum() == pytest.approx(1)

    def test_bivariate_probability_stat(self, x, y):

        h = Histogram(stat="probability")
        heights, _ = h(x, y)
        assert heights.sum() == 1

    def test_bivariate_frequency_stat(self, x, y):

        h = Histogram(stat="frequency")
        heights, (x_edges, y_edges) = h(x, y)
        area = np.outer(np.diff(x_edges), np.diff(y_edges))
        assert (heights * area).sum() == len(x)

    def test_bivariate_cumulative_count(self, x, y):

        h = Histogram(stat="count", cumulative=True)
        heights, _ = h(x, y)
        assert heights[-1, -1] == len(x)

    def test_bivariate_cumulative_density(self, x, y):

        h = Histogram(stat="density", cumulative=True)
        heights, _ = h(x, y)
        assert heights[-1, -1] == pytest.approx(1)

    def test_bivariate_cumulative_frequency(self, x, y):

        h = Histogram(stat="frequency", cumulative=True)
        heights, _ = h(x, y)
        assert heights[-1, -1] == len(x)

    def test_bivariate_cumulative_probability(self, x, y):

        h = Histogram(stat="probability", cumulative=True)
        heights, _ = h(x, y)
        assert heights[-1, -1] == pytest.approx(1)

    def test_bad_stat(self):

        with pytest.raises(ValueError):
            Histogram(stat="invalid")


class TestECDF(DistributionFixtures):

    def test_univariate_proportion(self, x):

        ecdf = ECDF()
        stat, vals = ecdf(x)
        assert_array_equal(vals[1:], np.sort(x))
        assert_array_almost_equal(stat[1:], np.linspace(0, 1, len(x) + 1)[1:])
        assert stat[0] == 0

    def test_univariate_count(self, x):

        ecdf = ECDF(stat="count")
        stat, vals = ecdf(x)

        assert_array_equal(vals[1:], np.sort(x))
        assert_array_almost_equal(stat[1:], np.arange(len(x)) + 1)
        assert stat[0] == 0

    def test_univariate_proportion_weights(self, x, weights):

        ecdf = ECDF()
        stat, vals = ecdf(x, weights=weights)
        assert_array_equal(vals[1:], np.sort(x))
        expected_stats = weights[x.argsort()].cumsum() / weights.sum()
        assert_array_almost_equal(stat[1:], expected_stats)
        assert stat[0] == 0

    def test_univariate_count_weights(self, x, weights):

        ecdf = ECDF(stat="count")
        stat, vals = ecdf(x, weights=weights)
        assert_array_equal(vals[1:], np.sort(x))
        assert_array_almost_equal(stat[1:], weights[x.argsort()].cumsum())
        assert stat[0] == 0

    @pytest.mark.skipif(smdist is None, reason="Requires statsmodels")
    def test_against_statsmodels(self, x):

        sm_ecdf = smdist.empirical_distribution.ECDF(x)

        ecdf = ECDF()
        stat, vals = ecdf(x)
        assert_array_equal(vals, sm_ecdf.x)
        assert_array_almost_equal(stat, sm_ecdf.y)

        ecdf = ECDF(complementary=True)
        stat, vals = ecdf(x)
        assert_array_equal(vals, sm_ecdf.x)
        assert_array_almost_equal(stat, sm_ecdf.y[::-1])

    def test_invalid_stat(self, x):

        with pytest.raises(ValueError, match="`stat` must be one of"):
            ECDF(stat="density")

    def test_bivariate_error(self, x, y):

        with pytest.raises(NotImplementedError, match="Bivariate ECDF"):
            ecdf = ECDF()
            ecdf(x, y)
