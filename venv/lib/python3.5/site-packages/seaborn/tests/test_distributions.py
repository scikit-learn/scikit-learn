import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from distutils.version import LooseVersion
import nose.tools as nt
import numpy.testing as npt
from numpy.testing.decorators import skipif

from . import PlotTestCase
from .. import distributions as dist

try:
    import statsmodels.nonparametric.api
    assert statsmodels.nonparametric.api
    _no_statsmodels = False
except ImportError:
    _no_statsmodels = True


_old_matplotlib = LooseVersion(mpl.__version__) < "1.5"


class TestKDE(PlotTestCase):

    rs = np.random.RandomState(0)
    x = rs.randn(50)
    y = rs.randn(50)
    kernel = "gau"
    bw = "scott"
    gridsize = 128
    clip = (-np.inf, np.inf)
    cut = 3

    def test_scipy_univariate_kde(self):
        """Test the univariate KDE estimation with scipy."""
        grid, y = dist._scipy_univariate_kde(self.x, self.bw, self.gridsize,
                                             self.cut, self.clip)
        nt.assert_equal(len(grid), self.gridsize)
        nt.assert_equal(len(y), self.gridsize)
        for bw in ["silverman", .2]:
            dist._scipy_univariate_kde(self.x, bw, self.gridsize,
                                       self.cut, self.clip)

    @skipif(_no_statsmodels)
    def test_statsmodels_univariate_kde(self):
        """Test the univariate KDE estimation with statsmodels."""
        grid, y = dist._statsmodels_univariate_kde(self.x, self.kernel,
                                                   self.bw, self.gridsize,
                                                   self.cut, self.clip)
        nt.assert_equal(len(grid), self.gridsize)
        nt.assert_equal(len(y), self.gridsize)
        for bw in ["silverman", .2]:
            dist._statsmodels_univariate_kde(self.x, self.kernel, bw,
                                             self.gridsize, self.cut,
                                             self.clip)

    def test_scipy_bivariate_kde(self):
        """Test the bivariate KDE estimation with scipy."""
        clip = [self.clip, self.clip]
        x, y, z = dist._scipy_bivariate_kde(self.x, self.y, self.bw,
                                            self.gridsize, self.cut, clip)
        nt.assert_equal(x.shape, (self.gridsize, self.gridsize))
        nt.assert_equal(y.shape, (self.gridsize, self.gridsize))
        nt.assert_equal(len(z), self.gridsize)

        # Test a specific bandwidth
        clip = [self.clip, self.clip]
        x, y, z = dist._scipy_bivariate_kde(self.x, self.y, 1,
                                            self.gridsize, self.cut, clip)

        # Test that we get an error with an invalid bandwidth
        with nt.assert_raises(ValueError):
            dist._scipy_bivariate_kde(self.x, self.y, (1, 2),
                                      self.gridsize, self.cut, clip)

    @skipif(_no_statsmodels)
    def test_statsmodels_bivariate_kde(self):
        """Test the bivariate KDE estimation with statsmodels."""
        clip = [self.clip, self.clip]
        x, y, z = dist._statsmodels_bivariate_kde(self.x, self.y, self.bw,
                                                  self.gridsize,
                                                  self.cut, clip)
        nt.assert_equal(x.shape, (self.gridsize, self.gridsize))
        nt.assert_equal(y.shape, (self.gridsize, self.gridsize))
        nt.assert_equal(len(z), self.gridsize)

    @skipif(_no_statsmodels)
    def test_statsmodels_kde_cumulative(self):
        """Test computation of cumulative KDE."""
        grid, y = dist._statsmodels_univariate_kde(self.x, self.kernel,
                                                   self.bw, self.gridsize,
                                                   self.cut, self.clip,
                                                   cumulative=True)
        nt.assert_equal(len(grid), self.gridsize)
        nt.assert_equal(len(y), self.gridsize)
        # make sure y is monotonically increasing
        npt.assert_((np.diff(y) > 0).all())

    def test_kde_cummulative_2d(self):
        """Check error if args indicate bivariate KDE and cumulative."""
        with npt.assert_raises(TypeError):
            dist.kdeplot(self.x, data2=self.y, cumulative=True)

    def test_bivariate_kde_series(self):
        df = pd.DataFrame({'x': self.x, 'y': self.y})

        ax_series = dist.kdeplot(df.x, df.y)
        ax_values = dist.kdeplot(df.x.values, df.y.values)

        nt.assert_equal(len(ax_series.collections),
                        len(ax_values.collections))
        nt.assert_equal(ax_series.collections[0].get_paths(),
                        ax_values.collections[0].get_paths())

    @skipif(_old_matplotlib)
    def test_bivariate_kde_colorbar(self):

        f, ax = plt.subplots()
        dist.kdeplot(self.x, self.y,
                     cbar=True, cbar_kws=dict(label="density"),
                     ax=ax)
        nt.assert_equal(len(f.axes), 2)
        nt.assert_equal(f.axes[1].get_ylabel(), "density")
