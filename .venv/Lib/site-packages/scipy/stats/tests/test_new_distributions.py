# file for distribution-specific tests with new infrastructure (UnivariateDistribution)
import numpy as np
from numpy.testing import assert_allclose
from scipy import stats

class TestBinomial:
    def test_gh23708_binomial_logcdf_method_complement(self):
        # gh-23708 found that `logcdf` method='complement' was inaccurate in the tails
        x = np.asarray([0., 18.])
        X = stats.Binomial(n=np.asarray([18.]), p=np.asarray(0.71022842))
        assert_allclose(X.logcdf(x, method='complement'), X.logcdf(x), rtol=1e-15)
        assert_allclose(X.logccdf(x, method='complement'), X.logccdf(x), rtol=1e-15)

        # going even deeper into the tails
        X = stats.Binomial(n=100, p=0.5)
        assert_allclose(X.logcdf(0, method='complement'), X.logpmf(0), rtol=1e-15)
        assert_allclose(X.logccdf(99, method='complement'), X.logpmf(100), rtol=1e-15)
