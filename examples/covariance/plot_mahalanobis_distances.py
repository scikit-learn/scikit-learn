"""
================================================================
Robust covariance estimation and Mahalanobis distances relevance
================================================================

For Gaussian ditributed data, the distance of an observation
:math:`x_i` to the mode of the distribution can be computed using its
Mahalanobis distance: :math:`d_{(\mu,\Sigma)}(x_i)^2 = (x_i -
\mu)'\Sigma^{-1}(x_i - \mu)` where :math:`\mu` and :math:`\Sigma` are
the location and the covariance of the underlying gaussian
distribution.

In practice, :math:`\mu` and :math:`\Sigma` are replaced by some
estimates.  The usual covariance maximum likelihood estimate is very
sensitive to the presence of outliers in the data set and therefor,
the corresponding Mahalanobis distances are. One would better have to
use a robust estimator of covariance to garanty that the estimation is
resistant to "errorneous" observations in the data set and that the
associated Mahalanobis distances accurately reflect the true
organisation of the observations.

The Minimum Covariance Determinant estimator is a robust,
high-breakdown point (i.e. it can be used to estimate the covariance
matrix of highly contaminated datasets, up to
:math:`\frac{n_samples-n_features-1}{2}` outliers) estimator of
covariance. The idea is to find :math:`\frac{n_samples+n_features+1}{2}`
observations whose empirical covariance has the smallest determinant,
yielding a "pure" subset of observations from which to compute
standards estimates of location and covariance.

The Minimum Covariance Determinant estimator (MCD) has been introduced
by P.J.Rousseuw in [1].

This example illustrates how the Mahalanobis distances are affected by
outlying data: observations drawn from a contaminating distribution
are not distinguishable from the observations comming from the real,
Gaussian distribution that one may want to work with. Using MCD-based
Mahalanobis distances, the two populations become
distinguishable. Associated applications are outliers detection,
observations ranking, clustering, ...
For vizualisation purpose, the cubique root of the Mahalanobis distances
are represented in the boxplot, as Wilson and Hilferty suggest [2]

[1] P. J. Rousseeuw. Least median of squares regression. J. Am
    Stat Ass, 79:871, 1984.
[2] Wilson, E. B., & Hilferty, M. M. (1931). The distribution of chi-square.
    Proceedings of the National Academy of Sciences of the United States
    of America, 17, 684-688.

"""
print __doc__

import numpy as np
import pylab as pl

from sklearn.covariance import EmpiricalCovariance, MinCovDet

n_samples = 125
n_outliers = 25
n_features = 2

# generate data
gen_cov = np.eye(n_features)
gen_cov[0, 0] = 2.
X = np.dot(np.random.randn(n_samples, n_features), gen_cov)
# add some outliers
outliers_cov = np.eye(n_features)
outliers_cov[np.arange(1, n_features), np.arange(1, n_features)] = 7.
X[-n_outliers:] = np.dot(np.random.randn(n_outliers, n_features), outliers_cov)

# fit a Minimum Covariance Determinant (MCD) robust estimator to data
robust_cov = MinCovDet().fit(X)

# compare estimators learnt from the full data set with true parameters
emp_cov = EmpiricalCovariance().fit(X)


# Display results
fig = pl.figure()

# Show data set
subfig1 = pl.subplot(3, 1, 1)
subfig1.scatter(X[:, 0], X[:, 1], color='black', label='inliers')
subfig1.scatter(X[:, 0][-n_outliers:], X[:, 1][-n_outliers:],
                color='red', label='outliers')
subfig1.set_xlim(subfig1.get_xlim()[0], 11.)
subfig1.set_title("Mahalanobis distances of a contaminated data set:")
subfig1.legend(loc="upper right")

emp_mahal = emp_cov.mahalanobis(X - np.mean(X, 0)) ** (0.33)
subfig2 = pl.subplot(2, 2, 3)
subfig2.boxplot([emp_mahal[:-n_outliers], emp_mahal[-n_outliers:]], widths=.25)
subfig2.plot(1.26 * np.ones(n_samples - n_outliers),
             emp_mahal[:-n_outliers], '+k', markeredgewidth=1)
subfig2.plot(2.26 * np.ones(n_outliers),
             emp_mahal[-n_outliers:], '+k', markeredgewidth=1)
subfig2.axes.set_xticklabels(('inliers', 'outliers'), size=11)
subfig2.set_ylabel(r"$\sqrt[3]{\rm{(Mahal. dist.)}}$")
subfig2.set_title("1. from non-robust estimates\n(Maximum Likelihood)")

robust_mahal = robust_cov.mahalanobis(X - robust_cov.location_) ** (0.33)
subfig3 = pl.subplot(2, 2, 4)
subfig3.boxplot([robust_mahal[:-n_outliers], robust_mahal[-n_outliers:]],
                widths=.25)
subfig3.plot(1.26 * np.ones(n_samples - n_outliers),
             robust_mahal[:-n_outliers], '+k', markeredgewidth=1)
subfig3.plot(2.26 * np.ones(n_outliers),
             robust_mahal[-n_outliers:], '+k', markeredgewidth=1)
subfig3.axes.set_xticklabels(('inliers', 'outliers'), size=11)
subfig3.set_ylabel(r"$\sqrt[3]{\rm{(Mahal. dist.)}}$")
subfig3.set_title("2. from robust estimates\n(Minimum Covariance Determinant)")

pl.show()
