"""
================================================================
Robust covariance estimation and Mahalanobis distances relevance
================================================================

For Gaussian ditributed data, the distance of an observation $x_i$ to
the mode of the distribution can be computed using its Mahalanobis
distance: $d_{(\mu,\Sigma)}(x_i)^2 = (x_i - \mu)'\Sigma^{-1}(x_i -
\mu)$ where $\mu$ and $\Sigma$ are the location and the covariance of
the underlying gaussian distribution.

In practice, $\mu$ and $\Sigma$ are replaced by some estimates.  The
usual covariance maximum likelihood estimate is very sensitive to the
presence of outliers in the data set and therefor, the corresponding
Mahalanobis distances are. One would better have to use a robust
estimator of covariance to garanty that the estimation is resistant to
"errorneous" observations in the data set and that the associated
Mahalanobis distances accurately reflect the true organisation of the
observations.

The Minimum Covariance Determinant estimator is a robust,
high-breakdown point (i.e. it can be used to estimate the covariance
matrix of highly contaminated datasets, up to
$\frac{n_samples-n_features-1}{2}$ outliers) estimator of
covariance. The idea is to find $\frac{n_samples+n_features+1}{2}$
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

[1] P. J. Rousseeuw. Least median of squares regression. J. Am
    Stat Ass, 79:871, 1984.

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
# variables and parameters for cosmetic
offset_left = fig.subplotpars.left
offset_bottom = fig.subplotpars.bottom
width = fig.subplotpars.right - offset_left
subfig1 = pl.subplot(3, 1, 1)
subfig2 = pl.subplot(3, 1, 2)
subfig3 = pl.subplot(3, 1, 3)

# Show data set
subfig1.scatter(X[:, 0], X[:, 1], color='black', label='inliers')
subfig1.scatter(X[:, 0][-n_outliers:], X[:, 1][-n_outliers:],
                color='red', label='outliers')
subfig1.set_xlim(subfig1.get_xlim()[0], 11.)
subfig1.set_title("Mahalanobis distances of a contaminated data set:")
subfig1.legend(loc="upper right")

# Empirical covariance -based Mahalanobis distances
subfig2.scatter(np.arange(n_samples), emp_cov.mahalanobis(X),
                color='black', label='inliers')
subfig2.scatter(np.arange(n_samples)[-n_outliers:],
                emp_cov.mahalanobis(X)[-n_outliers:],
                color='red', label='outliers')
subfig2.set_ylabel("Mahal. dist.")
subfig2.set_title("1. from empirical estimates")
subfig2.axes.set_position(pos=[offset_left, 0.39, width, .2])

# MCD-based Mahalanobis distances
subfig3.scatter(np.arange(n_samples), robust_cov.mahalanobis(X),
                color='black', label='inliers')
subfig3.scatter(np.arange(n_samples)[-n_outliers:],
                robust_cov.mahalanobis(X)[-n_outliers:],
                color='red', label='outliers')
subfig3.set_ylabel("Mahal. dist.")
subfig3.set_title("2. from robust estimates (Minimum Covariance Determinant)")
subfig3.axes.set_position(pos=[offset_left, offset_bottom, width, .2])

pl.show()
