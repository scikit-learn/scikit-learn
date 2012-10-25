"""
===========================================
Ledoit-Wolf vs Covariance simple estimation
===========================================

The usual covariance maximum likelihood estimator can be regularized
in order to reduce its variance. This, in turn, introduces some bias.
This example illustrates a regularization method called convex shrinkage that
consists in replacing the maximum likelihood estimator by a convex
combination of itself and an identity matrix (*). Various convex shrinkage-
based covariance estimators are shown in the example.

Ledoit and Wolf proposed a close formula to compute
the asymptotical optimal shrinkage parameter (minimizing a MSE
criterion), yielding the Ledoit-Wolf covariance estimate.

Chen et al. proposed an improvement of the Ledoit-Wolf shrinkage
parameter, the OAS coefficient, whose convergence is significantly
better under the assumption that the data are gaussian.

In this example, we compute the likelihood of unseen data for
different values of the shrinkage parameter, highlighting the LW and
OAS estimates. The Ledoit-Wolf estimate stays close to the likelihood
criterion optimal value, which is an artifact of the method since it
is asymptotic and we are working with a small number of observations.
The OAS estimate deviates from the likelihood criterion optimal value
but better approximate the MSE optimal value, especially for a small
number a observations.
We also show the best shrunk estimate that we obtained by cross-validating
the likelihood on three folds according to a grid of potential shrinkage
parameters.

(*) Although any other structured target could be used instead, the
identity matrix is the weakest assumption that can be made to regularize
a covariance estimator.

"""
print __doc__

import numpy as np
import pylab as pl
from scipy import linalg
from matplotlib.patches import Polygon

###############################################################################
# Generate sample data
n_features, n_samples = 30, 20
base_X_train = np.random.normal(size=(n_samples, n_features))
base_X_test = np.random.normal(size=(n_samples, n_features))

# Color samples
coloring_matrix = np.random.normal(size=(n_features, n_features))
X_train = np.dot(base_X_train, coloring_matrix)
X_test = np.dot(base_X_test, coloring_matrix)

###############################################################################
# Compute Ledoit-Wolf and Covariances on a grid of shrinkages

from sklearn.covariance import LedoitWolf, OAS, ShrunkCovariance, \
    log_likelihood, empirical_covariance
from sklearn.grid_search import GridSearchCV

# Ledoit-Wolf optimal shrinkage coefficient estimate
lw = LedoitWolf()
loglik_lw = lw.fit(X_train).score(X_test)

# OAS coefficient estimate
oa = OAS()
loglik_oa = oa.fit(X_train).score(X_test)

# spanning a range of possible shrinkage coefficient values
shrinkages = np.logspace(-2, 0, 30)
negative_logliks = [-ShrunkCovariance(shrinkage=s).fit(X_train).score(X_test)
                     for s in shrinkages]

# GridSearch for an optimal shrinkage coefficient
tuned_parameters = [{'shrinkage': shrinkages}]
cv = GridSearchCV(ShrunkCovariance(), tuned_parameters)
cv.fit(X_train)

# getting the likelihood under the real model
real_cov = np.dot(coloring_matrix.T, coloring_matrix)
emp_cov = empirical_covariance(X_train)
loglik_real = -log_likelihood(emp_cov, linalg.inv(real_cov))

###############################################################################
# Plot results
fig = pl.figure()
pl.title("Regularized covariance: likelihood and shrinkage coefficient")
pl.xlabel('Shrinkage')
pl.ylabel('Negative log-likelihood')
# range shrinkage curve
pl.loglog(shrinkages, negative_logliks, label="Negative log-likelihood")

# real likelihood reference
# BUG: hlines(..., linestyle='--') breaks on some older versions of matplotlib
#pl.hlines(loglik_real, pl.xlim()[0], pl.xlim()[1], color='red',
#          label="real covariance likelihood", linestyle='--')
pl.plot(pl.xlim(), 2 * [loglik_real], '--r',
        label="real covariance likelihood")

# adjust view
lik_max = np.amax(negative_logliks)
lik_min = np.amin(negative_logliks)
ymin = lik_min - 6. * np.log((pl.ylim()[1] - pl.ylim()[0]))
ymax = lik_max + 10. * np.log(lik_max - lik_min)
xmin = shrinkages[0]
xmax = shrinkages[-1]
# LW likelihood
pl.vlines(lw.shrinkage_, ymin, -loglik_lw, color='magenta',
          linewidth=3, label='Ledoit-Wolf estimate')
# OAS likelihood
pl.vlines(oa.shrinkage_, ymin, -loglik_oa, color='purple',
          linewidth=3, label='OAS estimate')
# best CV estimator likelihood
pl.vlines(cv.best_estimator_.shrinkage, ymin,
          -cv.best_estimator_.score(X_test), color='cyan',
          linewidth=3, label='Cross-validation best estimate')

pl.ylim(ymin, ymax)
pl.xlim(xmin, xmax)
pl.legend()

pl.show()
