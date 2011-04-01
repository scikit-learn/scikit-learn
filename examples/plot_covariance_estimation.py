"""
===========================================
Ledoit-Wolf vs Covariance simple estimation
===========================================

Covariance estimation can be regularized using a shrinkage parameter.
Ledoit-Wolf estimates automatically this parameter. In this example,
we compute the likelihood of unseen data for different values of
the shrinkage parameter. The Ledoit-Wolf estimate reaches an
almost optimal value.
"""
print __doc__

import numpy as np
import pylab as pl
from scipy import linalg

###############################################################################
# Generate sample data
n_features, n_samples = 30, 20
X_train = np.random.normal(size=(n_samples, n_features))
X_test = np.random.normal(size=(n_samples, n_features))

# Color samples
coloring_matrix = np.random.normal(size=(n_features, n_features))
X_train = np.dot(X_train, coloring_matrix)
X_test = np.dot(X_test, coloring_matrix)

###############################################################################
# Compute Ledoit-Wolf and Covariances on a grid of shrinkages

from scikits.learn.covariance import LedoitWolf, OAS, ShrunkCovariance
from scikits.learn.covariance import Covariance, log_likelihood

lw = LedoitWolf()
loglik_lw = lw.fit(X_train).score(X_test)

oa = OAS()
loglik_oa = oa.fit(X_train).score(X_test)

shrinkages = np.logspace(-2, 0, 30)
negative_logliks = [-ShrunkCovariance(shrinkage=s).fit(X_train).score(X_test) \
                                                        for s in shrinkages]

emp_cov = Covariance(store_precision=False).fit(X_test).covariance_
real_cov = np.dot(coloring_matrix.T, coloring_matrix)
loglik_real = -log_likelihood(emp_cov, linalg.inv(real_cov))

###############################################################################
# Plot results
pl.loglog(shrinkages, negative_logliks)
pl.xlabel('Shrinkage')
pl.ylabel('Negative log-likelihood')
pl.vlines(lw.shrinkage_, pl.ylim()[0], -loglik_lw, color='g',
          linewidth=3, label='Ledoit-Wolf estimate')
pl.vlines(oa.shrinkage_, pl.ylim()[0], -loglik_oa, color='orange',
          linewidth=3, label='OAS estimate')
pl.hlines(loglik_real, pl.xlim()[0], pl.xlim()[1], color='red',
          label="Real covariance")
pl.ylim(pl.ylim()[0] - (pl.ylim()[1]-pl.ylim()[0]), pl.ylim()[1])
pl.legend()
pl.show()
