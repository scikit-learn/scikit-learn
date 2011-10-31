"""
==================================
l1-penalized covariance estimation
==================================

"""
# author: Gael Varoquaux <gael.varoquaux@inria.fr>
# License: BSD Style
# Copyright: INRIA

import numpy as np
from scipy import linalg
from sklearn.datasets.samples_generator import make_sparse_spd_matrix
from sklearn.covariance import g_lasso
import pylab as pl

################################################################################
# Generate the data
N_SAMPLES = 40
DIM = 20

prng = np.random.RandomState(0)
prec = make_sparse_spd_matrix(DIM, alpha=.95, random_state=prng)
cov = linalg.inv(prec)
X = prng.multivariate_normal(np.zeros(DIM), cov, size=N_SAMPLES)
X -= X.mean(axis=0)

################################################################################
# Estimate the covariance
emp_cov = np.dot(X.T, X)/N_SAMPLES

alpha = .2
cov_, prec_ = g_lasso(X, alpha=alpha)

################################################################################
# Plot the results
pl.figure()

covs = {'True': cov, 'L1': cov_, 'Empirical': emp_cov}
vmin = np.min(covs.values())
vmax = np.max(covs.values())
vmax = max(-vmin, vmax)

for i, (name, this_cov) in enumerate(sorted(covs.iteritems())):
    pl.subplot(2, 3, i+1)
    pl.imshow(this_cov, interpolation='nearest', vmin=-vmax, vmax=vmax,
                cmap=pl.cm.RdBu_r)
    pl.xticks(())
    pl.yticks(())
    pl.title('%s covariance' % name)


precs = {'True': prec, 'L1': prec_, 'Empirical': linalg.inv(emp_cov)}
vmax = 3

for i, (name, this_prec) in enumerate(sorted(precs.iteritems())):
    pl.subplot(2, 3, i+4)
    pl.imshow(np.ma.masked_equal(this_prec, 0),
                interpolation='nearest', vmin=-vmax, vmax=vmax,
                cmap=pl.cm.Spectral_r)
    pl.xticks(())
    pl.yticks(())
    pl.title('%s precision' % name)


pl.show()

