"""
======================================
LASSO model selection with BIC and AIC
======================================

This examples illustrates how the BIC and AIC criterion
can be used to estimate the regularization i.e. the number
of predictive features when using the LASSO estimator.

The estimation of the number of degrees of freedom is given by:

"On the degrees of freedom of the lasso"
Hui Zou, Trevor Hastie, and Robert Tibshirani
Ann. Statist. Volume 35, Number 5 (2007), 2173-2192.
"""
print __doc__

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

from math import log
import numpy as np

from scikits.learn.linear_model import lars_path
from scikits.learn.datasets import load_diabetes

diabetes = load_diabetes()
X, y = diabetes.data, diabetes.target

# add garbage features
rng = np.random.RandomState(42)
X = np.c_[X, rng.randn(X.shape[0], 10)]
n_samples = X.shape[0]

# Standardize the data to avoid intercept problems
y -= np.mean(y)
X -= np.mean(X, axis=0)
X /= np.std(X, axis=0)

print "Computing regularization path using the LARS ..."
alphas, features, coefs = lars_path(X, y, method='lasso', verbose=True)

###############################################################################
# BIC and AIC
K_aic = 2  # AIC
K_bic = log(n_samples)  # BIC

R = y[:, np.newaxis] - np.dot(X, coefs)  # residuals
mse = np.sum(R ** 2, axis=0)  # MSE ie. mean square error

df = np.zeros(coefs.shape[1], dtype=np.int)  # Degrees of freedom
for k, coef in enumerate(coefs.T):
    mask = coef != 0
    if not np.any(mask):
        continue
    Xc = X[:, mask]
    # get the number of degrees of freedom equal to:
    # Trace(Xc * inv(Xc.T, Xc) * Xc.T) ie the number of non-zero coefs
    df[k] = np.sum(mask)

aic_criterion = np.log(mse) + K_aic / float(n_samples) * df
bic_criterion = np.log(mse) + K_bic / float(n_samples) * df
n_aic = np.argmin(aic_criterion)
n_bic = np.argmin(bic_criterion)

print "Optimal number of features:"
print "    AIC : %d" % n_aic
print "    BIC : %d" % n_bic
aic_idx = np.where(coefs[:, n_aic] != 0)[0]
bic_idx = np.where(coefs[:, n_bic] != 0)[0]
print "Features indices for AIC : %s" % aic_idx
print "Features indices for BIC : %s" % bic_idx

###############################################################################
# plot AIC and BIC criterion along the path
import pylab as pl
pl.clf()
pl.plot(aic_criterion, label='AIC crit.')
pl.plot(bic_criterion, label='BIC crit.')
pl.vlines(n_aic, pl.ylim()[0], aic_criterion[n_aic], color='b',
          linewidth=3, label='AIC estimate')
pl.vlines(n_bic, pl.ylim()[0], bic_criterion[n_bic], color='g',
          linewidth=3, label='BIC estimate')
pl.xlabel('Nb of features')
pl.ylabel('criterion')
pl.legend()
pl.show()
