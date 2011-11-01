"""
==================================
l1-penalized covariance estimation
==================================

Using the l1-penalized covariance estimator to learn a covariance and
precision from a small number of samples.

To estimate a probabilistic model (e.g. a Gaussian model), estimating the
precision matrix, that is the inverse covariance matrix, is as important
as estimating the covariance matrix. Indeed a Gaussian model is
parametrized by the precision matrix.

To be in favorable recovery conditions, we sample the data from a model
with a sparse inverse covariance matrix. In addition, we ensure that the
data is not too much correlated (limiting the largest coefficient of the
precision matrix) and that there a no small coefficients in the
precision matrix that cannot be recovered.

Here, the number of samples is slightly larger than the number of
dimensions, thus the empirical covariance is still invertible. However,
as the observations are strongly correlated, the empirical covariance
matrix is ill-conditioned and as a result its inverse --the empirical
precision matrix-- is very far from the ground truth.

If we use l2 shrinkage, as with the Ledoit-Wolf estimator, as the number
of samples is small, we need to shrink a lot. As a result, the
Ledoit-Wolf precision is fairly close to the ground truth precision, that
is not far from being diagonal, but the off-diagonal structure is lost.

The l1-penalized estimator can recover part of this off-diagonal
structure. It learns a sparse precision. It is not able to
recover the exact sparsity pattern: it detects too many non-zero
coefficients. However, the highest non-zero coefficients of the l1
estimated correspond to the non-zero coefficients in the ground truth.
Finally, the coefficients of the l1 estimate are biased toward zero:
because of the penalty, they are all smaller than the corresponding
ground truth value.

Note that, the color range of the precision matrices is tweeked to
improve readibility of the figure. The full range of values of the
empirical precision is not displayed.
"""
# author: Gael Varoquaux <gael.varoquaux@inria.fr>
# License: BSD Style
# Copyright: INRIA

import numpy as np
from scipy import linalg
from sklearn.datasets.samples_generator import make_sparse_spd_matrix
from sklearn.covariance import g_lasso, ledoit_wolf
import pylab as pl

################################################################################
# Generate the data
N_SAMPLES = 15
DIM = 10

prng = np.random.RandomState(1)
prec = make_sparse_spd_matrix(DIM, norm_diag=True, smallest_coef=.4,
                              largest_coef=.7,
                              random_state=prng)
cov = linalg.inv(prec)
X = prng.multivariate_normal(np.zeros(DIM), cov, size=N_SAMPLES)
X -= X.mean(axis=0)

################################################################################
# Estimate the covariance
emp_cov = np.dot(X.T, X)/N_SAMPLES

alpha = .3
cov_, prec_ = g_lasso(X, alpha=alpha)
lw_cov_, _ = ledoit_wolf(X)
lw_prec_ = linalg.inv(lw_cov_)

################################################################################
# Plot the results
pl.figure(figsize=(10, 6))
pl.subplots_adjust(left=0.02, right=0.98)

# plot the covariances
covs = [('Empirical', emp_cov), ('Ledoit-Wolf', lw_cov_),
        ('L1', cov_), ('True', cov)]
vmax = cov_.max()
for i, (name, this_cov) in enumerate(covs):
    pl.subplot(2, 4, i+1)
    pl.imshow(this_cov, interpolation='nearest', vmin=-vmax, vmax=vmax,
                cmap=pl.cm.RdBu_r)
    pl.xticks(())
    pl.yticks(())
    pl.title('%s covariance' % name)


# plot the precisions
precs = [('Empirical', linalg.inv(emp_cov)), ('Ledoit-Wolf', lw_prec_),
         ('L1', prec_), ('True', prec)]
vmax = .9*prec_.max()
for i, (name, this_prec) in enumerate(precs):
    ax = pl.subplot(2, 4, i+5)
    pl.imshow(np.ma.masked_equal(this_prec, 0),
                interpolation='nearest', vmin=-vmax, vmax=vmax,
                cmap=pl.cm.RdBu_r)
    pl.xticks(())
    pl.yticks(())
    pl.title('%s precision' % name)
    ax.set_axis_bgcolor('.7')


pl.show()

