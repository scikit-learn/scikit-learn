"""
Lasso and elastic net (L1 and L2 penalisation) implemented using a
coordinate descent.
"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id$

from itertools import cycle
import numpy as np
import pylab as pl

from scikits.learn.glm.cd import Lasso, ElasticNet, lasso_path, enet_path

n_samples, n_features, maxit = 5, 10, 30

np.random.seed(0)
y = np.random.randn(n_samples)
X = np.random.randn(n_samples, n_features)

################################################################################
# Fit models 
################################################################################

# Lasso
lasso = Lasso(alpha=1)
lasso.fit(X, y, maxit=maxit)

print "Duality gap Lasso (should be small): %f" % lasso.compute_gap(X, y)


# Elastic-Net 
enet = ElasticNet(alpha=1, beta=1)
enet.fit(X, y, maxit=maxit)

print "Duality gap (should be small): %f" % enet.compute_gap(X, y)

# Display results
pl.figure(-1, figsize=(8, 4))
pl.clf()
pl.subplots_adjust(wspace=.4, right=.95)
pl.subplot(1, 2, 1)
pl.plot(lasso.objective, label='Lasso')
pl.plot(enet.objective,  label='Elastic Net')
pl.xlabel('Iteration')
pl.ylabel('Cost function')
pl.legend()
pl.title('Convergence')

################################################################################
# Demo path functions
################################################################################

alphas_lasso, weights_lasso = lasso_path(X, y, factor=0.97, n_alphas = 100)
alphas_enet, weights_enet = enet_path(X, y, factor=0.97, n_alphas = 100, 
                                                beta=0.1)

# Display results
pl.subplot(1, 2, 2)
color_iter = cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
for color, weight_lasso, weight_enet in zip(color_iter, 
                            weights_lasso.T, weights_enet.T):
    pl.plot(-np.log(alphas_lasso), weight_lasso, color)
    pl.plot(-np.log(alphas_enet), weight_enet, color+'x')

pl.xlabel('-log(lambda)')
pl.ylabel('weights')
pl.title('Lasso and Elastic-Net Paths')
pl.legend(['Lasso','Elastic-Net'])
pl.show()

