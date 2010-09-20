"""
=====================
Lasso and Elastic Net
=====================

Lasso and elastic net (L1 and L2 penalisation) implemented using a
coordinate descent.
"""
print __doc__

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

from itertools import cycle
import numpy as np
import pylab as pl

from scikits.learn.glm import lasso_path, enet_path
from scikits.learn import datasets

diabetes = datasets.load_diabetes()
X = diabetes.data
y = diabetes.target

X /= X.std(0) # Standardize data (easier to set the rho parameter)

################################################################################
# Compute paths

eps = 5e-3 # the smaller it is the longer is the path

print "Computing regularization path using the lasso..."
models = lasso_path(X, y, eps=eps)
alphas_lasso = np.array([model.alpha for model in models])
coefs_lasso = np.array([model.coef_ for model in models])

print "Computing regularization path using the elastic net..."
models = enet_path(X, y, eps=eps, rho=0.8)
alphas_enet = np.array([model.alpha for model in models])
coefs_enet = np.array([model.coef_ for model in models])

################################################################################
# Display results

color_iter = cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
for color, coef_lasso, coef_enet in zip(color_iter,
                            coefs_lasso.T, coefs_enet.T):
    pl.plot(-np.log10(alphas_lasso), coef_lasso, color)
    pl.plot(-np.log10(alphas_enet), coef_enet, color + '--')

pl.xlabel('-Log(lambda)')
pl.ylabel('weights')
pl.title('Lasso and Elastic-Net Paths')
pl.legend(['Lasso','Elastic-Net'])
pl.axis('tight')
pl.show()

