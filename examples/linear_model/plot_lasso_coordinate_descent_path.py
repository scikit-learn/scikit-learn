"""
=====================
Lasso and Elastic Net
=====================

Lasso and elastic net (L1 and L2 penalisation) implemented using a
coordinate descent.

The coefficients can be forced to be positive.
"""
print(__doc__)

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD 3 clause

import numpy as np
import pylab as pl

from sklearn.linear_model import lasso_path, enet_path, lasso_path_cd, \
    enet_path_cd
from sklearn import datasets

diabetes = datasets.load_diabetes()
X = diabetes.data
y = diabetes.target

X /= X.std(0)  # Standardize data (easier to set the l1_ratio parameter)

###############################################################################
# Compute paths

eps = 5e-3  # the smaller it is the longer is the path

print("Computing regularization path using the lasso...")
#models = lasso_path(X, y, eps=eps)
#alphas_lasso = np.array([model.alpha for model in models])
#coefs_lasso = np.array([model.coef_ for model in models])

alphas_lasso, coefs_lasso = lasso_path_cd(X, y, eps)

print("Computing regularization path using the positive lasso...")
#models = lasso_path(X, y, eps=eps, positive=True)
#alphas_positive_lasso = np.array([model.alpha for model in models])
#coefs_positive_lasso = np.array([model.coef_ for model in models])

alphas_positve_lasso, coefs_positive_lasso = lasso_path_cd(X, y,
                                                                 eps,
                                                                 positive=True)
print("Computing regularization path using the elastic net...")
#models = enet_path(X, y, eps=eps, l1_ratio=0.8)
#alphas_enet = np.array([model.alpha for model in models])
#coefs_enet = np.array([model.coef_ for model in models])

alphas_enet, coefs_enet = enet_path_cd(X, y, eps=eps, l1_ratio=0.8)

print("Computing regularization path using the positve elastic net...")
#models = enet_path(X, y, eps=eps, l1_ratio=0.8, positive=True)
#alphas_positive_enet = np.array([model.alpha for model in models])
#coefs_positive_enet = np.array([model.coef_ for model in models])
alphas_positive_enet, coefs_positive_enet = enet_path_cd(X, y,
                                                         eps=eps,l1_ratio=0.8,
                                                         positive=True)

###############################################################################
# Display results

pl.figure(1)
ax = pl.gca()
ax.set_color_cycle(2 * ['b', 'r', 'g', 'c', 'k'])
l1 = pl.plot(coefs_lasso.T)
l2 = pl.plot(coefs_enet.T, linestyle='--')

pl.xlabel('-Log(lambda)')
pl.ylabel('coefficients')
pl.title('Lasso and Elastic-Net Paths')
pl.legend((l1[-1], l2[-1]), ('Lasso', 'Elastic-Net'), loc='lower left')
pl.axis('tight')


pl.figure(2)
ax = pl.gca()
ax.set_color_cycle(2 * ['b', 'r', 'g', 'c', 'k'])
l1 = pl.plot(coefs_lasso.T)
l2 = pl.plot(coefs_positive_lasso.T, linestyle='--')

pl.xlabel('-Log(lambda)')
pl.ylabel('coefficients')
pl.title('Lasso and positive Lasso')
pl.legend((l1[-1], l2[-1]), ('Lasso', 'positive Lasso'), loc='lower left')
pl.axis('tight')


pl.figure(3)
ax = pl.gca()
ax.set_color_cycle(2 * ['b', 'r', 'g', 'c', 'k'])
l1 = pl.plot(coefs_enet.T)
l2 = pl.plot(coefs_positive_enet.T, linestyle='--')

pl.xlabel('-Log(lambda)')
pl.ylabel('coefficients')
pl.title('Elastic-Net and positive Elastic-Net')
pl.legend((l1[-1], l2[-1]), ('Elastic-Net', 'positive Elastic-Net'),
          loc='lower left')
pl.axis('tight')
pl.show()
