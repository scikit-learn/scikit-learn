"""
=========================================
Compare kernel CCA with different kernels
=========================================

Compare the correlations found by the first componennt extracted from
CCA, kernel CCA with linear kernel, polynomial kernel and rbf kernel.

Given 2 multivariate covarying two-dimensional datasets, X, and Y,
CCA finds the directions to project X and Y separately to maxmize the
correlation between the projections. Kernel CCA finds the directions
to project the kernel matrices Kx (computed from X) and Ky (computed
from Y) to maximize the correlation between the projections. In this
example, the correlation found for the first component are shown in
the scatter plots.

Linear kernel CCA gives almost the same correlation as CCA; polynomial
kernel implicitly defines a richer feature space, which leads to higher
correlation; rbf kernel implicitly defines a feature space with infinite
dimensions, thus gives the highest correlation comparing with the others.
"""
print(__doc__)

# Authors: Huibin Shen
# License: BSD 3 clause

import numpy as np
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sklearn.cross_decomposition import CCA
from sklearn.cross_decomposition import KernelCCA

###############################################################################
# Dataset based latent variables model

n = 500

# 2 latents vars:
l1 = np.random.normal(size=n)
l2 = np.random.normal(size=n)

latents = np.array([l1, l1, l2, l2]).T
X = latents + np.random.normal(size=4 * n).reshape((n, 4))
Y = latents + np.random.normal(size=4 * n).reshape((n, 4))

###############################################################################
# Compare the projection on first component of CCA, kernel CCA
# with linear kernel, polynomial kernel and rbf kernel

cca = CCA(n_components=1)
cca.fit(X, Y)
r_cca = np.corrcoef(cca.x_scores_.T, cca.y_scores_.T)[0, 1]

# linear kernel CCA
kcca1 = KernelCCA(kernel="linear", n_components=1, kapa=0.1,
                  eta=0.1, pgso=True, center=True)
kcca1.fit(X, Y)
kx_linear_scores = np.dot(kcca1.KXc_, kcca1.alphas_)
ky_linear_scores = np.dot(kcca1.KYc_, kcca1.betas_)

# polynomial kernel CCA
kcca2 = KernelCCA(kernel="poly", n_components=1, kapa=0.1,
                  eta=0.1, pgso=True, center=True, coef0=0.1)
kcca2.fit(X, Y)
kx_poly_scores = np.dot(kcca2.KXc_, kcca2.alphas_)
ky_poly_scores = np.dot(kcca2.KYc_, kcca2.betas_)

# rbf kernel CCA
kcca3 = KernelCCA(kernel="rbf", n_components=1, kapa=0.1,
                  eta=0.1, pgso=True, center=True)
kcca3.fit(X, Y)
kx_rbf_scores = np.dot(kcca3.KXc_, kcca3.alphas_)
ky_rbf_scores = np.dot(kcca3.KYc_, kcca3.betas_)


# Scatter plot of first component
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1) CCA, scatter plot of first component for X and Y

plt.figure(figsize=(14, 10))
plt.subplot(221)
plt.plot(cca_x_scores, cca_y_scores, 'o')
plt.xlabel("x scores")
plt.ylabel("y scores")
plt.title('Comp. 1: CCA (corr = %.2f)' % r_cca)
plt.xticks(())
plt.yticks(())

# 2) Linear kernel CCA, scatter plot of first component for KX and KY
plt.subplot(222)
plt.plot(kx_linear_scores, ky_linear_scores, 'o')
plt.xlabel("x scores")
plt.ylabel("y scores")
plt.title('Comp. 1: linear kernel CCA (corr = %.2f)' % kcca1.lambdas_[0])
plt.xticks(())
plt.yticks(())

# 3) Ploynomial kernel CCA, scatter plot of first component for KX and KY
plt.subplot(223)
plt.plot(kx_poly_scores, ky_poly_scores, 'o')
plt.xlabel("x scores")
plt.ylabel("y scores")
plt.title('Comp. 1: polynomial kernel CCA (corr = %.2f)' % kcca2.lambdas_[0])
plt.xticks(())
plt.yticks(())

# 4) Rbf kernel CCA, scatter plot of first component for KX and KY
plt.subplot(224)
plt.plot(kx_rbf_scores, ky_rbf_scores, 'o')
plt.xlabel("x scores")
plt.ylabel("y scores")
plt.title('Comp. 1: rbf kernel CCA (corr = %.2f)' % kcca3.lambdas_[0])
plt.xticks(())
plt.yticks(())
plt.show()
#plt.savefig('kcca.pdf')
