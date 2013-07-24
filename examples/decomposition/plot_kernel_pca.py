"""
==========
Kernel PCA
==========

This example shows that Kernel PCA is able to find a projection of the data
that makes data linearly separable.
"""
print(__doc__)

# Authors: Mathieu Blondel
#          Andreas Mueller
# License: BSD 3 clause

import numpy as np
import pylab as pl

from sklearn.decomposition import PCA, KernelPCA
from sklearn.datasets import make_circles

np.random.seed(0)

X, y = make_circles(n_samples=400, factor=.3, noise=.05)

kpca = KernelPCA(kernel="rbf", fit_inverse_transform=True, gamma=10)
X_kpca = kpca.fit_transform(X)
X_back = kpca.inverse_transform(X_kpca)
pca = PCA()
X_pca = pca.fit_transform(X)

# Plot results

pl.figure()
pl.subplot(2, 2, 1, aspect='equal')
pl.title("Original space")
reds = y == 0
blues = y == 1

pl.plot(X[reds, 0], X[reds, 1], "ro")
pl.plot(X[blues, 0], X[blues, 1], "bo")
pl.xlabel("$x_1$")
pl.ylabel("$x_2$")

X1, X2 = np.meshgrid(np.linspace(-1.5, 1.5, 50), np.linspace(-1.5, 1.5, 50))
X_grid = np.array([np.ravel(X1), np.ravel(X2)]).T
# projection on the first principal component (in the phi space)
Z_grid = kpca.transform(X_grid)[:, 0].reshape(X1.shape)
pl.contour(X1, X2, Z_grid, colors='grey', linewidths=1, origin='lower')

pl.subplot(2, 2, 2, aspect='equal')
pl.plot(X_pca[reds, 0], X_pca[reds, 1], "ro")
pl.plot(X_pca[blues, 0], X_pca[blues, 1], "bo")
pl.title("Projection by PCA")
pl.xlabel("1st principal component")
pl.ylabel("2nd component")

pl.subplot(2, 2, 3, aspect='equal')
pl.plot(X_kpca[reds, 0], X_kpca[reds, 1], "ro")
pl.plot(X_kpca[blues, 0], X_kpca[blues, 1], "bo")
pl.title("Projection by KPCA")
pl.xlabel("1st principal component in space induced by $\phi$")
pl.ylabel("2nd component")

pl.subplot(2, 2, 4, aspect='equal')
pl.plot(X_back[reds, 0], X_back[reds, 1], "ro")
pl.plot(X_back[blues, 0], X_back[blues, 1], "bo")
pl.title("Original space after inverse transform")
pl.xlabel("$x_1$")
pl.ylabel("$x_2$")

pl.subplots_adjust(0.02, 0.10, 0.98, 0.94, 0.04, 0.35)

pl.show()
