"""
=================================
Conditional Gaussian Mixture Model
=================================

stuff to add

"""

import numpy as np
from matplotlib import pyplot as plt
from numpy.random import multivariate_normal as mvn
from sklearn import mixture

# Generate samples from a Gaussian Mixture Model
n_components=2
X1 = mvn(mean=[3, 3], cov=np.eye(2), size=50000)
X2 = mvn(mean=[-4, 0], cov=np.array([[1, 0.8],[0.8, 1]]), size=50000)
X = np.vstack([X1, X2])

# Fit a Gaussian mixture model
gmm = mixture.GaussianMixture(n_components=n_components)
gmm.fit(X)

# Our conditional value of x2
x2 = 1.0

# Histogram approximation of conditional distribution
indices = np.where((X[:, 1] > x2 - 0.1) &
                   (X[:, 1] < x2 + 0.1))

# Create histogram
fig, ax = plt.subplots(nrows=2, ncols=1)
ax[0].plot(X[:, 0], X[:, 1], 'k x')
ax[0].plot(X[indices, 0], X[indices, 1], 'r o')
ax[0].set_xlabel('x1')
ax[0].set_ylabel('x2')
ax[0].set_xlim([-8, 8])
ax[1].hist(X[indices, 0].T, density=True, bins=50)
ax[1].set_xlim([-8, 8])

# Create our GMM_Conditional object
gmm_cond = mixture.ConditionalGaussianMixture(gmm=gmm,
                                              i_cond=np.array([False, True]))

# Plot conditional distribution over range of x1 values
x1_range = np.linspace(-7, 7, 1000)
pdf = np.zeros(1000)
for i in range(1000):
    pdf[i] = gmm_cond.pdf_xa_cond_xb(x1_range[i], x2)
ax[1].plot(x1_range, pdf, 'k')
ax[1].set_xlim([-8, 8])
ax[1].set_xlabel('x1')
ax[1].set_ylabel('p(x1 | x2)')

# Tidy-up and display plots
plt.tight_layout()
plt.show()
