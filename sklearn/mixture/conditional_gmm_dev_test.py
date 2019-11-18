import numpy as np
from matplotlib import pyplot as plt
from numpy.random import multivariate_normal as mvn
from sklearn.mixture import GaussianMixture
from sklearn.mixture import ConditionalGaussianMixture

# Author: Peter L Green <p.l.green@liverpool.ac.uk>

"""
This script provides a simple visual check for the class
ConditionalGaussianMixture that we're developing for 
_gaussian_mixture.py

"""

# Generate samples from a Gaussian Mixture Model
n_components=2
X1 = mvn(mean=[3, 3], cov=np.eye(2), size=10000)
X2 = mvn(mean=[-5, 1], cov=np.array([[1, 0.8],[0.8, 1]]), size=10000)
X = np.vstack([X1, X2])

# Fit a Gaussian mixture model
gmm = GaussianMixture(n_components=n_components)
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
ax[0].set_xlim([-9, 9])
ax[1].hist(X[indices, 0].T, density=True, bins=50)
ax[1].set_xlim([-9, 9])

# Create our GMM_Conditional object
i_cond = np.array([False, True])
cond_gmm = ConditionalGaussianMixture(gmm, i_cond)

# Plot conditional distribution over range of x1 values
x1_range = np.linspace(-9, 9, 1000)
pdf = np.zeros(1000)
for i in range(1000):
    pdf[i] = cond_gmm.pdf_xa_cond_xb(x1_range[i], x2)
ax[1].plot(x1_range, pdf, 'k')
ax[1].set_xlim([-9, 9])
ax[1].set_xlabel('x1')
ax[1].set_ylabel('p(x1 | x2)')

# Tidy-up and display plots
plt.tight_layout()
plt.show()
