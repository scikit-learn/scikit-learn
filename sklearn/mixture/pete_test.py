import numpy as np
from numpy.random import multivariate_normal as mvn
from sklearn.mixture import GaussianMixture
from sklearn.mixture import ConditionalGaussianMixture

# Generate samples from a Gaussian Mixture Model
n_components=2
X1 = mvn(mean=[3, 3], cov=np.eye(2), size=1000)
X2 = mvn(mean=[-4, 0], cov=np.array([[1, 0.8],[0.8, 1]]), size=1000)
X = np.vstack([X1, X2])

# Fit a Gaussian mixture model
gmm = GaussianMixture(n_components=n_components)
gmm.fit(X)

cond_gmm = ConditionalGaussianMixture(gmm)
