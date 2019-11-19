"""
=================================
Conditional Gaussian Mixture Model
=================================

The ConditionalGaussianMixture class takes an (already trained) Gaussian 
mixture model as an input, and can then be used to evaluate conditional
distributions of the mixture. In the following we generate data from a
2D Gaussian mixture before using the samples to try and estimate the 
resulting 1D conditional distribution. We then fit our conditional 
Gaussian mixture and show that it matches closely with the approximate 
solution. 
"""
print(__doc__)

# Author: P L Green <p.l.green@liverpool.ac.uk> 

import numpy as np
from matplotlib import pyplot as plt
from numpy.random import multivariate_normal as mvn
from sklearn import mixture

# Generate samples from a 2D Gaussian Mixture Model
n_components=2
X1 = mvn(mean=[3, 3], cov=np.eye(2), size=50000)
X2 = mvn(mean=[-4, 0], cov=np.array([[1, 0.8],[0.8, 1]]), size=50000)
X = np.vstack([X1, X2])

# Fit a standard Gaussian mixture model
gmm = mixture.GaussianMixture(n_components=n_components)
gmm.fit(X)

# We are going to hold x2 equal to 1 and evaluate p(x1|x2) for different
# values of x1. 
x2 = 1.0

# Here we use the existing samples to generate a histogram of the conditional
# distribution p(x1|x2)
indices = np.where((X[:, 1] > x2 - 0.1) &
                   (X[:, 1] < x2 + 0.1))
                   
fig, ax = plt.subplots(nrows=2, ncols=1)
ax[0].plot(X[:, 0], X[:, 1], 'k x')
ax[0].plot(X[indices, 0], X[indices, 1], 'r o')
ax[0].set_xlabel('x1')
ax[0].set_ylabel('x2')
ax[0].set_xlim([-8, 8])
ax[0].set_title('Gaussian mixture')
ax[1].hist(X[indices, 0].T, density=True, bins=50)
ax[1].set_xlim([-8, 8])

# Create our GMM_Conditional object. Note that the indices where
# i_cond is True represent the indices of x that will be held at a 
# conditional value. 
gmm_cond = mixture.ConditionalGaussianMixture(gmm=gmm,
                                              i_cond=np.array([False, True]))

# Plot conditional distribution over range of x1 values
x1_range = np.linspace(-7, 7, 1000)
cond_pdf = np.zeros(1000)
for i in range(1000):
    cond_pdf[i] = gmm_cond.pdf_xa_cond_xb(x1_range[i], x2)
ax[1].plot(x1_range, cond_pdf, 'k')
ax[1].set_xlim([-8, 8])
ax[1].set_title('Conditional of Gaussian mixture')
ax[1].set_xlabel('x1')
ax[1].set_ylabel('p(x1 | x2)')

# Tidy-up and display plots
plt.tight_layout()
plt.show()
