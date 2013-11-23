"""
Multiclass Generative Sampling
==============================
This example shows the use of the Generative Bayesian classifier for sampling
from a multi-class distribution.

The first figure shows a simple 2D distribution, overlaying the input points
and new points generated from the class-wise model.

The second figure extends this to a higher dimension.  A generative Bayes
classifier based on kernel density estimation is fit to the handwritten digits
data, and a new sample is drawn from each of the class-wise generative
models.
"""
import matplotlib.pyplot as plt
import numpy as np
from sklearn.naive_bayes import GenerativeBayes
from sklearn.decomposition import PCA
from sklearn.grid_search import GridSearchCV
from sklearn.neighbors import KernelDensity
from sklearn.datasets import make_blobs, load_digits

#----------------------------------------------------------------------
# First figure: two-dimensional blobs

# Make 4 blobs with different numbers of points
np.random.seed(0)
X1, y1 = make_blobs(50, 2, centers=2)
X2, y2 = make_blobs(200, 2, centers=2)

X = np.vstack([X1, X2])
y = np.concatenate([y1, y2 + 2])

# Fit a generative Bayesian model to the data
clf = GenerativeBayes('norm_approx')
clf.fit(X, y)

# Sample new data from the generative Bayesian model
X_new, y_new = clf.sample(200)

# Plot the input data and the sampled data
fig, ax = plt.subplots()
ax.scatter(X[:, 0], X[:, 1], c=y, alpha=0.2)
ax.scatter(X_new[:, 0], X_new[:, 1], c=y_new)

# Create the legend by plotting some empty data
ax.scatter([], [], c='w', alpha=0.2, label="Training (input) data")
ax.scatter([], [], c='w', label="Samples from Model")
ax.legend()

ax.set_xlim(-4, 10)
ax.set_ylim(-8, 8)


#----------------------------------------------------------------------
# Second figure: sampling from digits digits

# load the digits data
digits = load_digits()
data = digits.data
labels = digits.target

# project the 64-dimensional data to a lower dimension
pca = PCA(n_components=15, whiten=False)
data = pca.fit_transform(digits.data)

# use grid search cross-validation to optimize the bandwidth
params = {'bandwidth': np.logspace(-1, 1, 20)}
grid = GridSearchCV(KernelDensity(), params)
grid.fit(data)

print "best bandwidth: {0}".format(grid.best_estimator_.bandwidth)

# train the model with this bandwidth
clf = GenerativeBayes('kde',
                      model_kwds={'bandwidth':grid.best_estimator_.bandwidth})
clf.fit(data, labels)

new_data, new_labels = clf.sample(44, random_state=0)
new_data = pca.inverse_transform(new_data)

# turn data into a 4x11 grid
new_data = new_data.reshape((4, 11, -1))
real_data = digits.data[:44].reshape((4, 11, -1))

# plot real digits and resampled digits
fig, ax = plt.subplots(9, 11, subplot_kw=dict(xticks=[], yticks=[]))
for j in range(11):
    ax[4, j].set_visible(False)
    for i in range(4):
        im = ax[i, j].imshow(real_data[i, j].reshape((8, 8)),
                             cmap=plt.cm.binary, interpolation='nearest')
        im.set_clim(0, 16)
        im = ax[i + 5, j].imshow(new_data[i, j].reshape((8, 8)),
                                 cmap=plt.cm.binary, interpolation='nearest')
        im.set_clim(0, 16)

ax[0, 5].set_title('Selection from the input data')
ax[5, 5].set_title('"New" digits drawn from the class-wise kernel density model')

plt.show()
