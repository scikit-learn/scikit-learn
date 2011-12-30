"""
============================================================
Plot the decision surfaces on the Gaussian quantiles dataset
============================================================

This plot shows the decision surfaces learned by an AdaBoosted decision tree
classifier on the Gaussian quantiles dataset.
"""
print __doc__

import numpy as np
import pylab as pl

from sklearn.datasets import make_gaussian_quantiles
from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier

# Parameters
n_classes = 3
plot_colors = "bry"
plot_step = 0.02

model = AdaBoostClassifier(DecisionTreeClassifier(min_samples_leaf=10),
                           n_estimators=200)

# Load data
X, y = make_gaussian_quantiles(n_samples=500, n_features=2,
                               n_classes=n_classes)

# Train
model.fit(X, y)

# Plot the decision boundary
pl.subplot(111)

x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, plot_step),
                     np.arange(y_min, y_max, plot_step))

norm = sum(model.boost_weights_)
for weight, tree in zip(model.boost_weights_, model.estimators_):
    Z = tree.predict(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
    cs = pl.contourf(xx, yy, Z, alpha=weight / norm, cmap=pl.cm.Paired)

pl.axis("tight")

# Plot the training points
for i, c in zip(xrange(n_classes), plot_colors):
    idx = np.where(y == i)
    pl.scatter(X[idx, 0], X[idx, 1], c=c, cmap=pl.cm.Paired)

pl.axis("tight")

pl.suptitle("Decision surfaces of a boosted decision tree.")
pl.show()
