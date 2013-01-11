"""
==============================
Adaboosted Class Probabilities
==============================

This example fits an Adaboosted decision tree on a classification dataset and
plots the decision surfaces and class probabilities.
"""
print __doc__

import pylab as pl
import numpy as np

from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.datasets import make_classification

X, y = make_classification(n_samples=1000,
                           n_features=2,
                           n_classes=2,
                           n_informative=2,
                           n_redundant=0,
                           random_state=0)

bdt = AdaBoostClassifier()

bdt.fit(X, y)

plot_colors = "br"
plot_step = 0.02

pl.figure(figsize=(15, 5))

# plot the decision boundaries
pl.subplot(121)
x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, plot_step),
                     np.arange(y_min, y_max, plot_step))
norm = sum(bdt.weights_)
for weight, tree in zip(bdt.weights_, bdt.estimators_):
    Z = tree.predict(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)
    cs = pl.contourf(xx, yy, Z, alpha=weight / norm, cmap=pl.cm.Paired)
pl.axis("tight")

# Plot the training points
for i, c in zip(xrange(2), plot_colors):
    idx = np.where(y == i)
    pl.scatter(X[idx, 0], X[idx, 1], c=c, cmap=pl.cm.Paired)
pl.axis("tight")
pl.xlabel("Decision Surfaces")

a_prob = bdt.predict_proba(X[y == 0])[:, -1]
b_prob = bdt.predict_proba(X[y == 1])[:, -1]
hist_range = (min(a_prob.min(), b_prob.min()),
              max(a_prob.max(), b_prob.max()))

# plot the class probabilities
pl.subplot(122)
pl.hist(a_prob,
        bins=20,
        range=hist_range,
        facecolor=plot_colors[0],
        label='Class A')
pl.hist(b_prob,
        bins=20,
        range=hist_range,
        facecolor=plot_colors[1],
        label='Class B')
pl.legend()
pl.ylabel('Samples')
pl.xlabel('Class Probability')

pl.show()
