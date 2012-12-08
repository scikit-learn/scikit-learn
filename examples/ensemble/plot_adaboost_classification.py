"""
==============================
Adaboosted Class Probabilities
==============================

This example fits an Adaboosted decision tree on a classification dataset and
plots the decision surfaces and class probabilities.
"""
print __doc__

from itertools import izip

import pylab as pl
import numpy as np

from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.datasets import make_classification

X, y = make_classification(n_samples=1000,
                           n_features=2,
                           n_classes=2,
                           n_informative=2,
                           n_redundant=0)

bdt = AdaBoostClassifier(DecisionTreeClassifier(min_samples_leaf=100),
                         n_estimators=50,
                         learn_rate=.5)

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
norm = sum(bdt.boost_weights_)
for weight, tree in zip(bdt.boost_weights_, bdt.estimators_):
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

# plot the class probabilities
pl.subplot(122)
pl.hist(bdt.predict_proba(X[y==0])[:,-1], bins=20, range=(0, 1),
        facecolor=plot_colors[0],
        label='Class A')
pl.hist(bdt.predict_proba(X[y==1])[:,-1], bins=20, range=(0, 1),
        facecolor=plot_colors[1],
        label='Class B')
pl.legend()
pl.ylabel('Samples')
pl.xlabel('Class Probability')

pl.show()
