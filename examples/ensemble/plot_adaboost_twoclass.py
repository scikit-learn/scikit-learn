"""
==================
Two-class AdaBoost
==================

This example fits an AdaBoosted decision stump on a classification dataset and
plots the decision boundary, class probabilities, and continuous two-class
output value.

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
                           n_clusters_per_class=1,
                           n_informative=2,
                           n_redundant=0,
                           random_state=1)

bdt = AdaBoostClassifier(DecisionTreeClassifier(max_depth=1), real=False)

bdt.fit(X, y)

plot_colors = "br"
plot_step = 0.02
class_names = "AB"

pl.figure(figsize=(15, 5))

# Plot the decision boundaries
pl.subplot(131)
x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, plot_step),
                     np.arange(y_min, y_max, plot_step))

Z = bdt.predict(np.c_[xx.ravel(), yy.ravel()])
Z = Z.reshape(xx.shape)
cs = pl.contourf(xx, yy, Z, cmap=pl.cm.Paired)
pl.axis("tight")

# Plot the training points
for i, n, c in zip(xrange(2), class_names, plot_colors):
    idx = np.where(y == i)
    pl.scatter(X[idx, 0], X[idx, 1],
               c=c, cmap=pl.cm.Paired,
               label="Class %s" % n)
pl.axis("tight")
pl.legend(loc='upper right')
pl.xlabel("Decision Boundary")

# Plot the class probabilities
class_proba = bdt.predict_proba(X)[:, 0]
pl.subplot(132)
for i, n, c in zip(xrange(2), class_names, plot_colors):
    pl.hist(class_proba[y == i],
            bins=20,
            range=(0, 1),
            facecolor=c,
            label='Class %s' % n)
pl.legend(loc='upper center')
pl.ylabel('Samples')
pl.xlabel('Class Probability')

# Plot the two-class output
twoclass_output = bdt.predict_twoclass(X)
pl.subplot(133)
for i, n, c in zip(xrange(2), class_names, plot_colors):
    pl.hist(twoclass_output[y == i],
            bins=20,
            range=(0, 1),
            facecolor=c,
            label='Class %s' % n)
pl.legend(loc='upper right')
pl.ylabel('Samples')
pl.xlabel('Two-Class Output')

pl.subplots_adjust(wspace=0.25)
pl.show()
