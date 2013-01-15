"""
=====================================
Multi-Class AdaBoosted Decision Trees
=====================================

This example reproduces Figure 1 of Zhu et al [1] and shows how boosting can
improve prediction accuracy on a multi-class problem. The classification
dataset is constructed by taking a ten-dimensional standard normal distribution
and defining three classes separated by nested concentric ten-dimensional
spheres such that roughly equal numbers of samples are in each class (quantiles
of the :math:`\Chi^2` distribution).

The performance of the SAMME and SAMME.R [1] algorithms are compared.
The error of each algorithm on the test set after each boosting iteration is
shown on the left, the classification error on the test set of each tree is
shown in the middle, and the boost weight of each tree is shown on the right.

.. [1] J. Zhu, H. Zou, S. Rosset, T. Hastie, "Multi-class AdaBoost", 2009.

"""
print __doc__

# Author: Noel Dawe <noel.dawe@gmail.com>
#
# License: BSD

from itertools import izip

import pylab as pl

from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.datasets.samples_generator import make_gaussian_quantiles
from sklearn.metrics import accuracy_score


X, y = make_gaussian_quantiles(n_samples=13000, n_features=10,
                               n_classes=3, random_state=1)

n_split = 3000

X_train, X_test = X[:n_split], X[n_split:]
y_train, y_test = y[:n_split], y[n_split:]

bdt_real = AdaBoostClassifier(
    DecisionTreeClassifier(max_depth=2),
    n_estimators=600,
    learning_rate=1)

bdt_discrete = AdaBoostClassifier(
    DecisionTreeClassifier(max_depth=2),
    n_estimators=600,
    learning_rate=1.5,
    real=False)

bdt_real.fit(X_train, y_train)
bdt_discrete.fit(X_train, y_train)

real_test_errors = []
discrete_test_errors = []

for real_test_predict, discrete_train_predict in izip(
        bdt_real.staged_predict(X_test), bdt_discrete.staged_predict(X_test)):
    real_test_errors.append(1. -
        accuracy_score(real_test_predict, y_test))
    discrete_test_errors.append(1. -
        accuracy_score(discrete_train_predict, y_test))

n_trees = xrange(1, len(bdt_discrete) + 1)

pl.figure(figsize=(15, 5))

pl.subplot(131)
pl.plot(n_trees, discrete_test_errors, c='black', label='SAMME')
pl.plot(n_trees, real_test_errors, c='black', linestyle='dashed', label='SAMME.R')
pl.legend()
pl.ylim(0.18, 0.62)
pl.ylabel('Test Error')
pl.xlabel('Number of Trees')

pl.subplot(132)
pl.plot(n_trees, bdt_discrete.errors_, "b", label='SAMME')
pl.plot(n_trees, bdt_real.errors_, "r", label='SAMME.R')
pl.legend()
pl.ylabel('Error')
pl.xlabel('Tree')
pl.ylim((.2, max(bdt_real.errors_.max(), bdt_discrete.errors_.max()) * 1.2))
pl.xlim((-20, len(bdt_discrete) + 20))

pl.subplot(133)
pl.plot(n_trees, bdt_discrete.weights_, "b", label='SAMME')
pl.legend()
pl.ylabel('Weight')
pl.xlabel('Tree')
pl.ylim((0, bdt_discrete.weights_.max() * 1.2))
pl.xlim((-20, len(bdt_discrete) + 20))

# prevent overlapping y-axis labels
pl.subplots_adjust(wspace=0.25)
pl.show()
