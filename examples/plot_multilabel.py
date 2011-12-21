"""
=========================
Multilabel classification
=========================

This example simulates a multi-label document classification problem. The
dataset is generated randomly based on the following process:

    - pick the number of labels: n ~ Poisson(n_labels)
    - n times, choose a class c: c ~ Multinomial(theta)
    - pick the document length: k ~ Poisson(length)
    - k times, choose a word: w ~ Multinomial(theta_c)

In the above process, rejection sampling is used to make sure that
n is never zero or more than 2, and that the document length
is never zero. Likewise, we reject classes which have already been chosen.
The documents that are assigned to both classes are plotted surrounded by
two colored circles.

The classification is performed by projecting to the first two principal
components for visualisation purposes, followed by using the
:class:`sklearn.multiclass.OneVsRestClassifier` metaclassifier using two SVCs
with linear kernels to learn a discriminative model for each class.
"""
print __doc__

import numpy as np
import matplotlib.pylab as pl

from sklearn.datasets import make_multilabel_classification
from sklearn.multiclass import OneVsRestClassifier
from sklearn.svm import SVC
from sklearn.decomposition import PCA


def plot_hyperplane(clf, min_x, max_x, linestyle, label):
    # get the separating hyperplane
    w = clf.coef_[0]
    a = -w[0] / w[1]
    xx = np.linspace(min_x, max_x)
    yy = a * xx - (clf.intercept_[0]) / w[1]
    pl.plot(xx, yy, linestyle, label=label)


X, Y = make_multilabel_classification(n_classes=2, n_labels=1, random_state=42)
X = PCA(n_components=2).fit_transform(X)
min_x = np.min(X[:, 0])
max_x = np.max(X[:, 0])

classif = OneVsRestClassifier(SVC(kernel='linear'))
classif.fit(X, Y)

pl.figure()
pl.title('Multilabel classification example')
pl.xlabel('First principal component')
pl.ylabel('Second principal component')

zero_class = np.where([0 in y for y in Y])
one_class = np.where([1 in y for y in Y])
pl.scatter(X[:, 0], X[:, 1], s=40, c='gray')
pl.scatter(X[zero_class, 0], X[zero_class, 1], s=160, edgecolors='b',
           facecolors='none', linewidths=2, label='Class 1')
pl.scatter(X[one_class, 0], X[one_class, 1], s=80, edgecolors='orange',
           facecolors='none', linewidths=2, label='Class 2')
pl.axis('tight')

plot_hyperplane(classif.estimators_[0], min_x, max_x, 'k--',
                'Boundary\nfor class 1')
plot_hyperplane(classif.estimators_[1], min_x, max_x, 'k-.',
                'Boundary\nfor class 2')

pl.legend()

pl.show()
