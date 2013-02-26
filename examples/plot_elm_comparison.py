#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
======================
ELM Classifiers Comparison
======================
A comparison of a several ELMClassifiers with different types of hidden
layer activations.

ELMClassifier is a classifier based on the Extreme Learning Machine,
a single layer feedforward network with random hidden layer components
and least squares fitting of the hidden->output weights by default [1][2]

The point of this example is to illustrate the nature of decision boundaries
with different hidden layer activation types and regressors.

This should be taken with a grain of salt, as the intuition conveyed by
these examples does not necessarily carry over to real datasets.

In particular in high dimensional spaces data can more easily be separated
linearly and the simplicity of classifiers such as naive Bayes and linear SVMs
might lead to better generalization.

The plots show training points in solid colors and testing points
semi-transparent. The lower right shows the classification accuracy on the test
set.

References
__________
.. [1] http://www.extreme-learning-machines.org
.. [2] G.-B. Huang, Q.-Y. Zhu and C.-K. Siew, "Extreme Learning Machine:
          Theory and Applications", Neurocomputing, vol. 70, pp. 489-501,
          2006.

===============================================================================
Basis Functions:
  rbf = exp(-gamma * (||x-c||/r)^2)
  tanh = np.tanh
  sinsq = (lambda x: np.power(np.sin(x), 2.0))
  tribas = (lambda x: np.clip(1.0 - np.fabs(x), 0.0, 1.0))
  hardlim = (lambda x: np.array(x > 0.0, dtype=float))

Label Legend:
  ELM(10,tanh)      :10 tanh units
  ELM(10,tanh,LR)   :10 tanh units, LogisticRegression
  ELM(10,sinsq)     :10 sin*sin units
  ELM(10,tribas)    :10 tribas units
  ELM(10,hardlim)   :10 hardlim units
  ELM(20,rbf(0.1))  :20 rbf units gamma=0.1

"""
print __doc__


# Code source: Gael Varoqueux
#              Andreas Mueller
# Modified for Documentation merge by Jaques Grobler
# Modified for Extreme Learning Machine Classifiers by David Lambert
# License: BSD

import numpy as np
import pylab as pl

from matplotlib.colors import ListedColormap
from sklearn.datasets import make_moons, make_circles, make_classification
from sklearn.preprocessing import StandardScaler
from sklearn.cross_validation import train_test_split
from sklearn.linear_model import LogisticRegression

from sklearn.elm import ELMClassifier
from sklearn.random_hidden_layer import (SimpleRandomHiddenLayer,
                                         RBFRandomHiddenLayer)


def get_data_bounds(X):
    h = .02  # step size in the mesh

    x_min, x_max = X[:, 0].min() - .5, X[:, 0].max() + .5
    y_min, y_max = X[:, 1].min() - .5, X[:, 1].max() + .5

    xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
                         np.arange(y_min, y_max, h))

    return (x_min, x_max, y_min, y_max, xx, yy)


def plot_data(ax, X_train, y_train, X_test, y_test, xx, yy):
    cm = ListedColormap(['#FF0000', '#0000FF'])
    # Plot the training points
    ax.scatter(X_train[:, 0], X_train[:, 1], c=y_train, cmap=cm)
    # and testing points
    ax.scatter(X_test[:, 0], X_test[:, 1], c=y_test, cmap=cm, alpha=0.6)
    ax.set_xlim(xx.min(), xx.max())
    ax.set_ylim(yy.min(), yy.max())
    ax.set_xticks(())
    ax.set_yticks(())


def plot_contour(ax, X_train, y_train, X_test, y_test, xx, yy, Z):
    cm = pl.cm.RdBu
    cm_bright = ListedColormap(['#FF0000', '#0000FF'])

    ax.contourf(xx, yy, Z, cmap=cm, alpha=.8)

    # Plot also the training points
    ax.scatter(X_train[:, 0], X_train[:, 1], c=y_train, cmap=cm_bright)
    # and testing points
    ax.scatter(X_test[:, 0], X_test[:, 1], c=y_test, cmap=cm_bright, alpha=0.6)

    ax.set_xlim(xx.min(), xx.max())
    ax.set_ylim(yy.min(), yy.max())
    ax.set_xticks(())
    ax.set_yticks(())

    ax.set_title(name)
    ax.text(xx.max() - 0.3, yy.min() + 0.3, ('%.2f' % score).lstrip('0'),
            size=13, horizontalalignment='right')


def make_datasets():
    return [make_moons(n_samples=200, noise=0.3, random_state=0),
            make_circles(n_samples=200, noise=0.2, factor=0.5, random_state=1),
            make_linearly_separable()]


def make_classifiers():
    sinsq = (lambda x: np.power(np.sin(x), 2.0))
    tribas = (lambda x: np.clip(1.0 - np.fabs(x), 0.0, 1.0))
    hardlim = (lambda x: np.array(x > 0.0, dtype=float))

    names = ["ELM(10,tanh)", "ELM(10,tanh,LR)", "ELM(10,sinsq)",
             "ELM(10,tribas)", "ELM(hardlim)", "ELM(20,rbf(0.1))"]

    nh = 10
    srhl_rbf = RBFRandomHiddenLayer(n_hidden=nh*2, gamma=0.1, random_state=0)
    srhl_tanh = SimpleRandomHiddenLayer(n_hidden=nh, user_func=np.tanh, random_state=0)
    srhl_sinsq = SimpleRandomHiddenLayer(n_hidden=nh, user_func=sinsq, random_state=0)
    srhl_tribas = SimpleRandomHiddenLayer(n_hidden=nh, user_func=tribas, random_state=0)
    srhl_hardlim = SimpleRandomHiddenLayer(n_hidden=nh, user_func=hardlim, random_state=0)

    log_reg = LogisticRegression()

    classifiers = [ELMClassifier(srhl_tanh),
                   ELMClassifier(srhl_tanh, regressor=log_reg),
                   ELMClassifier(srhl_sinsq),
                   ELMClassifier(srhl_tribas),
                   ELMClassifier(srhl_hardlim),
                   ELMClassifier(srhl_rbf)]

    return names, classifiers


def make_linearly_separable():
    X, y = make_classification(n_samples=200, n_features=2, n_redundant=0,
                               n_informative=2, random_state=1,
                               n_clusters_per_class=1)
    rng = np.random.RandomState(2)
    X += 2 * rng.uniform(size=X.shape)
    return (X, y)

###############################################################################

datasets = make_datasets()
names, classifiers = make_classifiers()

i = 1
figure = pl.figure(figsize=(18, 9))

# iterate over datasets
for ds in datasets:
    # preprocess dataset, split into training and test part
    X, y = ds
    X = StandardScaler().fit_transform(X)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.4, random_state=0)

    x_min, x_max, y_min, y_max, xx, yy = get_data_bounds(X)

    # plot dataset first
    ax = pl.subplot(len(datasets), len(classifiers) + 1, i)
    plot_data(ax, X_train, y_train, X_test, y_test, xx, yy)

    i += 1

    # iterate over classifiers
    for name, clf in zip(names, classifiers):
        ax = pl.subplot(len(datasets), len(classifiers) + 1, i)
        clf.fit(X_train, y_train)
        score = clf.score(X_test, y_test)

        # Plot the decision boundary. For that, we will asign a color to each
        # point in the mesh [x_min, m_max]x[y_min, y_max].
        Z = clf.decision_function(np.c_[xx.ravel(), yy.ravel()])

        # Put the result into a color plot
        Z = Z.reshape(xx.shape)

        plot_contour(ax, X_train, y_train, X_test, y_test, xx, yy, Z)

        i += 1

figure.subplots_adjust(left=.02, right=.98)
pl.show()
