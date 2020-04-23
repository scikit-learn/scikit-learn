#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
=========================================================
Logistic Regression 3-class Classifier
=========================================================

Show below is a logistic-regression classifiers decision boundaries on the
first two dimensions (sepal length and width) of the `iris
<https://en.wikipedia.org/wiki/Iris_flower_data_set>`_ dataset. The datapoints
are colored according to their labels.

"""
print(__doc__)

# Code source: Gaël Varoquaux
# Modified for documentation by Jaques Grobler
# License: BSD 3 clause

import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn import datasets
from sklearn.inspection import plot_decision_boundary

# import some data to play with
iris = datasets.load_iris()
X = iris.data[:, :2]  # we only take the first two features.
Y = iris.target

logreg = LogisticRegression(C=1e5)

# Create an instance of Logistic Regression Classifier and fit the data.
logreg.fit(X, Y)

_, ax = plt.subplots(figsize=(4, 3))
plot_decision_boundary(logreg, X, cmap=plt.cm.Paired, ax=ax,
                       response_method='predict', plot_method='pcolormesh')

# Plot also the training points
plt.scatter(X[:, 0], X[:, 1], c=Y, edgecolors='k', cmap=plt.cm.Paired)
plt.xlabel('Sepal length')
plt.ylabel('Sepal width')

plt.xticks(())
plt.yticks(())

plt.show()
