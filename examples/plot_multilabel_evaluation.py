# Authors: Arnaud Joly
# License: BSD 3 clause
"""
=================================
Multilabel classifier and metrics
=================================


TODO:
    - Note on multiclass vs multilabel => explosion of label 2 ** n_classes

"""
import numpy as np
import matplotlib.pylab as pl

from sklearn.datasets import make_multilabel_classification
from sklearn.cross_validation import train_test_split
from sklearn.datasets import make_multilabel_classification
from sklearn.decomposition import PCA

from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier

X, y = make_multilabel_classification(return_indicator=True,
                                      n_features=20,
                                      n_samples=50,
                                      n_classes=2,
                                      random_state=5)


def plot_subfigure(X, y, subplot, title, transform):
    X = PCA(n_components=2).fit_transform(X)
    min_x = np.min(X[:, 0])
    max_x = np.max(X[:, 0])

    pl.subplot(2, 2, subplot)
    pl.title(title)

    zero_class = y[:, 0] == 1
    one_class = y[:, 1] == 1
    pl.scatter(X[:, 0], X[:, 1], s=40, c='gray')
    pl.scatter(X[zero_class, 0], X[zero_class, 1], s=160, edgecolors='b',
               facecolors='none', linewidths=2, label='Class 1')
    pl.scatter(X[one_class, 0], X[one_class, 1], s=80, edgecolors='orange',
               facecolors='none', linewidths=2, label='Class 2')
    pl.axis('tight')

    pl.xticks(())
    pl.yticks(())

    pl.savefig("temp.png")

plot_subfigure(X, y, 1, "test", "pca")
