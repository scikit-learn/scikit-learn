"""
=================================================
Plot multi-class DecisionTree on the iris dataset
=================================================

Plot decision surface of multi-class DecisionTree on iris dataset on pairwise
selection of features.
 
"""
print __doc__

import numpy as np
import pylab as pl
from scikits.learn import datasets
from scikits.learn.tree_model import DecisionTreeClassifier

# import some data to play with
iris = datasets.load_iris()
for pairidx, pair in enumerate([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]]):
    X = iris.data[:, pair]  # we only take the first two features. We could
    # avoid this ugly slicing by using a two-dim dataset
    y = iris.target
    K=3
    colors = "bry"

    # shuffle
    idx = np.arange(X.shape[0])
    np.random.seed(13)
    np.random.shuffle(idx)
    X = X[idx]
    y = y[idx]

    # standardize
    mean = X.mean(axis=0)
    std = X.std(axis=0)
    X = (X - mean) / std

    h = .02  # step size in the mesh

    clf = DecisionTreeClassifier(K).fit(X, y)

    # create a mesh to plot in
    x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
    y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
                         np.arange(y_min, y_max, h))

    pl.subplot(2,3,pairidx+1)
    pl.set_cmap(pl.cm.Paired)

    # Plot the decision boundary. For that, we will asign a color to each
    # point in the mesh [x_min, m_max]x[y_min, y_max].
    Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
    # Put the result into a color plot
    Z = Z.reshape(xx.shape)
    pl.set_cmap(pl.cm.Paired)
    cs = pl.contourf(xx, yy, Z)
    pl.xlabel(iris.feature_names[pair[0]])
    pl.ylabel(iris.feature_names[pair[1]])
    pl.axis('tight')

    # Plot also the training points
    for i, color in zip(xrange(K), colors):
        idx = np.where(y == i)
        pl.scatter(X[idx, 0], X[idx, 1], c=color, label=iris.target_names[i])
    pl.axis('tight')

pl.suptitle("Decision surface of multi-class decision tree using paired features")

# TODO: Add tree diagram

pl.legend()
pl.show()

