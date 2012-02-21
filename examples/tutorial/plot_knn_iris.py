#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
=========================================================
KNN (k-nearest neighbors) classification example
=========================================================
The K-Nearest-Neighbors algorithm is used below as a
classification tool. The data set 
(`Iris <http://en.wikipedia.org/wiki/Iris_flower_data_set>`_)
is first cut into a training set, as done in the interactive
python code.

The decision boundaries, whom are obtained by using the training
set, are shown with all the points in the training-set.


"""
print __doc__


# Code source: Gael Varoqueux
# Modified for Documentation merge by Jaques Grobler
# License: BSD

import numpy as np
import pylab as pl
from sklearn import neighbors, datasets

# import some data to play with
iris = datasets.load_iris()
X = iris.data[:, :2] # we only take the first two features. 
Y = iris.target

# Split iris data in train and test data
# A random permutation, to split the data randomly
np.random.seed(0)
indices = np.random.permutation(len(X))
X_train = X[indices[:-10]]
Y_train = Y[indices[:-10]]
X_test = X[indices[-10:]]
Y_test = Y[indices[-10:]]

#The splitting above serves no purpose in this program
#it's just to follow the ipython example in the tutorial.

h = .02 # step size in the mesh

knn=neighbors.KNeighborsClassifier()

# we create an instance of Neighbours Classifier and fit the data.
knn.fit(X_train, Y_train)

# Plot the decision boundary. For that, we will asign a color to each
# point in the mesh [x_min, m_max]x[y_min, y_max].
x_min, x_max = X_test[:,0].min() - .5, X_test[:,0].max() + .5
y_min, y_max = X_test[:,1].min() - .5, X_test[:,1].max() + .5
xx, yy = np.meshgrid(np.arange(x_min, x_max, h), np.arange(y_min, y_max, h))
Z = knn.predict(np.c_[xx.ravel(), yy.ravel()])

# Put the result into a color plot
Z = Z.reshape(xx.shape)
pl.figure(1, figsize=(4, 3))
pl.set_cmap(pl.cm.Paired)
pl.pcolormesh(xx, yy, Z)

# Plot also the training points
pl.scatter(X_train[:,0], X_train[:,1],facecolors='none', edgecolors='k' )
pl.xlabel('Sepal length')
pl.ylabel('Sepal width')

pl.xlim(xx.min(), xx.max())
pl.ylim(yy.min(), yy.max())
pl.xticks(())
pl.yticks(())

pl.show()
