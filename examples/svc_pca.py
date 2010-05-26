# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

import pylab as pl
import numpy as np

from scikits.learn.svm import SVC
from scikits.learn import datasets

iris = datasets.load_iris()
y = iris.target

# Some noisy data not correlated
np.random.seed(0)
E = np.random.normal(size=(len(iris.data), 200))

# Add the noisy data to the informative features
X = np.hstack((iris.data, E))

n_features = X.shape[1]

C = 1.0

n_samples = X.shape[0]

classifier = SVC(kernel='linear', C=C)

import random
random.seed(0)
order = range(n_samples)
random.shuffle(order)
X = X[order]
y = y[order]

# split data in 2 folds
X_train, X_test = X[:n_samples/2], X[n_samples/2:]
y_train, y_test = y[:n_samples/2], y[n_samples/2:]

print "** no PCA"
y_pred = classifier.fit(X_train, y_train).predict(X_test)

classif_rate = np.mean(y_pred.ravel() == y_test.ravel()) * 100
print  "classif_rate : %f " % classif_rate

# print "** with PCA"
# # Do PCA
# import scipy.linalg as linalg
# U, s, Vh = linalg.svd(X, full_matrices=False)
# K = 10
# X = np.dot(U[:,:K], np.diag(s[:K]))
# 
# print X.shape
# y_pred = classifier.fit(X_train, y_train).predict(X_test)
# 
# classif_rate = np.mean(y_pred.ravel() == y_test.ravel()) * 100
# print  "classif_rate : %f " % classif_rate

print "** with Precomputed kernel"
# Do PCA
classifier = SVC(kernel='precomputed', C=C)

y_pred = classifier.fit(np.dot(X_train, X_train.T), y_train).predict(X_test)

classif_rate = np.mean(y_pred.ravel() == y_test.ravel()) * 100
print  "classif_rate : %f " % classif_rate
