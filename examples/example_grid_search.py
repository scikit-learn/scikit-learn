"""Run parameter estimation using grid search
in a nested cross-validation setting.
"""

import numpy as np
from scikits.learn.svm import SVC
from scikits.learn.cross_val import LeaveOneOut
from scikits.learn import datasets
from scikits.learn.grid_search import GridSearch

# The IRIS dataset
iris = datasets.load_iris()

X = iris.data
y = iris.target

parameters = {'kernel':('linear', 'rbf'), 'C':[0.1, 1]}

def loss_func(y1, y2):
    return np.mean(y1 != y2)

clf = GridSearch(SVC, parameters, LeaveOneOut, loss_func, n_jobs=2)

"""
Run crossvalidation (different than the nested crossvalidation that is used
to select the best classifier). The classifier is optimized by "nested"
crossvalidation
"""
n_samples, n_features = X.shape
y_pred = np.zeros_like(y)
for train, test in LeaveOneOut(n_samples):
    y_pred[test] = clf.fit(X[train], y[train]).predict(X[test]).astype(np.int)

classif_rate = np.mean(y_pred == y) * 100
print "Classification rate : %f" % classif_rate
