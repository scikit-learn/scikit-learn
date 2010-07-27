"""
Parameter estimation using grid search with a nested cross-validation
=======================================================================

The classifier is optimized by "nested" cross-validation using the
GridSearchCV object.

The performance of the selected parameters is evaluated using
cross-validation (different than the nested crossvalidation that is used
to select the best classifier). 

"""

import numpy as np
from scikits.learn.svm import SVC
from scikits.learn.cross_val import StratifiedKFold, GridSearchCV
from scikits.learn import datasets

# The Digits dataset
digits = datasets.load_digits()

# To apply an classifier on this data, we need to flatten the image, to
# turn the data in a (samples, feature) matrix:
n_samples = len(digits.images)
X = digits.images.reshape((n_samples, -1))
y = digits.target

parameters = {'kernel':('linear', 'rbf'), 'C':[0.1, 1]}

def loss_func(y1, y2):
    return np.mean(y1 != y2)

clf = GridSearchCV(SVC(), parameters, loss_func, n_jobs=2)

n_samples, n_features = X.shape
y_pred = []
y_true = []
for train, test in StratifiedKFold(y, 2):
    cv = StratifiedKFold(y[train], 2)
    y_pred.append(clf.fit(X[train], y[train], cv=cv).predict(X[test]))
    y_true.append(y[test])

y_pred = np.concatenate(y_pred)
y_true = np.concatenate(y_true)
classif_rate = np.mean(y_pred == y_true) * 100
print "Classification rate : %f" % classif_rate
