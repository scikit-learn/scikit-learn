"""
Parameter estimation using grid search with a nested cross-validation
=======================================================================

The classifier is optimized by "nested" cross-validation using the
GridSearchCV object.

The performance of the selected parameters is evaluated using
cross-validation (different than the nested cross-validation that is used
to select the best classifier). 

"""

import numpy as np
from scikits.learn.svm import SVC
from scikits.learn.cross_val import StratifiedKFold
from scikits.learn.grid_search import GridSearchCV
from scikits.learn import datasets
from scikits.learn.metrics import zero_one

################################################################################
# Loading the Digits dataset
digits = datasets.load_digits()

# To apply an classifier on this data, we need to flatten the image, to
# turn the data in a (samples, feature) matrix:
n_samples = len(digits.images)
X = digits.images.reshape((n_samples, -1))
y = digits.target

################################################################################
# Set the parameters by cross-validation
tuned_parameters = [{'kernel':('rbf', ), 'gamma':[1e-3, 1e-4]},
                    {'kernel':('linear', )}]

clf = GridSearchCV(SVC(C=1), tuned_parameters, n_jobs=2)

y_pred = []
y_true = []
for train, test in StratifiedKFold(y, 2):
    clf.fit(X[train], y[train], cv=StratifiedKFold(y[train], 5))
    y_pred = np.append(y_pred, clf.predict(X[test]))
    y_true = np.append(y_true, y[test])

classif_rate = np.mean(y_pred == y_true) * 100
print "Classification rate : %f" % classif_rate
