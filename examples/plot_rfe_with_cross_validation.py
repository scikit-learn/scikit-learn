"""
===================================================
Recursive feature elimination with cross-validation
===================================================

A recursive feature elimination example with automatic tuning of the
number of features selected with cross-validation.
"""
print __doc__

import numpy as np
from sklearn.svm import SVC
from sklearn.cross_validation import StratifiedKFold
from sklearn.feature_selection import RFECV
from sklearn.datasets import make_classification
from sklearn.metrics import zero_one

# Build a classification task using 3 informative features
X, y = make_classification(n_samples=1000,
                           n_features=25,
                           n_informative=3,
                           n_redundant=2,
                           n_repeated=0,
                           n_classes=8,
                           n_clusters_per_class=1,
                           random_state=0)

# Create the RFE object and compute a cross-validated score.
svc = SVC(kernel="linear")
rfecv = RFECV(estimator=svc,
              step=1,
              cv=StratifiedKFold(y, 2),
              loss_func=zero_one)
rfecv.fit(X, y)

print "Optimal number of features : %d" % rfecv.n_features_

# Plot number of features VS. cross-validation scores
import pylab as pl
pl.figure()
pl.xlabel("Number of features selected")
pl.ylabel("Cross validation score (nb of misclassifications)")
pl.plot(xrange(1, len(rfecv.cv_scores_) + 1), rfecv.cv_scores_)
pl.show()
