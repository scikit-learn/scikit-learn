"""
===================================================
Recursive feature elimination with cross-validation
===================================================

Recursive feature elimination with automatic tuning of the
number of features selected with cross-validation
"""
print __doc__

from scikits.learn.svm import SVC
from scikits.learn.cross_val import StratifiedKFold
from scikits.learn.feature_selection import RFECV
from scikits.learn.datasets import samples_generator
from scikits.learn.metrics import zero_one

################################################################################
# Loading a dataset

X, y = samples_generator.test_dataset_classif(n_features=500, k=5, seed=0)

################################################################################
# Create the RFE object and compute a cross-validated score

svc = SVC(kernel='linear')
rfecv = RFECV(estimator=svc, n_features=2, percentage=0.1, loss_func=zero_one)
rfecv.fit(X, y, cv=StratifiedKFold(y, 2))

print 'Optimal number of features : %d' % rfecv.support_.sum()

import pylab as pl
pl.figure()
pl.plot(rfecv.cv_scores_)
pl.show()

