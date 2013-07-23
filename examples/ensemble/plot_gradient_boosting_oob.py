"""
===============================
Gradient Boosting OOB estimates
===============================

Out-of-bag estimates can be used to estimate the "optimal" number of
boosting iterations. To compute OOB estimates the ``subsample``
argument has to be smaller than 1.
The argument ``oob_improvement_[i]`` holds the improvement in loss on
the out-of-bag samples for the i-th iteration.

The generated plot shows the cumulative sum of the negative OOB improvements
as a function of the boosting iteration. This function should track the
loss on the test set.
The plot also shows the performance of 3-fold cross validation which
usually gives a better estimate but is computationally more demanding.
"""
print(__doc__)

# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD 3 clause

import numpy as np
import pylab as pl

from sklearn import ensemble
from sklearn.cross_validation import KFold

###############################################################################
# Generate data (adapted from G. Ridgeway's gbm example)

n = 1000
rs = np.random.RandomState(13)
x1 = rs.uniform(size=n)
x2 = rs.uniform(size=n)
x3 = rs.randint(0, 4, size=n)

p = 1 / (1.0 + np.exp(-(np.sin(3 * x1) - 4 * x2 + x3)))
y = rs.binomial(1, p, size=n)

X = np.c_[x1, x2, x3]

X = X.astype(np.float32)
offset = int(X.shape[0] * 0.8)
X_train, y_train = X[:offset], y[:offset]
X_test, y_test = X[offset:], y[offset:]

###############################################################################
# Fit regression model
params = {'n_estimators': 2000, 'max_depth': 4, 'subsample': 0.5,
          'learning_rate': 0.01, 'min_samples_leaf': 1, 'random_state': 3}
clf = ensemble.GradientBoostingClassifier(**params)

clf.fit(X_train, y_train)
acc = clf.score(X_test, y_test)
print("ACC: %.4f" % acc)

n_estimators = params['n_estimators']
x = np.arange(n_estimators) + 1


def heldout_score(clf, X_test, y_test):
    score = np.zeros((n_estimators,), dtype=np.float64)
    for i, y_pred in enumerate(clf.staged_decision_function(X_test)):
        score[i] = clf.loss_(y_test, y_pred)
    return score


def cv_estimate(n_folds=3):
    cv = KFold(n=X_train.shape[0], n_folds=n_folds)
    cv_clf = ensemble.GradientBoostingClassifier(**params)
    val_scores = np.zeros((n_estimators,), dtype=np.float64)
    for train, test in cv:
        cv_clf.fit(X_train[train], y_train[train])
        val_scores += heldout_score(cv_clf, X_train[test], y_train[test])
    val_scores /= n_folds
    return val_scores


test_score = heldout_score(clf, X_test, y_test)
cv_score = cv_estimate(3)


#pl.plot(x, -np.cumsum(clf.oob_score_))
cumsum = -np.cumsum(clf.oob_improvement_)
oob_best_iter = x[np.argmin(cumsum)]
test_score -= test_score[0]
test_best_iter = np.argmin(test_score)
cv_score -= cv_score[0]
cv_best_iter = np.argmin(cv_score)

oob_color = map(lambda x: x / 256.0, (190, 174, 212))
test_color = map(lambda x: x / 256.0, (127, 201, 127))
cv_color = map(lambda x: x / 256.0, (253, 192, 134))

pl.plot(x, cumsum, label='OOB loss', color=oob_color)
pl.plot(x, test_score, label='Test loss', color=test_color)
pl.plot(x, cv_score, label='CV loss', color=cv_color)
pl.axvline(x=oob_best_iter, color=oob_color)
pl.axvline(x=test_best_iter, color=test_color)
pl.axvline(x=cv_best_iter, color=cv_color)

pl.legend(loc='upper right')
pl.ylabel('normalized loss')

pl.show()
