"""
========================================
Testing and Training Error with Boosting
========================================

This example shows the use of boosting to improve prediction accuracy. The
error on the test and training sets after each boosting iteration is plotted on
the left. The classification error of each tree is shown on the right.
"""
print __doc__

from itertools import izip

import pylab as pl

from sklearn.ensemble import AdaBoostClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.datasets.samples_generator import make_gaussian_quantiles

X, y = make_gaussian_quantiles(n_samples=13000, n_features=10,
                               n_classes=3, random_state=1)

n_split = 3000

X_train, X_test = X[:n_split], X[n_split:]
y_train, y_test = y[:n_split], y[n_split:]

test_errors = []
train_errors = []

bdt = AdaBoostClassifier(
    DecisionTreeClassifier(max_depth=2),
    n_estimators=600,
    learning_rate=1)

bdt.fit(X_train, y_train)


for y_test_predict, y_train_predict in izip(bdt.staged_predict(X_test),
                                            bdt.staged_predict(X_train)):
    test_errors.append((y_test_predict != y_test).sum() /
                       float(y_test.shape[0]))
    train_errors.append((y_train_predict != y_train).sum() /
                        float(y_train.shape[0]))

n_trees = xrange(1, len(bdt) + 1)

pl.figure(figsize=(10, 5))

pl.subplot(1, 2, 1)
pl.plot(n_trees, test_errors, "b", label='test')
pl.plot(n_trees, train_errors, "r", label='train')
pl.legend()
pl.ylim(0., 0.62)
pl.ylabel('Error')
pl.xlabel('Number of Trees')

pl.subplot(1, 2, 2)
pl.plot(n_trees, bdt.errors_, "b")
pl.ylabel('Error')
pl.xlabel('Tree')
pl.ylim((.2, max(bdt.errors_) * 1.2))
pl.xlim((-20, len(bdt) + 20))

# prevent overlapping y-axis labels
pl.subplots_adjust(wspace=0.25)
pl.show()
