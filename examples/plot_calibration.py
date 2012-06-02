"""
================
Calibration plot
================

"""
print __doc__

# Author: Mathieu Blondel <mathieu@mblondel.org>
# License: BSD Style.

from sklearn.datasets import make_classification, load_digits
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import brier_score, calibration_plot

X, y = make_classification(n_samples=5000, random_state=42)
X_train, y_train = X[:1000], y[:1000]
X_test, y_test = X[1000:], y[1000:]

clf = LogisticRegression()
clf.fit(X_train, y_train)
prob_pos = clf.predict_proba(X_test)[:, 1]
print "Brier scores:"
print "LogisticRegression", brier_score(y_test, prob_pos)
pt_lr, pp_lr = calibration_plot(y_test, prob_pos)

clf = SGDClassifier(loss="log", seed=0)
clf.fit(X_train, y_train)
prob_pos = clf.predict_proba(X_test)
print "SGDClassifier (log)", brier_score(y_test, prob_pos)

pt_sgd, pp_sgd = calibration_plot(y_test, prob_pos)

clf = SGDClassifier(loss="modified_huber", seed=0)
clf.fit(X_train, y_train)
prob_pos = clf.predict_proba(X_test)
print "SGDClassifier (modified_huber)", brier_score(y_test, prob_pos)

pt_sgd2, pp_sgd2 = calibration_plot(y_test, prob_pos)

import pylab as pl

pl.figure()
pl.xlabel("Predicted probability")
pl.ylabel("True probability")
pl.plot([0, 1], [0, 1], "b", label="Perfectly calibrated")
pl.plot(pp_sgd, pt_sgd, "rs-", label="SGDClassifier (log)")
pl.plot(pp_sgd2, pt_sgd2, "gs-", label="SGDClassifier (modified_huber)")
pl.plot(pp_lr, pt_lr, "ms-", label="LogisticRegression")
pl.legend(loc="lower right")
pl.show()
