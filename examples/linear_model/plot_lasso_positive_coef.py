"""
===========================================
Lasso with positive constraint coefficients
===========================================

"""
print __doc__

import numpy as np
import pylab as pl

from sklearn.linear_model import Lasso
from sklearn.metrics import r2_score
from sklearn.datasets.samples_generator import make_regression
from sklearn.cross_validation import train_test_split

###############################################################################
# generate some sparse data to play with

X, y, coef = make_regression(n_samples=100, n_features=1000, noise=0.5,
                     coef=True, random_state=3)

X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

###############################################################################
# Lasso
alpha = 0.1
lasso = Lasso(alpha=alpha, max_iter=1000)

y_pred_lasso = \
    lasso.fit(X_train, y_train, positive_constraint=False).predict(X_test)

r2_score_lasso = r2_score(y_test, y_pred_lasso)
print lasso
print "r^2 on test data : %f" % r2_score_lasso

###############################################################################
# positive constraint Lasso
alpha = 0.1
pos_lasso = Lasso(alpha=alpha, max_iter=1000)

pos_y_pred_lasso = \
    pos_lasso.fit(X_train, y_train, positive_constraint=True).predict(X_test)

pos_r2_score_lasso = r2_score(y_test, pos_y_pred_lasso)
print pos_lasso
print "r^2 on test data : %f" % pos_r2_score_lasso


pl.plot(pos_lasso.coef_, label='positive Lasso coefficients')
pl.plot(lasso.coef_, label='Lasso coefficients')
pl.plot(coef, '--', label='original coefficients')
pl.legend(loc='best')
pl.title("Lasso R^2: %f, positive Lasso R^2: %f" % (r2_score_lasso,
    pos_r2_score_lasso))
pl.show()
