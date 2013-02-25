#!/usr/bin/python
# -*- coding: utf-8 -*-

"""SimpleELM Examples"""

# Author: David C. Lambert <dcl@panix.com>
# License: Simple BSD

from sklearn.cross_validation import train_test_split
from sklearn.datasets import make_regression, load_digits

from sklearn.elm import SimpleELMRegressor, SimpleELMClassifier

mrx, mry = make_regression(n_samples=2000, n_targets=1, random_state=0)
X_train, X_test, y_train, y_test = train_test_split(mrx, mry, test_size=0.3, random_state=0)
elmr = SimpleELMRegressor(n_hidden=200, random_state=0)
elmr.fit(X_train, y_train)

print
print "SimpleELMRegressor(n_hidden=200) on make_regression(n_samples=2000, n_targets=1, random_state=0)"
print "TrainAcc: %.3f TestAcc: %.3f" %  (elmr.score(X_train, y_train), elmr.score(X_test, y_test))
print

digits = load_digits()
digits_X, digits_y = digits.data, digits.target
X_train, X_test, y_train, y_test = train_test_split(digits_X, digits_y, test_size=0.3, random_state=0)
elmc = SimpleELMClassifier(n_hidden=200, random_state=0)
elmc.fit(X_train, y_train)

print "SimpleELMClassifier(n_hidden=500) on digits data"
print "TrainAcc: %.3f TestAcc: %.3f" %  (elmc.score(X_train, y_train), elmc.score(X_test, y_test))
print
