"""
===============================
Gradient Boosting OOB estimates
===============================


"""
print(__doc__)

# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD 3 clause

import numpy as np
import pylab as plt

from scipy.stats.distributions import binom

from sklearn import ensemble

###############################################################################
# Generate data

n = 1000
rs = np.random.RandomState(13)
x1 = rs.uniform(size=n)
x2 = rs.uniform(size=n)
x3 = rs.randint(0, 4, size=n)

p = 1 / (1.0 + np.exp(-(np.sin(3 * x1) - 4 * x2 + x3)))
y = binom.rvs(1, p, size=n)
X = np.c_[x1, x2, x3]

X = X.astype(np.float32)
offset = int(X.shape[0] * 0.8)
X_train, y_train = X[:offset], y[:offset]
X_test, y_test = X[offset:], y[offset:]

###############################################################################
# Fit regression model
params = {'n_estimators': 3000, 'max_depth': 4, 'subsample': 0.5,
          'learning_rate': 0.01, 'min_samples_leaf': 1, 'random_state': 3}
clf = ensemble.GradientBoostingClassifier(**params)

clf.fit(X_train, y_train)
acc = clf.score(X_test, y_test)
print("ACC: %.4f" % acc)

n_estimators = params['n_estimators']
x = np.arange(n_estimators) + 1

test_score = np.zeros((n_estimators,), dtype=np.float64)

for i, y_pred in enumerate(clf.staged_decision_function(X_test)):
    test_score[i] = clf.loss_(y_test, y_pred)


#plt.plot(x, -np.cumsum(clf.oob_score_))
cumsum = np.cumsum(clf.oob_improvement_)
best_iter = x[np.argmax(cumsum)]
true_best_iter = np.argmin(test_score)
print("best_iter: %d" % best_iter)
print("true_best_iter: %d" % true_best_iter)

plt.plot(x, cumsum, label='oob', color='cyan')
plt.plot(x, test_score, label='test', color='red')

plt.plot([best_iter, best_iter], [0.0, 1.0], color='cyan')
plt.plot([true_best_iter, true_best_iter], [0.0, 1.0], color='red')

plt.legend(loc='upper right')

plt.show()
