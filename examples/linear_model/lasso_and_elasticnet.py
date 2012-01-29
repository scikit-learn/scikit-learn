"""
========================
Lasso regression example
========================

"""
print __doc__

import numpy as np

###############################################################################
# generate some sparse data to play with

n_samples, n_features = 50, 200
X = np.random.randn(n_samples, n_features)
coef = 3 * np.random.randn(n_features)
coef[10:] = 0  # sparsify coef
y = np.dot(X, coef)

# add noise
y += 0.01 * np.random.normal((n_samples,))

# Split data in train set and test set
n_samples = X.shape[0]
X_train, y_train = X[:n_samples / 2], y[:n_samples / 2]
X_test, y_test = X[n_samples / 2:], y[n_samples / 2:]

###############################################################################
# Lasso
from sklearn.linear_model import Lasso

alpha = 0.1
lasso = Lasso(alpha=alpha)

y_pred_lasso = lasso.fit(X_train, y_train).predict(X_test)
print lasso
print "r^2 on test data : %f" % (1 - np.linalg.norm(y_test - y_pred_lasso) ** 2
                                      / np.linalg.norm(y_test) ** 2)

###############################################################################
# ElasticNet
from sklearn.linear_model import ElasticNet

enet = ElasticNet(alpha=alpha, rho=0.7)

y_pred_enet = enet.fit(X_train, y_train).predict(X_test)
print enet
print "r^2 on test data : %f" % (1 - np.linalg.norm(y_test - y_pred_enet) ** 2
                                      / np.linalg.norm(y_test) ** 2)
