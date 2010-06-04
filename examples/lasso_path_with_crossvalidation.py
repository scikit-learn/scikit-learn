"""
=========================================================
Lasso parameter estimation with path and cross-validation
=========================================================

"""

import numpy as np

################################################################################
# generate some sparse data to play with

n_samples, n_features = 60, 100

np.random.seed(1)
X = np.random.randn(n_samples, n_features)
coef = 3*np.random.randn(n_features)
coef[10:] = 0 # sparsify coef
y = np.dot(X, coef)

# add noise
y += 0.01 * np.random.normal((n_samples,))

# Split data in train set and test set
X_train, y_train = X[:n_samples/2], y[:n_samples/2]
X_test, y_test = X[n_samples/2:], y[n_samples/2:]


################################################################################
# Lasso with path and cross-validation using optimized_lasso
from scikits.learn.cross_val import KFold
from scikits.learn.glm.coordinate_descent import optimized_lasso

# Instanciate cross-validation generator
cv = KFold(n_samples/2, 5)

# Estimate optimized lasso model
lasso_opt = optimized_lasso(X_train, y_train, cv, n_alphas=100, eps=1e-3, maxit=100)
y_ = lasso_opt.predict(X_test)

print lasso_opt

# Compute explained variance on test data
print "r^2 on test data : %f" % (1 - np.linalg.norm(y_test - y_)**2
                                      / np.linalg.norm(y_test)**2)

################################################################################
# Lasso with path and cross-validation using LassoPath path
from scikits.learn.glm.coordinate_descent import LassoPath
lasso_path = LassoPath()

y_pred = lasso_path.fit(X_train, y_train).predict(X_test)

print lasso_path

# Compute explained variance on test data
print "r^2 on test data : %f" % (1 - np.linalg.norm(y_test - y_)**2
                                      / np.linalg.norm(y_test)**2)

