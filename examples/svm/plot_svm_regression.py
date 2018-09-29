"""
===================================================================
Support Vector Regression (SVR) using linear and non-linear kernels
===================================================================

Toy example of 1D regression using linear, polynomial and RBF kernels.

"""
print(__doc__)

import numpy as np
from sklearn.svm import SVR
import matplotlib.pyplot as plt

# #############################################################################
# Generate sample data
X = np.sort(5 * np.random.rand(40, 1), axis=0)
y = np.sin(X).ravel()

# #############################################################################
# Add noise to targets
y[::5] += 3 * (0.5 - np.random.rand(8))

# #############################################################################
# Fit regression model
svr_rbf = SVR(kernel='rbf', C=100, gamma=0.1, epsilon=.1)
svr_lin = SVR(kernel='linear', C=100)
svr_poly = SVR(kernel='poly', C=100, degree=3, epsilon=.1, coef0=0)
y_rbf = svr_rbf.fit(X, y).predict(X)
y_lin = svr_lin.fit(X, y).predict(X)
y_poly = svr_poly.fit(X, y).predict(X)

# #############################################################################
# Look at the results
lw = 2

plt.scatter(X[svr_rbf.support_], y[svr_rbf.support_], facecolor="none",
            edgecolor="m", marker='8', label='rbf support vectors', s=50)
plt.scatter(X[svr_lin.support_], y[svr_lin.support_], facecolor="none",
            edgecolor="c", marker='^', label='linear support vectors', s=50)
plt.scatter(X[svr_poly.support_], y[svr_poly.support_], facecolor="none",
            edgecolor="g", marker='s', label='poly support vectors', s=50)

plt.plot(X, y_rbf, color='m', lw=lw, label='RBF model')
plt.plot(X, y_lin, color='c', lw=lw, label='Linear model')
plt.plot(X, y_poly, color='g', lw=lw, label='Polynomial model')
plt.xlabel('data')
plt.ylabel('target')
plt.title('Support Vector Regression')
plt.legend()
plt.show()
