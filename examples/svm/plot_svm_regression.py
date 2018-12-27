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
svr_poly = SVR(kernel='poly', C=100, degree=3, epsilon=.1, coef0=1)
y_rbf = svr_rbf.fit(X, y).predict(X)
y_lin = svr_lin.fit(X, y).predict(X)
y_poly = svr_poly.fit(X, y).predict(X)

# #############################################################################
# Look at the results
lw = 2
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 10), sharey=True)

ax1.plot(X, y_rbf, color='m', lw=lw, label='RBF model')
ax1.scatter(X[svr_rbf.support_], y[svr_rbf.support_], facecolor="none",
            edgecolor="m", label='rbf support vectors', s=50)
ax1.scatter(X[np.setdiff1d(np.arange(len(X)), svr_rbf.support_)],
            y[np.setdiff1d(np.arange(len(X)), svr_rbf.support_)],
            facecolor="none",
            edgecolor="k", label='other training data', s=50)
ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
           ncol=1, fancybox=True, shadow=True)


ax2.plot(X, y_lin, color='c', lw=lw, label='Linear model')
ax2.scatter(X[svr_lin.support_], y[svr_lin.support_], facecolor="none",
            edgecolor="c", label='linear support vectors', s=50)
ax2.scatter(X[np.setdiff1d(np.arange(len(X)), svr_lin.support_)],
            y[np.setdiff1d(np.arange(len(X)), svr_lin.support_)],
            facecolor="none",
            edgecolor="k", label='other training data', s=50)
ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
           ncol=1, fancybox=True, shadow=True)


ax3.plot(X, y_poly, color='g', lw=lw, label='Polynomial model')
ax3.scatter(X[svr_poly.support_], y[svr_poly.support_], facecolor="none",
            edgecolor="g", label='poly support vectors', s=50)
ax3.scatter(X[np.setdiff1d(np.arange(len(X)), svr_poly.support_)],
            y[np.setdiff1d(np.arange(len(X)), svr_poly.support_)],
            facecolor="none",
            edgecolor="k", label='other training data', s=50)
ax3.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1),
           ncol=1, fancybox=True, shadow=True)

fig.text(0.5, 0.04, 'data', ha='center', va='center')
fig.text(0.06, 0.5, 'target', ha='center', va='center', rotation='vertical')
fig.suptitle("Support Vector Regression", fontsize=14)
plt.show()