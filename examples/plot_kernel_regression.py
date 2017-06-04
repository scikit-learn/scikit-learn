"""
========================================================================
Comparison of kernel regression (KR) and support vector regression (SVR)
========================================================================

Toy example of 1D regression using kernel regression (KR) and support vector
regression (SVR). KR provides an efficient way of selecting a kernel's
bandwidth via leave-one-out cross-validation, which is considerably faster
that an explicit grid-search as required by SVR. The main disadvantages are
that it does not support regularization and is not robust to outliers.
"""
print(__doc__)

import time

import numpy as np
from sklearn.svm import SVR
from sklearn.grid_search import GridSearchCV
from sklearn.learning_curve import learning_curve
from sklearn.kernel_regression import KernelRegression
import matplotlib.pyplot as plt

np.random.seed(0)


###############################################################################
# Generate sample data
X = np.sort(5 * np.random.rand(100, 1), axis=0)
y = np.sin(X).ravel()

###############################################################################
# Add noise to targets
y += 0.5 * (0.5 - np.random.rand(y.size))

###############################################################################
# Fit regression models
svr = GridSearchCV(SVR(kernel='rbf'), cv=5,
                   param_grid={"C": [1e-1, 1e0, 1e1, 1e2],
                               "gamma": np.logspace(-2, 2, 10)})
kr = KernelRegression(kernel="rbf", gamma=np.logspace(-2, 2, 10))
t0 = time.time()
y_svr = svr.fit(X, y).predict(X)
print "SVR complexity and bandwidth selected and model fitted in %.3f s" \
    % (time.time() - t0)
t0 = time.time()
y_kr = kr.fit(X, y).predict(X)
print "KR including bandwith fitted in %.3f s" \
    % (time.time() - t0)

###############################################################################
# Visualize models
plt.scatter(X, y, c='k', label='data')
plt.hold('on')
plt.plot(X, y_kr, c='g', label='Kernel Regression')
plt.plot(X, y_svr, c='r', label='SVR')
plt.xlabel('data')
plt.ylabel('target')
plt.title('Kernel regression versus SVR')
plt.legend()

# Visualize learning curves
plt.figure()
train_sizes, train_scores_svr, test_scores_svr = \
    learning_curve(svr, X, y, train_sizes=np.linspace(0.1, 1, 10),
                   scoring="mean_squared_error", cv=10)
train_sizes_abs, train_scores_kr, test_scores_kr = \
    learning_curve(kr, X, y, train_sizes=np.linspace(0.1, 1, 10),
                   scoring="mean_squared_error", cv=10)
plt.plot(train_sizes, test_scores_svr.mean(1), 'o-', color="r",
         label="SVR")
plt.plot(train_sizes, test_scores_kr.mean(1), 'o-', color="g",
         label="Kernel Regression")
plt.yscale("symlog", linthreshy=1e-7)
plt.ylim(-10, -0.01)
plt.xlabel("Training size")
plt.ylabel("Mean Squared Error")
plt.title('Learning curves')
plt.legend(loc="best")
plt.show()
