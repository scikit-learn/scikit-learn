"""
=================================
Poisson Regression with B-Splines
=================================

As in the :ref:`sphx_glr_auto_examples_ensemble_plot_adaboost_regression.py`
example, a Poisson regression with penalized B-splines (P-splines) [1]_ is
fitted on slightly different sinusoidal, Poisson distributed data and
compared to an AdaBoost model with decision trees.
One can see, that this is a hard problem for both estimators.

.. [1] Eilers, Paul H. C.; Marx, Brian D. "Flexible smoothing with B -splines
       and penalties". Statist. Sci. 11 (1996), no. 2, 89--121.
       `doi:10.1214/ss/1038425655
       <https://projecteuclid.org/euclid.ss/1038425655>`_

"""
print(__doc__)

# Author: Christian Lorentzen <lorentzen.ch@gmail.com>
# based on the AdaBoost regression example from Noel Dawe <noel.dawe@gmail.com>
# License: BSD 3 clause

# importing necessary libraries
import numpy as np
from scipy.linalg import toeplitz
# from scipy.interpolate import BSpline
from scipy.interpolate import splev
import matplotlib.pyplot as plt
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.linear_model import GeneralizedLinearRegressor


# Create the dataset
xmin, xmax = 0, 6
rng = np.random.RandomState(1)
X = np.linspace(xmin, xmax, 500)[:, np.newaxis]
y_true = 0.5 * (2.1 + np.sin(X).ravel() + np.sin(6 * X).ravel())
y = rng.poisson(y_true, X.shape[0])

# b-spline basis
nknots, degree = 40, 3
ns = nknots - degree - 1  # number of base spline functions
dx = (xmax - xmin) / (nknots - 1 - 2 * degree)
knots = np.linspace(xmin - degree * dx, 6 + degree * dx, nknots)
coef = np.zeros(ns)
splineBasis = np.empty((X.shape[0], ns), dtype=float)
for i in range(ns):
    coef[i] = 1
#    splineBasis[:, i] = BSpline(knots, coef, degree, extrapolate=False)(X) \
#        .ravel()
    splineBasis[:, i] = splev(X, (knots, coef, degree)).ravel()
    coef[i] = 0

# second order difference matrix
P2 = toeplitz([2, -1] + [0] * (ns - 2)).astype(float)
P2[0, 0] = P2[-1, -1] = 1

# Fit regression model
regr_1 = AdaBoostRegressor(DecisionTreeRegressor(max_depth=4),
                           n_estimators=10, random_state=rng)

regr_2 = GeneralizedLinearRegressor(family='poisson', link='log',
                                    fit_intercept=True, alpha=0.02,
                                    l1_ratio=0.1, P2=P2)

regr_1.fit(X, y)
regr_2.fit(splineBasis, y)

# Predict
y_1 = regr_1.predict(X)
y_2 = regr_2.predict(splineBasis)

# Plot the results
plt.figure()
plt.plot(X, y_true, c="b", label="true mean")
plt.scatter(X, y, c="k", marker='.', label="training samples")
plt.plot(X, y_1, c="g", label="AdaBoost n_estimator=10", linewidth=2)
plt.plot(X, y_2, c="r", label="Poisson GLM with B-splines", linewidth=2)
plt.xlabel("data")
plt.ylabel("target")
plt.title("Regression Comparison")
plt.legend()
plt.show()
