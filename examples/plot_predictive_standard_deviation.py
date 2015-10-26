"""
==============================================================
Comparison of predictive distributions of different regressors
==============================================================

A simple one-dimensional, noisy regression problem addressed by two different
regressors:

1. A Gaussian Process
2. Bagging with extra-trees

The regressors are fitted based on noisy observations where the magnitude of
the noise at the different training point is constant and known. Plotted are
both the mean and the pointwise 95% confidence interval of the predictions.
The mean predictions are evaluated on noise-less test data using the mean
squared error. The mean log probabilities of the noise-less test data are used
to evaluate the predictive distributions (a normal distribution with the
predicted mean and standard deviation) of the three regressors.

The mean predictions of the Gaussian Process are slightly better than those of
Random Forest. The predictive distribution (taking into account
also the predictive variance) of Bagging is however
more likely for this example.
"""
print(__doc__)

# Authors: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#          Gilles Louppe <g.louppe@gmail.com>
# Licence: BSD 3 clause

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

from sklearn.ensemble import BaggingRegressor
from sklearn.tree import ExtraTreeRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from sklearn.metrics import mean_squared_error
from sklearn.utils import check_random_state

rng = check_random_state(1)

# Observations and noise
def f(x):
    """The function to predict."""
    return x * np.sin(x)

X_train = rng.rand(50) * 10.0
X_train = np.atleast_2d(X_train).T

y_train = f(X_train).ravel()
dy = np.ones_like(y_train)
noise = rng.normal(0, dy)
y_train += noise

# Mesh the input space for evaluations of the real function, the prediction and
# its standard deviation
X_test = np.atleast_2d(np.linspace(0, 10, 1000)).T

regressors = [
    ("Gaussian Process", GaussianProcessRegressor(
        alpha=(dy / y_train) ** 2,
        kernel=1.0 * RBF() + WhiteKernel(),
        n_restarts_optimizer=10,
        random_state=rng), "b"),
    ("Bagging", BaggingRegressor(
        base_estimator=ExtraTreeRegressor(),
        n_estimators=250,
        random_state=rng), "g")
]

# Plot the function and the observations
fig = plt.figure()
plt.plot(X_test, f(X_test), 'r', label=u'$f(x) = x\,\sin(x)$')
plt.fill(np.concatenate([X_test, X_test[::-1]]),
         np.concatenate([f(X_test) - 1.9600, (f(X_test) + 1.9600)[::-1]]),
         alpha=.3, fc='r', ec='None')
plt.plot(X_train.ravel(), y_train, 'ko', zorder=5, label=u'Observations')

# Plot predictive distibutions of GP and Bagging
mse = {}
log_pdf_loss = {}

for name, regressor, color in regressors:
    regressor.fit(X_train, y_train)

    # Make the prediction on the meshed x-axis (ask for standard deviation
    # as well)
    y_pred, sigma = regressor.predict(X_test, return_std=True)

    # Compute mean-squared error and log predictive loss
    mse[name] = mean_squared_error(f(X_test), y_pred)
    log_pdf_loss[name] = \
        norm(y_pred, sigma).logpdf(f(X_test)).mean()

    # Plot 95% confidence interval based on the predictive standard deviation
    plt.plot(X_test, y_pred, color, label=name)
    plt.fill(np.concatenate([X_test, X_test[::-1]]),
             np.concatenate([y_pred - 1.9600 * sigma,
                             (y_pred + 1.9600 * sigma)[::-1]]),
             alpha=.3, fc=color, ec='None')

plt.xlabel('$x$')
plt.ylabel('$f(x)$')
plt.ylim(-10, 20)
plt.legend(loc='upper left')

print("Mean-squared error of predictors on 1000 equidistant noise-less test "
      "datapoints:\n\tBagging: %.2f\n\tGaussian Process: %.2f\n"
      % (mse["Bagging"], mse["Gaussian Process"]))

print("Mean log-probability of 1000 equidistant noise-less test datapoints\n"
      "under the (normal) predictive distribution of the predictors, i.e.,\n"
      "log N(y_true| y_pred_mean, y_pred_std) [less is better]:"
      "\n\tBagging: %.2f\n\tGaussian Process: %.2f\n"
      % (log_pdf_loss["Bagging"],
         log_pdf_loss["Gaussian Process"]))

plt.show()
