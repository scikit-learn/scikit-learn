"""
=====================================================
Prediction Intervals for Support Vector Regression
=====================================================

This example shows how quantile regression can be used to create prediction
intervals for non-linear functions.

"""

# %%
# Generate some data for a synthetic regression problem by applying the
# function f to uniformly sampled random inputs.
import numpy as np
from sklearn.model_selection import train_test_split


def f(x):
    """The function to predict."""
    return x * np.sin(x)


rng = np.random.RandomState(42)
X = np.atleast_2d(rng.uniform(0, 10.0, size=1000)).T
expected_y = f(X).ravel()

# %%
# To make the problem interesting, we generate observations of the target y as
# the sum of a deterministic term computed by the function f and a random noise
# term that follows a centered `log-normal
# <https://en.wikipedia.org/wiki/Log-normal_distribution>`_. To make this even
# more interesting we consider the case where the amplitude of the noise
# depends on the input variable x (heteroscedastic noise).
#
# The lognormal distribution is non-symmetric and long tailed: observing large
# outliers is likely but it is impossible to observe small outliers.
sigma = 0.5 + X.ravel() / 10
noise = rng.lognormal(sigma=sigma) - np.exp(sigma**2 / 2)
y = expected_y + noise

# %%
# Split into train, test datasets:
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

# %%
# Fitting non-linear quantile regressors
# --------------------------------------------------------
#
# Fit Quantile SVR with quantiles
# quantile=0.05, 0.5, 0.95.
#
# The models obtained for quantile=0.05 and alpha=0.95 produce a 90% confidence
# interval (95% - 5% = 90%).
#
# The model trained with alpha=0.5 produces a regression of the median: on
# average, there should be the same number of target observations above and
# below the predicted values.
from sklearn.svm import QuantileSVR
from sklearn.kernel_ridge import KernelRidge
from sklearn.metrics import mean_pinball_loss, mean_squared_error


all_models = {}
for alpha in [0.05, 0.5, 0.95]:
    qsvr = QuantileSVR(quantile=alpha, C=100, gamma=0.2)
    all_models["q %1.2f" % alpha] = qsvr.fit(X_train, y_train)

# %%
# For the sake of comparison, we also fit a kernel ridge regression model that
# optimizes the mean square error (MSE) instead of the pinball loss. The
# parameter `alpha` corresponds to `1 / (2C)` in support vector regression.
kr = KernelRidge(alpha=0.005, kernel="rbf", gamma=0.2)
all_models["mse"] = kr.fit(X_train, y_train)

# %%
# Create an evenly spaced evaluation set of input values spanning the [0, 10]
# range.
xx = np.atleast_2d(np.linspace(0, 10, 1000)).T

# %%
# Plot the true conditional mean function f, the conditional median and the conditional
# 90% interval (from 5th to 95th conditional percentiles).
import matplotlib.pyplot as plt

y_pred = all_models["mse"].predict(xx)
y_lower = all_models["q 0.05"].predict(xx)
y_upper = all_models["q 0.95"].predict(xx)
y_med = all_models["q 0.50"].predict(xx)

fig = plt.figure()
plt.plot(xx, f(xx), color="k", linestyle="--", label=r"$f(x) = x\,\sin(x)$")
plt.plot(xx, y_pred, color="tab:orange", label="Predicted mean")
plt.plot(xx, y_med, color="tab:blue", label="Predicted median")
plt.fill_between(
    xx.ravel(),
    y_lower,
    y_upper,
    alpha=0.35,
    color="tab:blue",
    label="Predicted 90% interval",
)
plt.scatter(X_test, y_test, color="k", s=5, label="Test observations", alpha=0.5)
plt.xlabel("$x$")
plt.ylabel("$f(x)$")
plt.ylim(-10, 25)
plt.legend(loc="upper left")
plt.show()


# %%
# Comparing the predicted median with the predicted mean, we note that the
# median is on average below the mean as the noise is skewed towards high
# values (large outliers).
#
# Analysis of the error metrics
# -----------------------------
#
# Measure the models with :func:`mean_squared_error` and
# :func:`mean_pinball_loss` metrics on the test dataset.
import pandas as pd


def highlight_min(x):
    x_min = x.min()
    return ["font-weight: bold" if v == x_min else "" for v in x]


results = []
for name, mod in sorted(all_models.items()):
    metrics = {"model": name}
    y_pred = mod.predict(X_test)
    for alpha in [0.05, 0.5, 0.95]:
        metrics["pbl=%1.2f" % alpha] = mean_pinball_loss(y_test, y_pred, alpha=alpha)
    metrics["MSE"] = mean_squared_error(y_test, y_pred)
    results.append(metrics)

pd.DataFrame(results).set_index("model").style.apply(highlight_min)

# %%
# One column shows all models evaluated by the same metric. The minimum number
# on a column should be obtained when the model is trained and measured with
# the same metric. This should be always the case on the training set if the
# training converged.
#
# Calibration of the confidence interval
# --------------------------------------
#
# We can also evaluate the ability of the two extreme quantile estimators at
# producing a well-calibrated conditational 90%-confidence interval.
#
# To do this we can compute the fraction of observations that fall between the
# predictions:
def coverage_fraction(y, y_low, y_high):
    return np.mean(np.logical_and(y >= y_low, y <= y_high))


coverage_fraction(
    y_train,
    all_models["q 0.05"].predict(X_train),
    all_models["q 0.95"].predict(X_train),
)

# %%
# On the training set the calibration is very close to the expected coverage
# value for a 90% confidence interval.
coverage_fraction(
    y_test, all_models["q 0.05"].predict(X_test), all_models["q 0.95"].predict(X_test)
)

# %%
# On the test set, the estimated confidence interval is slightly too wide.
# Note, however, that we would need to wrap those metrics in a cross-validation
# loop to assess their variability under data resampling.
