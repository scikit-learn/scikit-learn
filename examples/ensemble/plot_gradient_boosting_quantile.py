"""
=====================================================
Prediction Intervals for Gradient Boosting Regression
=====================================================

This example shows how quantile regression can be used to create prediction
intervals.
"""
# %%
# Generate some data for a synthetic regression problem by applying the f
# function to uniformly sampled random inputs.
import numpy as np
from sklearn.model_selection import train_test_split


def f(x):
    """The function to predict."""
    return x * np.sin(x)


rng = np.random.RandomState(42)
X = np.atleast_2d(rng.uniform(0, 10.0, size=1000)).T
X = X.astype(np.float32)
y = f(X).ravel()

# %%
# To make the problem interesting, add centered `log-normal distributed
# <https://en.wikipedia.org/wiki/Log-normal_distribution>`_ random noise to the
# target variable.
#
# The lognormal distribution is very skewed, meaning that it is likely to get
# large outliers but impossible to observe small outliers.
sigma = 1.2
noise = rng.lognormal(sigma=sigma, size=y.shape) - np.exp(sigma ** 2 / 2)
y += noise

# %%
# Split into train, test datasets:
X_train, X_test, y_train, y_test = train_test_split(X, y)

# %%
# Fit gradient boosting models trained with the quantile loss and
# alpha=0.05, 0.5, 0.95.
#
# The models obtained for alpha=0.05 and alpha=0.95 produce a 90% confidence
# interval (95% - 5% = 90%).
#
# The model trained with alpha=0.5 produces a regression of the median: on
# average, there should be the same number of target observations above and
# below the predicted values.
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import pinball_loss, mean_squared_error


all_models = {}
common_params = dict(
    learning_rate=0.05,
    n_estimators=250,
    max_depth=2,
    min_samples_leaf=9,
    min_samples_split=9,
)
for alpha in [0.05, 0.5, 0.95]:
    gbr = GradientBoostingRegressor(loss='quantile', alpha=alpha,
                                    **common_params)
    all_models["q %1.2f" % alpha] = gbr.fit(X_train, y_train)

# %%
# For the sake of comparison, also fit a baseline model trained with the usual
# least squares loss (ls), also known as the mean squared error (MSE).
gbr_ls = GradientBoostingRegressor(loss='ls', **common_params)
all_models["ls"] = gbr_ls.fit(X_train, y_train)

# %%
# Create an evenly spaced set of input values spanning the [0, 10] range.
xx = np.atleast_2d(np.linspace(0, 10, 1000)).T
xx = xx.astype(np.float32)

# %%
# Plot the true function (expected mean), the prediction of the conditional
# mean (least squares loss), the conditional median and the conditional 90%
# interval (from 5th to 95th conditional percentiles).
import matplotlib.pyplot as plt


y_pred = all_models['ls'].predict(xx)
y_lower = all_models['q 0.05'].predict(xx)
y_upper = all_models['q 0.95'].predict(xx)
y_med = all_models['q 0.50'].predict(xx)

fig = plt.figure()
plt.plot(xx, f(xx), 'g:', label=r'$f(x) = x\,\sin(x)$')
plt.plot(X_test, y_test, 'b.', markersize=10, label='Test observations')
plt.plot(xx, y_med, 'r-', label='Predicted median', color="orange")
plt.plot(xx, y_pred, 'r-', label='Predicted mean')
plt.plot(xx, y_upper, 'k-')
plt.plot(xx, y_lower, 'k-')
plt.fill(np.concatenate([xx, xx[::-1]]),
         np.concatenate([y_upper, y_lower[::-1]]),
         alpha=.5, fc='b', ec='None', label='Predicted 90% interval')
plt.xlabel('$x$')
plt.ylabel('$f(x)$')
plt.ylim(-10, 20)
plt.legend(loc='upper left')
plt.show()

# %%
# Note that the predicted median is on average below the predicted mean as the
# noise is skewed towards high values (large outliers). Also note that the
# median estimate is smoother because of its natural robustness to outliers.
#
# Also observe that the inductive bias of gradient boosting trees is
# unfortunately preventing our 0.05 quantile to fully capture the sinoisoidal
# shape of the signal, in particular around x=8. Tuning hyper-parameters could
# potentially reduce this effect a bit.
#
# Measure the models with :func:`mean_square_error` and :func:`pinball_loss`
# metrics on the training dataset.
from pandas import DataFrame

results = []
for name, gbr in sorted(all_models.items()):
    metrics = {'model': name}
    y_pred = gbr.predict(X_train)
    for alpha in [0.05, 0.5, 0.95]:
        metrics["pbl=%1.2f" % alpha] = pinball_loss(
            y_train, y_pred, alpha=alpha)
    metrics['MSE'] = mean_squared_error(y_train, y_pred)
    results.append(metrics)
DataFrame(results).set_index('model')

# %%
# One column shows all models evaluated by the same metric. The minimum number
# on a column should be obtained when the model is trained and measured with
# the same metric. This should be always the case on the training set if the
# training converged.
#
# Note that because the target noise is skewed by the presence of large
# outliers, the expected conditional mean and conditional median are
# signficiantly different and therefore one could not use the least squares
# model get a good estimation of the conditional median nor the converse.
#
# If the target distribution had not been skewed and had no ouliers (e.g. with
# a Gaussian noise), then median estimator and the least squares estimator
# would have yielded similar predictions.
#
# We then do the same on the test set.
results = []
for name, gbr in sorted(all_models.items()):
    metrics = {'model': name}
    y_pred = gbr.predict(X_test)
    for alpha in [0.05, 0.5, 0.95]:
        metrics["pbl=%1.2f" % alpha] = pinball_loss(
            y_test, y_pred, alpha=alpha)
    metrics['MSE'] = mean_squared_error(y_test, y_pred)
    results.append(metrics)
DataFrame(results).set_index('model')

# %%
# Errors are higher meaning the models slightly overfitted the data. It still
# shows the minimum of a metric is obtained when the model is trained by
# minimizing this same metric.
