"""
===================
Quantile regression
===================

This example illustrates how quantile regression can predict non-trivial
conditional quantiles.

The left figure shows the case when the error distribution is normal,
but has non-constant variance, i.e. with heteroscedasticity.

The right figure shows an example of an asymmetric error distribution,
namely the Pareto distribution.
"""

print(__doc__)

# Authors: David Dale <dale.david@mail.ru>
#          Christian Lorentzen <lorentzen.ch@gmail.com>
#          Guillaume Lemaitre <glemaitre58@gmail.com>
# License: BSD 3 clause

# %%
# Dataset generation
# ------------------
#
# To illustrate the behaviour of quantile regression, we will generate two
# synthetic datasets. The true generative random processess for both datasets will be
# composed by an expected value with a linear relationship with a single feature `x`.
import numpy as np

rng = np.random.RandomState(42)
x = np.linspace(start=0, stop=10, num=100)
X = x[:, np.newaxis]
y_true_mean = 10 + 0.5 * x

# %%
# We will create two subsequent problems by changing the distribution of the target `y`
# while keeping the same expected value:
#
# - in the first case, a heteroscedastic Normal noise is added;
# - in the second case, an asymmetric Pareto noise is added.
y_normal = y_true + rng.normal
    loc=0, scale=0.5 + 0.5 * x, size=x.shape[0]
)
a = 5
y_pareto = y_true + 10 * (rng.pareto(a, size=x.shape[0]) - 1 / (a - 1))

# %%
# Let's first visualize the datasets as well as the distribution of the residuals `y - mean(y)`.
import matplotlib.pyplot as plt

_, axs = plt.subplots(
    nrows=2, ncols=2, figsize=(15, 10), sharex="col", sharey="col"
)

axs[0, 0].plot(x, y_true, label="True generative process")
axs[0, 0].scatter(
    x, y_noisy_normal, color="black", alpha=0.5, label="Noisy observations"
)
axs[0, 1].hist(y_true - y_noisy_normal, edgecolor="black")


axs[1, 0].plot(x, y_true, label="True generative process")
axs[1, 0].scatter(
    x, y_noisy_pareto, color="black", alpha=0.5, label="Noisy observations"
)
axs[1, 1].hist(y_true - y_noisy_pareto, edgecolor="black")

axs[0, 0].set_title("Dataset with heteroscedastic Normal distributed targets")
axs[1, 0].set_title("Dataset with asymmetric Pareto distributed target")
axs[0, 1].set_title(
    "Residuals distribution for heteroscedastic Normal distributed targets"
)
axs[1, 1].set_title(
    "Residuals distribution for asymmetric Pareto distributed target"
)
axs[0, 0].legend()
axs[1, 0].legend()
axs[0, 0].set_ylabel("y")
axs[1, 0].set_ylabel("y")
axs[1, 0].set_xlabel("x")
_ = axs[1, 1].set_xlabel("Residuals")

# %%
# With the heteroscedastic Normal distributed target, we observe that the
# variance of the noise is increasing when the value of the feature `x` is
# increasing.
#
# With the asymmetric Pareto distributed target, we observe that the positive
# residuals are bounded.
#
# These types of noisy targets violate one of the key assumption made when one
# wants to use a :class:`~sklearn.linear_model.LinearRegression`: the noise is
# expected to be normally distributed with a constant variance.
#
# Thus, :class:`~sklearn.linear_model.QuantileRegressor` can be used in this
# case to estimate the conditional quantile (e.g. median) instead of the
# conditional mean as in :class:`sklearn.linear_model.QuantileRegressor`.
#
# In the remainder of this tutorial, we will show how
# :class:`~sklearn.linear_model.QuantileRegressor` can be used in practice and
# give the intuition into the properties of the fitted models. Finally,
# we will compare the both :class:`~sklearn.linear_model.QuantileRegressor`
# and :class:`~sklearn.linear_model.LinearRegression`.
#
# Fitting a `QuantileRegressor`
# -----------------------------
#
# In this section, we want to estimate the conditional median as well as
# a low and high quantiles fixed at 5% and 95%, respectively. Thus, we will get
# three linear models, one for each quantile.
#
# We will use the quantile at 5% and 95% to find the training samples that
# are outside of this interval.
from sklearn.linear_model import QuantileRegressor

quantiles = [0.05, 0.5, 0.95]
predictions = {}
out_bounds_predictions = np.zeros_like(y_true, dtype=np.bool_)
for quantile in quantiles:
    qr = QuantileRegressor(quantile=quantile, alpha=0)
    y_pred = qr.fit(X, y_noisy_normal).predict(X)
    predictions[quantile] = y_pred

    if quantile == min(quantiles):
        out_bounds_predictions = np.logical_or(
            out_bounds_predictions, y_pred >= y_noisy_normal
        )
    elif quantile == max(quantiles):
        out_bounds_predictions = np.logical_or(
            out_bounds_predictions, y_pred <= y_noisy_normal
        )

# %%
# Now, we can plot the three linear models and the distinguished samples that
# are within the 5%-95% quantiles from samples that are outside these bounds.
plt.plot(X, y_true, color="black", linestyle="dashed", label="True mean")

for quantile, y_pred in predictions.items():
    plt.plot(X, y_pred, label=f"Quantile: {quantile}")

plt.scatter(
    x[out_bounds_predictions],
    y_noisy_normal[out_bounds_predictions],
    color="black",
    marker="+",
    alpha=0.5,
    label="Outside interval",
)
plt.scatter(
    x[~out_bounds_predictions],
    y_noisy_normal[~out_bounds_predictions],
    color="black",
    alpha=0.5,
    label="Inside interval",
)

plt.legend()
plt.xlabel("x")
plt.ylabel("y")
_ = plt.title("Quantiles of heteroscedastic Normal distributed target")

# %%
# Since the noise is still Normally distributed, in particular is symmetric, the true conditional mean
# and the true conditional median coincide. Indeed, we see that the estimated median almost
# hits the true mean. We observe the effect of having an
# increasing noise variance on the 5% and 95% quantiles: the slopes of those quantiles is
# very different and the interval between them becomes wider with increasing `x`.
#
# We can repeat the same experiment using the asymmetric Pareto distributed
# target.
quantiles = [0.05, 0.5, 0.95]
predictions = {}
out_bounds_predictions = np.zeros_like(y_true, dtype=np.bool_)
for quantile in quantiles:
    qr = QuantileRegressor(quantile=quantile, alpha=0)
    y_pred = qr.fit(X, y_noisy_pareto).predict(X)
    predictions[quantile] = y_pred

    if quantile == min(quantiles):
        out_bounds_predictions = np.logical_or(
            out_bounds_predictions, y_pred >= y_noisy_pareto
        )
    elif quantile == max(quantiles):
        out_bounds_predictions = np.logical_or(
            out_bounds_predictions, y_pred <= y_noisy_pareto
        )

# %%
plt.plot(X, y_true, color="black", linestyle="dashed", label="True mean")

for quantile, y_pred in predictions.items():
    plt.plot(X, y_pred, label=f"Quantile: {quantile}")

plt.scatter(
    x[out_bounds_predictions],
    y_noisy_pareto[out_bounds_predictions],
    color="black",
    marker="+",
    alpha=0.5,
    label="Outside interval",
)
plt.scatter(
    x[~out_bounds_predictions],
    y_noisy_pareto[~out_bounds_predictions],
    color="black",
    alpha=0.5,
    label="Inside interval",
)

plt.legend()
plt.xlabel("x")
plt.ylabel("y")
_ = plt.title("Quantiles of asymmetric Pareto distributed target")


# %%
# Due to the asymmetry distribution of the noise, we observe that the true mean
# and estimated conditional median are different. We also observe that each
# quantile model has different parameters to better fit the desired
# quantile.
#
# Comparing `QuantileRegressor` and `LinearRegression`
# ----------------------------------------------------
#
# In this section, we will linger on the difference regarding the error that
# :class:`~sklearn.linear_model.QuantileRegressor` and
# :class:`~sklearn.linear_model.LinearRegression` are minimizing.
#
# Indeed, :class:`~sklearn.linear_model.LinearRegression` is a least squares
# approach minimizing the mean squared error between the training and predicted
# targets. In contrast, :class:`~sklearn.linear_model.QuantileRegressor`
# with `quantile=0.5` minimizes the mean absolute error instead.
#
# The second part of the example shows that
# :class:`~sklearn.linear_model.LinearRegression` minimizes the mean squared
# error (MSE) in order to predict the mean, while
# :class:`~sklearn.linear_model.QuantileRegressor` with `quantile=0.5`
# minimizes the mean absolute error (MAE) in order to predict the median.
#
# Let's first compute the training errors of such models in terms of mean
# squared error and mean absolute error. We will use the asymmetric Pareto
# distributed target to make it more interesting as mean and median are not equal.
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error

linear_regression = LinearRegression()
quantile_regression = QuantileRegressor(quantile=0.5, alpha=0)

y_pred_lr = linear_regression.fit(X, y_noisy_pareto).predict(X)
y_pred_qr = quantile_regression.fit(X, y_noisy_pareto).predict(X)

print(
    f"""Training error (in-sample performance)
    {linear_regression.__class__.__name__}:
    MAE = {mean_absolute_error(y_noisy_pareto, y_pred_lr):.3f}
    MSE = {mean_squared_error(y_noisy_pareto, y_pred_lr):.3f}
    {quantile_regression.__class__.__name__}:
    MAE = {mean_absolute_error(y_noisy_pareto, y_pred_qr):.3f}
    MSE = {mean_squared_error(y_noisy_pareto, y_pred_qr):.3f}
    """
)

# %%
# On the training set, we see that MAE is lower for
# :class:`~sklearn.linear_model.QuantileRegressor` than
# :class:`~sklearn.linear_model.LinearRegression`. In contrast to that, MSE is lower
# for :class:`~sklearn.linear_model.LinearRegressor` than
# :class:`~sklearn.linear_model.QuantileRegressor`. These results confirms that
# MAE is the loss minimized by
# :class:`~sklearn.linear_model.QuantileRegressor` while MSE is the loss
# minimized :class:`~sklearn.linear_model.LinearRegression`.
#
# We can make a similar evaluation but looking a the test error obtained by
# cross-validation.
from sklearn.model_selection import cross_validate

cv_results_lr = cross_validate(
    linear_regression,
    X,
    y_noisy_pareto,
    cv=3,
    scoring=["neg_mean_absolute_error", "neg_mean_squared_error"],
)
cv_results_qr = cross_validate(
    quantile_regression,
    X,
    y_noisy_pareto,
    cv=3,
    scoring=["neg_mean_absolute_error", "neg_mean_squared_error"],
)
print(
    f"""Test error (cross-validated performance)
    {linear_regression.__class__.__name__}:
    MAE = {-cv_results_lr["test_neg_mean_absolute_error"].mean():.3f}
    MSE = {-cv_results_lr["test_neg_mean_squared_error"].mean():.3f}
    {quantile_regression.__class__.__name__}:
    MAE = {-cv_results_qr["test_neg_mean_absolute_error"].mean():.3f}
    MSE = {-cv_results_qr["test_neg_mean_squared_error"].mean():.3f}
    """
)

# %%
# We reach similar conclusions on out-of-samples evaluation.
