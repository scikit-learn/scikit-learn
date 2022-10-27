"""
=====================================================
Probabilistic Gradient Boosting Regression
=====================================================

This example shows how probalistic gradient boosting regression can be
used to create a rich distribution of regression outputs with proper
confidence intervals.

We will compare our probabilistic model to quantile
regression models, demonstrating that we can achieve comparable performance
training only a single model compared to training models for multiple
quantiles.

The key benefit of this approach is that a rich distribution of outputs
can be obtained with only a single model.

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
# Fitting least squares regressors and quantile regressors
# --------------------------------------------------------
#
# We fit a baseline model trained with the (mean) squared error (MSE)
from sklearn.ensemble import HistGradientBoostingRegressor

all_models = {}
common_params = dict(
    learning_rate=0.05,
    max_iter=200,
    max_depth=2,
    min_samples_leaf=9,
)
gbr_ls = HistGradientBoostingRegressor(loss="squared_error", **common_params)
all_models["mse"] = gbr_ls.fit(X_train, y_train)
# %%
# For the sake of comparison, we train multiple models with the
# the quantile loss and quantiles=[0.05, 0.10, ..., 0.90, 0.95]
quantiles = np.linspace(0.05, 0.95, 19)
for quantile in quantiles:
    gbr = HistGradientBoostingRegressor(
        loss="quantile", quantile=quantile, **common_params
    )
    all_models["q %1.2f" % quantile] = gbr.fit(X_train, y_train)

# %%
# Create an evenly spaced evaluation set of input values spanning the [0, 10]
# range.
xx = np.atleast_2d(np.linspace(0, 10, 1000)).T

# %%
# We obtain the empirical mean and standard deviation from our single least
# squares model and use these statistics to sample from a normal distribution.
# Based on many estimates, we can obtain the empirical estimated quantiles.
#
# We plot the true conditional mean function f, the predictions of the conditional
# mean (loss equals squared error), and the conditional 90% interval
# (from 5th to 95th conditional percentiles).
import matplotlib.pyplot as plt

# First, obtain the outputs of the probabilistic gradient boosting model
# We both the mean and standard deviation by setting `return_std=True` in the predict
# method of our model.
y_mean, y_std = all_models["mse"].predict(xx, return_std=True)
# Using the mean and standard deviation, we can sample from a normal distribution
y_dist = all_models["mse"].sample(
    y_mean, y_std, n_estimates=1_000, distribution="normal"
)
# Based on our 1_000 estimates, we can obtain quantile estimates for every sample.
y_lower = np.quantile(y_dist, q=0.05, axis=0)
y_upper = np.quantile(y_dist, q=0.95, axis=0)
# Plot the outputs of the probabilistic gradient boosting model
fig, ax = plt.subplots(1, 2, figsize=(20, 10))
ax[0].plot(xx, f(xx), "g:", linewidth=3, label=r"$f(x) = x\,\sin(x)$")
ax[0].plot(X_test, y_test, "b.", markersize=10, label="Test observations")
ax[0].plot(xx, y_mean, "r-", label="Predicted mean")
ax[0].plot(xx, y_upper, "k-")
ax[0].plot(xx, y_lower, "k-")
ax[0].fill_between(
    xx.ravel(), y_lower, y_upper, alpha=0.4, label="Predicted 90% interval"
)
ax[0].set_xlabel("$x$")
ax[0].set_ylabel("$f(x)$")
ax[0].set_ylim(-10, 25)
ax[0].legend(loc="upper left")
ax[0].set_title("Probabilistic Gradient Boosting")

# Second, obtain the 90% interval from the quantile regression models
y_lower = all_models["q 0.05"].predict(xx)
y_upper = all_models["q 0.95"].predict(xx)

ax[1].plot(xx, f(xx), "g:", linewidth=3, label=r"$f(x) = x\,\sin(x)$")
ax[1].plot(X_test, y_test, "b.", markersize=10, label="Test observations")
ax[1].plot(xx, y_mean, "r-", label="Predicted mean")
ax[1].plot(xx, y_upper, "k-")
ax[1].plot(xx, y_lower, "k-")
ax[1].fill_between(
    xx.ravel(), y_lower, y_upper, alpha=0.4, label="Predicted 90% interval"
)
ax[1].set_xlabel("$x$")
ax[1].set_ylabel("$f(x)$")
ax[1].set_ylim(-10, 25)
ax[1].legend(loc="upper left")
ax[1].set_title("Quantile Gradient Boosting")
plt.show()
# %%
# Comparing the probabilistic single model with the separately trained quantile
# regression models we see that the first plot shows a symmetric shape
# of the confidence interval around the mean, which is caused by the fact that
# we sample from a Normal distribution, which is symmetric. The quantile
# models seem to track the assymmetry a bit better, but the overall confidence
# interval of the quantile models also seems to be a bit too wide. Let's evaluate
# the confidence interval at several sizes.
#
# Calibration of the confidence interval
# --------------------------------------
#
# We can evaluate the ability of the probabilistic model at
# producing a well-calibrated confidence interval across various sizes of the
# confidence interval.
#
# To do this we can compute the fraction of observations that fall between the
# predictions:


def coverage_fraction(y, y_low, y_high):
    return np.mean(np.logical_and(y >= y_low, y < y_high))


# We compute the coverage fraction across the 90%, 80%, ..., 20%, 10% confidence
# interval for the probabilistic model and for all the quantile models.
y_mean, y_std = all_models["mse"].predict(X_train, return_std=True)
y_dist_train = all_models["mse"].sample(
    y_mean, y_std, n_estimates=1_000, distribution="normal"
)
lower_quantiles = quantiles[: len(quantiles) // 2]
upper_quantiles = quantiles[len(quantiles) // 2 + 1 :][::-1]
n_quantile_pairs = len(lower_quantiles)

for i in range(n_quantile_pairs):
    lower_quantile = lower_quantiles[i]
    upper_quantile = upper_quantiles[i]
    interval_size = 100 * (upper_quantile - lower_quantile)
    cov_fraction_prob_model = coverage_fraction(
        y_train,
        y_low=np.quantile(y_dist_train, q=lower_quantile, axis=0),
        y_high=np.quantile(y_dist_train, q=upper_quantile, axis=0),
    )
    cov_fraction_quantile_models = coverage_fraction(
        y_train,
        all_models["q %1.2f" % lower_quantile].predict(X_train),
        all_models["q %1.2f" % upper_quantile].predict(X_train),
    )
    print(
        f"Coverage @{interval_size:0.0f}% confidence interval: \n             "
        f" Probabilistic model: {cov_fraction_prob_model*100:0.0f}% \n             "
        f" Quantile models    : {cov_fraction_quantile_models*100:0.0f}%"
    )
# %%
# As we can see, our probabilistic model is quite well calibrated on the training set
# and performs only slightly worse than the quantile models, even though the quantile
# models have been specifically trained on every quantile.

y_mean, y_std = all_models["mse"].predict(X_test, return_std=True)
y_dist_test = all_models["mse"].sample(
    y_mean, y_std, n_estimates=1_000, distribution="normal"
)

for i in range(n_quantile_pairs):
    lower_quantile = lower_quantiles[i]
    upper_quantile = upper_quantiles[i]
    interval_size = 100 * (upper_quantile - lower_quantile)
    cov_fraction_quantile_models = coverage_fraction(
        y_test,
        all_models["q %1.2f" % lower_quantile].predict(X_test),
        all_models["q %1.2f" % upper_quantile].predict(X_test),
    )
    cov_fraction_prob_model = coverage_fraction(
        y_test,
        y_low=np.quantile(y_dist_test, q=lower_quantile, axis=0),
        y_high=np.quantile(y_dist_test, q=upper_quantile, axis=0),
    )
    print(
        f"Coverage @{interval_size:0.0f}% confidence interval: \n             "
        f" Probabilistic model: {cov_fraction_prob_model*100:0.0f}% \n             "
        f" Quantile models    : {cov_fraction_quantile_models*100:0.0f}%"
    )
# %%
# On the test set, our single probabilistic model is better or equal across
# most confidence intervals. This is great: it means that on this task, we could
# just revert to having a single probabilistic model and sample from that
# model, rather than requiring the need to train separate models for every quantile.
#
# Trying different distributions for our probabilistic model
# ----------------------------------------------------------
#
# A key benefit of the probabilistic approach is that we can try out different
# distributions and find out which distribution fits our data the best
#
# We already saw that the Normal distribution gave quite good results on the training
# set, but perhaps if we use a Student's-t distribution we can capture the fat tail
# of the target at the top a bit better. We use a Student's-t distribution with
# 3-degrees of freedom: this results in a distribution similar to the Normal
# distribution, but with fatter tails.
y_mean, y_std = all_models["mse"].predict(X_test, return_std=True)
y_dist_test_normal = all_models["mse"].sample(
    y_mean, y_std, n_estimates=1_000, distribution="normal"
)
y_dist_test_studentt = all_models["mse"].sample(
    y_mean,
    y_std,
    n_estimates=1_000,
    distribution="studentt",
    studentt_degrees_of_freedom=3,
)

lower_quantiles = quantiles[: len(quantiles) // 2]
upper_quantiles = quantiles[len(quantiles) // 2 + 1 :][::-1]
n_quantile_pairs = len(lower_quantiles)

for i in range(n_quantile_pairs):
    lower_quantile = lower_quantiles[i]
    upper_quantile = upper_quantiles[i]
    interval_size = 100 * (upper_quantile - lower_quantile)
    cov_fraction_normal = coverage_fraction(
        y_test,
        y_low=np.quantile(y_dist_test_normal, q=lower_quantile, axis=0),
        y_high=np.quantile(y_dist_test_normal, q=upper_quantile, axis=0),
    )
    cov_fraction_studentt = coverage_fraction(
        y_test,
        y_low=np.quantile(y_dist_test_studentt, q=lower_quantile, axis=0),
        y_high=np.quantile(y_dist_test_studentt, q=upper_quantile, axis=0),
    )
    print(
        f"Coverage @{interval_size:0.0f}% confidence interval: \n              Normal"
        f" distribution        : {cov_fraction_normal*100:0.0f}% \n             "
        f" Student's-t(3) distribution: {cov_fraction_studentt*100:0.0f}%"
    )
# %%
# The calibration of the Student's-t(3) distribution is unfortunately not better
# than the Normal distribution: across every confidence interval it performs
# worse than the Normal distribution.
