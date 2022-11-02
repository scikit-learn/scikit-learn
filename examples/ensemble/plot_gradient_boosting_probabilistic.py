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
# We need to set the flag `with_variance=True` in order to be able
# to get probabilistic outputs from our model.
from sklearn.ensemble import HistGradientBoostingRegressor

all_models = {}
common_params = dict(
    learning_rate=0.05,
    max_iter=200,
    max_depth=2,
    min_samples_leaf=9,
)
gbr_ls = HistGradientBoostingRegressor(
    loss="squared_error", with_variance=True, **common_params
)
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
y_dist = all_models["mse"].sample(y_mean, y_std, n_estimates=10_000)
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
# producing well-calibrated confidence intervals across various sizes of the
# confidence interval.
#
# To do this we can compute the fraction of observations that fall between the
# predictions:


def coverage_fraction(y, y_low, y_high):
    return np.mean(np.logical_and(y >= y_low, y < y_high))


# We compute the coverage fraction across the 90%, 80%, ..., 20%, 10% confidence
# interval for the probabilistic model and for all the quantile models.
y_mean, y_std = all_models["mse"].predict(X_train, return_std=True)
y_dist_train = all_models["mse"].sample(y_mean, y_std, n_estimates=1_000)
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
y_dist_test = all_models["mse"].sample(y_mean, y_std, n_estimates=1_000)

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
# distributions and find out which distribution fits our data the best.
#
# We already saw that the Normal distribution gave quite good results on the training
# set, but we've seen that the data is relatively skewed, so perhaps a distribution
# that allows skewed output response could be more appropriate. One such distribution
# is provided by default, the Gumbel distribution
# (https://en.wikipedia.org/wiki/Gumbel_distribution). This is a distribution commonly
# use in Extreme Value Theory to model extreme values. We don't need to retrain
# our model, we can just set the distribution parameter directly and create a
# new sample of outputs.
y_mean, y_std = all_models["mse"].predict(xx, return_std=True)
y_dist_test_normal = all_models["mse"].sample(y_mean, y_std, n_estimates=1_000)
all_models["mse"].distribution = "gumbel"
y_dist_test_gumbel = all_models["mse"].sample(y_mean, y_std, n_estimates=1_000)
# Based on our 1_000 estimates, we can obtain quantile estimates for every sample.
y_lower_normal = np.quantile(y_dist_test_normal, q=0.05, axis=0)
y_upper_normal = np.quantile(y_dist_test_normal, q=0.95, axis=0)
y_lower_gumbel = np.quantile(y_dist_test_gumbel, q=0.05, axis=0)
y_upper_gumbel = np.quantile(y_dist_test_gumbel, q=0.95, axis=0)

# Plot the outputs of the model with Normal distribution
fig, ax = plt.subplots(1, 2, figsize=(20, 10))
ax[0].plot(xx, f(xx), "g:", linewidth=3, label=r"$f(x) = x\,\sin(x)$")
ax[0].plot(X_test, y_test, "b.", markersize=10, label="Test observations")
ax[0].plot(xx, y_mean, "r-", label="Predicted mean")
ax[0].plot(xx, y_upper_normal, "k-")
ax[0].plot(xx, y_lower_normal, "k-")
ax[0].fill_between(
    xx.ravel(),
    y_lower_normal,
    y_upper_normal,
    alpha=0.4,
    label="Predicted 90% interval",
)
ax[0].set_xlabel("$x$")
ax[0].set_ylabel("$f(x)$")
ax[0].set_ylim(-10, 25)
ax[0].legend(loc="upper left")
ax[0].set_title("Normal distribution")

# Plot the outputs of the model with Gumbel distribution
ax[1].plot(xx, f(xx), "g:", linewidth=3, label=r"$f(x) = x\,\sin(x)$")
ax[1].plot(X_test, y_test, "b.", markersize=10, label="Test observations")
ax[1].plot(xx, y_mean, "r-", label="Predicted mean")
ax[1].plot(xx, y_upper_gumbel, "k-")
ax[1].plot(xx, y_lower_gumbel, "k-")
ax[1].fill_between(
    xx.ravel(),
    y_lower_gumbel,
    y_upper_gumbel,
    alpha=0.4,
    label="Predicted 90% interval",
)
ax[1].set_xlabel("$x$")
ax[1].set_ylabel("$f(x)$")
ax[1].set_ylim(-10, 25)
ax[1].legend(loc="upper left")
ax[1].set_title("Gumbel distribution")
plt.show()
# %%
# Comparing the Normal distribution with the Gumbel distribution we
# find that the Gumbel seems to track the shape of the skewed distribution
# a bit better, however it also seems to miss some of the points at the bottom
# of the curve.
#
# Let's check the calibration of the confidence intervals on the test set
y_mean, y_std = all_models["mse"].predict(X_test, return_std=True)
all_models["mse"].distribution = "normal"
y_dist_test_normal = all_models["mse"].sample(y_mean, y_std, n_estimates=1_000)
all_models["mse"].distribution = "gumbel"
y_dist_test_gumbel = all_models["mse"].sample(y_mean, y_std, n_estimates=1_000)
for i in range(n_quantile_pairs):
    lower_quantile = lower_quantiles[i]
    upper_quantile = upper_quantiles[i]
    interval_size = 100 * (upper_quantile - lower_quantile)
    cov_fraction_normal = coverage_fraction(
        y_test,
        y_low=np.quantile(y_dist_test_normal, q=lower_quantile, axis=0),
        y_high=np.quantile(y_dist_test_normal, q=upper_quantile, axis=0),
    )
    cov_fraction_gumbel = coverage_fraction(
        y_test,
        y_low=np.quantile(y_dist_test_gumbel, q=lower_quantile, axis=0),
        y_high=np.quantile(y_dist_test_gumbel, q=upper_quantile, axis=0),
    )
    print(
        f"Coverage @{interval_size:0.0f}% confidence interval: \n              Normal"
        f" distribution   : {cov_fraction_normal*100:0.0f}% \n             "
        f" Gumbel distribution   : {cov_fraction_gumbel*100:0.0f}%"
    )
# %%
# The calibration of the Gumbel is slightly worse than the Normal. Let's see
# whether we can further improve performance by tuning hyperparameters of
# our model.
#
# Tuning hyperparameters of our probabilistic predictions
# -------------------------------------------------------
#
# The probabilistic gradient boosting model introduces two key hyperparameters:
# 1. `tree_correlation` This value embodies the correlation we assume to exist
# between subsequent trees in the tree ensemble, and it is used to combine the
# variances of all the leaf values of every separate tree. By default, it is
# set at np.log10(n_training_samples) / 100.
# 2. `distribution` This is our distribution of choice that we choose to
# sample from.
#
# We can optimize these parameters in a CV procedure. In order to optimize
# the probabilistic estimate rather than the point estimate, we have to
# create a number of helper functions. We will evaluate the probabilistic
# forecast using Continuous Ranked Probability Score (CRPS). This score
# is a continous generalization of the Pinball loss, and can be used
# to evaluate probabilistic estimates (i.e., a cloud of estimates).
from sklearn.model_selection import GridSearchCV
from sklearn.metrics._scorer import _BaseScorer


def make_probabilistic_scorer(
    score_func,
    greater_is_better=True,
    **kwargs,
):

    sign = 1 if greater_is_better else -1
    return _ProbabilisticScorer(score_func, sign, kwargs)


class _ProbabilisticScorer(_BaseScorer):
    def _score(self, method_caller, estimator, X, y_true, sample_weight=None):
        y_pred, y_std = method_caller(estimator, "predict", X, return_std=True)
        yhat_dist = method_caller(estimator, "sample", y_pred, y_std, n_estimates=1_000)

        return self._sign * self._score_func(y_true, yhat_dist)


def crps_score(y, yhat_dist):
    """Calculate the empirical Continuously Ranked Probability Score (CRPS)
    for a set of forecasts for a number of samples (lower is better).

    Based on `crps_ensemble` from `properscoring`
    (https://pypi.org/project/properscoring/)

    :param yhat_dist: forecasts for each sample of size [n_forecasts x n_samples].
    :type yhat_dist: np.array
    :param y: ground truth value of each sample of size [n_samples].
    :type y: np.array

    :return: CRPS score for each sample
    :rtype: np.array

    See also:
    :arxiv:`O.Sprangers, S. Schelter, M. de Rijke, (2021) Probabilistic
    Gradient Boosting Machines, <2106.01682>.`

    """
    n_forecasts = yhat_dist.shape[0]
    n_samples = yhat_dist.shape[1]
    crps = np.zeros_like(y)
    # Loop over the samples
    for sample in range(n_samples):
        # Sort the forecasts in ascending order
        yhat_dist_sorted = np.sort(yhat_dist[:, sample])
        y_cdf = 0.0
        yhat_cdf = 0.0
        yhats_prev = 0.0
        ys = y[sample]
        # Loop over the forecasts per sample
        for yhats in yhat_dist_sorted:
            flag = (y_cdf == 0) * (ys < yhats) * 1.0
            crps[sample] += flag * (
                ((ys - yhats_prev) * yhat_cdf**2)
                + ((yhats - ys) * (yhat_cdf - 1) ** 2)
            )
            y_cdf += flag
            crps[sample] += (1 - flag) * (
                (yhats - yhats_prev) * (yhat_cdf - y_cdf) ** 2
            )
            yhat_cdf += 1 / n_forecasts
            yhats_prev = yhats

        # In case y_cdf == 0 after the loop
        flag = (y_cdf == 0) * 1.0
        crps[sample] += flag * (ys - yhats)

    return np.mean(crps)


score_func = make_probabilistic_scorer(crps_score, greater_is_better=False)

# Define the parameter grid.
param_grid = dict(
    learning_rate=[0.05, 0.1, 0.2],
    max_depth=[2, 5, 10],
    min_samples_leaf=[1, 5, 10, 20],
    tree_correlation=np.linspace(-0.05, 0.05, 10),
    distribution=["normal", "gumbel"],
)
gbr = HistGradientBoostingRegressor(random_state=0, with_variance=True)
search = GridSearchCV(
    gbr,
    param_grid,
    n_jobs=-1,
    scoring=score_func,
).fit(X_train, y_train)

# We subsequently fit and predict using the best found parameters
all_models["mse"].distribution = "normal"
y_mean, y_std = all_models["mse"].predict(xx, return_std=True)
y_dist_unoptimized = all_models["mse"].sample(y_mean, y_std, n_estimates=1_000)
y_mean, y_std = search.best_estimator_.predict(xx, return_std=True)
y_dist_optimized = search.best_estimator_.sample(y_mean, y_std, n_estimates=1_000)
# Based on our 1_000 estimates, we can obtain quantile estimates for every sample.
y_lower_unoptimized = np.quantile(y_dist_unoptimized, q=0.05, axis=0)
y_upper_unoptimized = np.quantile(y_dist_unoptimized, q=0.95, axis=0)
y_lower_optimized = np.quantile(y_dist_optimized, q=0.05, axis=0)
y_upper_optimized = np.quantile(y_dist_optimized, q=0.95, axis=0)

# Plot the outputs of the model with the default parameters
fig, ax = plt.subplots(1, 2, figsize=(20, 10))
ax[0].plot(xx, f(xx), "g:", linewidth=3, label=r"$f(x) = x\,\sin(x)$")
ax[0].plot(X_test, y_test, "b.", markersize=10, label="Test observations")
ax[0].plot(xx, y_mean, "r-", label="Predicted mean")
ax[0].plot(xx, y_upper_unoptimized, "k-")
ax[0].plot(xx, y_lower_unoptimized, "k-")
ax[0].fill_between(
    xx.ravel(),
    y_lower_unoptimized,
    y_upper_unoptimized,
    alpha=0.4,
    label="Predicted 90% interval",
)
ax[0].set_xlabel("$x$")
ax[0].set_ylabel("$f(x)$")
ax[0].set_ylim(-10, 25)
ax[0].legend(loc="upper left")
ax[0].set_title("Unoptimized hyperparameters")

# Plot the outputs of the model with the optimized parameters
ax[1].plot(xx, f(xx), "g:", linewidth=3, label=r"$f(x) = x\,\sin(x)$")
ax[1].plot(X_test, y_test, "b.", markersize=10, label="Test observations")
ax[1].plot(xx, y_mean, "r-", label="Predicted mean")
ax[1].plot(xx, y_upper_optimized, "k-")
ax[1].plot(xx, y_lower_optimized, "k-")
ax[1].fill_between(
    xx.ravel(),
    y_lower_optimized,
    y_upper_optimized,
    alpha=0.4,
    label="Predicted 90% interval",
)
ax[1].set_xlabel("$x$")
ax[1].set_ylabel("$f(x)$")
ax[1].set_ylim(-10, 25)
ax[1].legend(loc="upper left")
ax[1].set_title("Optimized hyperparameters")
plt.show()
# %%
# We can see that our optimized probabilistic forecast tracks the
# distribution of output values a bit better.
#
# Let's check the calibration of the confidence intervals on the test set
y_mean, y_std = all_models["mse"].predict(X_test, return_std=True)
y_dist_test_unoptimized = all_models["mse"].sample(y_mean, y_std, n_estimates=1_000)
y_mean, y_std = search.best_estimator_.predict(X_test, return_std=True)
y_dist_test_optimized = search.best_estimator_.sample(y_mean, y_std, n_estimates=1_000)
for i in range(n_quantile_pairs):
    lower_quantile = lower_quantiles[i]
    upper_quantile = upper_quantiles[i]
    interval_size = 100 * (upper_quantile - lower_quantile)
    cov_fraction_normal = coverage_fraction(
        y_test,
        y_low=np.quantile(y_dist_test_unoptimized, q=lower_quantile, axis=0),
        y_high=np.quantile(y_dist_test_unoptimized, q=upper_quantile, axis=0),
    )
    cov_fraction_gumbel = coverage_fraction(
        y_test,
        y_low=np.quantile(y_dist_test_optimized, q=lower_quantile, axis=0),
        y_high=np.quantile(y_dist_test_optimized, q=upper_quantile, axis=0),
    )
    print(
        f"Coverage @{interval_size:0.0f}% confidence interval: \n             "
        f" Unoptimized hyperparameters : {cov_fraction_normal*100:0.0f}% \n            "
        f"  Optimized hyperparameters   : {cov_fraction_gumbel*100:0.0f}%"
    )

# The coverage overall seems quite similar compared to our base case, except
# at the smaller confidence intervals, where the optimized version seems
# better.
#
# Let's check the overall CRPS score of our baseline model compared to our
# optimized model

crps_unoptimized = crps_score(y_test, y_dist_test_unoptimized)
crps_optimized = crps_score(y_test, y_dist_test_optimized)

# We find that the optimized CRPS is lower than the unoptimized CRPS, which
# means that our optimized model better captures the true distribution than
# the unoptimized version in terms of CRPS, which is what we aimed for!
