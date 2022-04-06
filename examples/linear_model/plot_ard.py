"""
====================================
Comparing Linear Bayesian Regressors
====================================

This example compares two different bayesian regressors:

 - a :ref:`automatic_relevance_determination`
 - a :ref:`bayesian_ridge_regression`

In addition, we will also use an :ref:`ordinary_least_squares` (OLS) model as a
baseline for comparing the models' coefficients with respect to the true
coefficients.

The estimation of the model is done in both cases by iteratively maximizing the
marginal log-likelihood of the observations.

We also plot predictions and uncertainties for the ARD and the Bayesian Ridge
regressions using a polynomial feature expansion to fit a non-linear
relationship between `X` and `y`. Notice that the ARD regression captures the
ground truth the best, but due to the limitations of a polynomial regression,
both models fail when extrapolating.

"""

# %%
# Generate synthetic dataset
# --------------------------
#
# We generate a dataset where `X` and `y` are linearly linked: 10 of the
# features of `X` will be used to generate `y`. The other features are not
# useful at predicting `y`. In addition, we built a system where n_samples ==
# n_features. Such setting is challenging for an OLS model and leads to
# arbitrary large weights. Having a prior on the weights and a penalty
# alleviates the problem. Finally, gaussian noise is added.

from sklearn.datasets import make_regression

X, y, true_weights = make_regression(
    n_samples=100,
    n_features=100,
    n_informative=10,
    noise=8,
    coef=True,
    random_state=42,
)

# %%
# Fit the regressors
# ------------------
#
# We now fit both Bayesian models and the OLS to later compare the models'
# coefficients.

import pandas as pd
from sklearn.linear_model import ARDRegression, LinearRegression, BayesianRidge

olr = LinearRegression().fit(X, y)
brr = BayesianRidge(compute_score=True, n_iter=30).fit(X, y)
ard = ARDRegression(compute_score=True, n_iter=30).fit(X, y)
df = pd.DataFrame(
    {
        "LinearRegression": olr.coef_,
        "BayesianRidge": brr.coef_,
        "ARDRegression": ard.coef_,
        "Weights of true generative process": true_weights,
    }
)

# %%
# Plot the true and estimated coefficients
# ----------------------------------------
#
# Now we compared the coefficients of each model with the weights of
# the true generative model.
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm

plt.figure(figsize=(10, 6))
ax = sns.heatmap(
    df.T,
    norm=LogNorm(),
    cbar_kws={"label": "coefficients' values"},
    cmap="viridis_r",
)
plt.ylabel("linear model")
plt.xlabel("coefficients")
plt.tight_layout(rect=(0, 0, 1, 0.95))
_ = plt.title("Models' coefficients")

# %%
# Compared to the OLS estimator, the coefficients using a Bayesian Ridge
# regression are slightly shifted toward zeros, which stabilises them. The ARD
# regression provides a sparser solution: some of the non-informative
# coefficients are set exactly to zero, while shifting others closer to zero.
# Some non-informative coefficients are still present and retain the order of
# magnitude of the ground truth coefficients. Due to the added noise, none of
# the models recover the true weights. Indeed, all models always have more than
# 10 non-zero coefficients.

# %%
# Plot the marginal log-likelihood
# --------------------------------
import numpy as np

ard_scores = -np.array(ard.scores_)
brr_scores = -np.array(brr.scores_)
plt.plot(ard_scores, color="navy", label="ARD with polynomial features")
plt.plot(brr_scores, color="red", label="BayesianRidge with polynomial features")
plt.ylabel("Log-likelihood")
plt.xlabel("Iterations")
plt.xlim(1, 30)
plt.legend()
_ = plt.title("Models log-likelihood")

# %%
# Plotting polynomial regressions with std errors of the scores
# -------------------------------------------------------------

from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler

rng = np.random.RandomState(0)
n_samples = 110

# sort the data to make plotting easier later
X = np.sort(-10 * rng.rand(n_samples) + 10)
noise = rng.normal(0, 1, n_samples) * 1.35
y = np.sqrt(X) * np.sin(X) + noise
full_data = pd.DataFrame({"input_feature": X, "target": y})
X = X.reshape((-1, 1))

# `fit_intercept=True` for `ARDRegression()` and `BayesianRidge()`,
# then `PolynomialFeatures` should not introduce the bias feature
ard_poly = make_pipeline(
    PolynomialFeatures(degree=10, include_bias=False),
    StandardScaler(),
    ARDRegression(),
)
brr_poly = make_pipeline(
    PolynomialFeatures(degree=10, include_bias=False),
    StandardScaler(),
    BayesianRidge(),
)
ard_poly.fit(X, y)
brr_poly.fit(X, y)
y_ard, y_ard_std = ard_poly.predict(X, return_std=True)
y_brr, y_brr_std = brr_poly.predict(X, return_std=True)

ax = sns.scatterplot(
    data=full_data, x="input_feature", y="target", color="black", alpha=0.75
)
ax.plot(X, y - noise, color="black", label="Ground Truth")
ax.plot(X, y_brr, color="red", label="BayesianRidge with polynomial features")
ax.plot(X, y_ard, color="navy", label="ARD with polynomial features")
ax.fill_between(
    full_data["input_feature"],
    y_ard - y_ard_std,
    y_ard + y_ard_std,
    color="navy",
    alpha=0.3,
)
ax.fill_between(
    full_data["input_feature"],
    y_brr - y_brr_std,
    y_brr + y_brr_std,
    color="red",
    alpha=0.3,
)
ax.legend()
_ = ax.set_title("Polynomial fit of a non-linear feature")
