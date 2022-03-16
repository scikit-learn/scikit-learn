"""
==========================================
Comparison of Linear Models for Regression
==========================================

This example compares three kind of linear models for regression:

 - a :ref:`ordinary-least-squares`
 - a :ref:`automatic_relevance_determination`
 - a :ref:`bayesian_ridge_regression`

The comparisons are based on the models' coefficients with respect to the true
coefficients.

The estimation of the model is done in both cases by iteratively maximizing the
marginal log-likelihood of the observations.

We also plot predictions and uncertainties for the ARD and the Bayesian Ridge
regressions using a polynomial feature expansion to fit a non-linear feature.
Notice that the ARD regression captures the ground truth the best, but due to
the limitations of a polynomial regression, both models fail when extrapolating.

"""

# %%
# Generate synthetic dataset
# --------------------------
from sklearn.datasets import make_regression

X, y, w = make_regression(
    n_samples=100,
    n_features=100,
    n_informative=10,
    noise=8,
    coef=True,
    random_state=42,
)

# %%
# Fit the regressions
# -------------------
import pandas as pd
from sklearn.linear_model import ARDRegression, LinearRegression, BayesianRidge

olr = LinearRegression()
brr = BayesianRidge(compute_score=True, n_iter=30)
ard = ARDRegression(compute_score=True, n_iter=30)

olr.fit(X, y)
brr.fit(X, y)
ard.fit(X, y)

df = pd.DataFrame(
    {
        "LinearRegression": olr.coef_,
        "BayesianRidge": brr.coef_,
        "ARDRegression": ard.coef_,
        "ground truth": w,
    }
)

# %%
# Plot the true and estimated coefficients
# ----------------------------------------
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
# Compared to the OLS (ordinary least squares) estimator, the coefficients using
# a Bayesian Ridge regression are slightly shifted toward zeros, which
# stabilises them. The ARD regression sets some of the non-informative
# coefficients exactly to zero, while shifting some them closer to zero. Some
# non-informative coefficients are still present and retain the order of
# magnitude of the ground truth coefficients.

# %%
# Plot the marginal log-likelihood
# --------------------------------
plt.plot(ard.scores_, color="navy", label="ARD with polynomial features")
plt.plot(brr.scores_, color="red", label="BayesianRidge with polynomial features")
plt.ylabel("Score")
plt.xlabel("Iterations")
plt.xlim(1, 30)
plt.legend()
_ = plt.title("Models log-likelihood")

# %%
# Plotting polynomial regressions with std errors of the scores
# -------------------------------------------------------------
import numpy as np
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler

rng = np.random.RandomState(0)
n_samples = 100

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
