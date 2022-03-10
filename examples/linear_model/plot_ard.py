"""
=======================================
Performance of Linear Regression Models
=======================================

Compare the performance of :ref:`automatic_relevance_determination` and
:ref:`bayesian_ridge_regression`.

Compared to the OLS (ordinary least squares) estimator, the coefficient weights
using a Bayesian Ridge are slightly shifted toward zeros, which stabilises them.
The ARD regression sets some of the non-informative cofficients exactly to zero,
while shifting some them closer to zero. Some non-informative coefficients are
still present and retain the order of magnitude of the ground truth weights.

The estimation of the model is done in both cases by iteratively maximizing the
marginal log-likelihood of the observations.

We also plot predictions and uncertainties for ARD and Bayesian Ridge
regressions using a polynomial feature expansion to fit a non-linear feature.
Notice that the ARD regression captures the ground truth the best, but due to
the limitations of a polynomial regression, both models fail when extrapolating.

"""

# %%
# Generating simulated data with Gaussian weights
from sklearn.datasets import make_regression

n_features = 100
X, y, w = make_regression(
    n_samples=100,
    n_features=n_features,
    n_informative=10,
    noise=8,
    coef=True,
    random_state=42,
)

# %%
# Fit the regressions
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
        "true weights": w,
    }
)
# %%
# Plot the true and estimated weights
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm

plt.figure(figsize=(10, 6))
ax = sns.heatmap(
    df.T,
    norm=LogNorm(),
    cbar_kws={"label": "weight value"},
    cmap="rocket_r",
)
plt.ylabel("predictor")
plt.xlabel("coeficients")
plt.tight_layout(rect=(0, 0, 1, 0.95))
_ = plt.title("Weights as a function of the model")
# %%
# Plot the marginal log-likelihood
plt.plot(ard.scores_, color="navy", label="Polynomial ARD")
plt.plot(brr.scores_, color="red", label="Polynomial BayesianRidge")
plt.ylabel("Score")
plt.xlabel("Iterations")
plt.xlim(1, 30)
plt.legend()
_ = plt.title("Convergence of marginal log-likelihood")
# %%
# Plotting polynomial regressions with standard deviations
import numpy as np
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import PolynomialFeatures, StandardScaler

rng = np.random.RandomState(0)
n_samples = 100
data_max, data_min = 0, 10
len_data = data_max - data_min
# sort the data to make plotting easier later
data = np.sort(rng.rand(n_samples) * len_data - len_data)
noise = rng.normal(0, 1, n_samples) * 1.35
target = np.sqrt(data) * np.sin(data) + noise
full_data = pd.DataFrame({"input_feature": data, "target": target})
# X should be 2D for sklearn: (n_samples, n_features)
data = data.reshape((-1, 1))

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
ard_poly.fit(data, target)
brr_poly.fit(data, target)
y_ard, y_ard_std = ard_poly.predict(data, return_std=True)
y_brr, y_brr_std = brr_poly.predict(data, return_std=True)

ax = sns.scatterplot(
    data=full_data, x="input_feature", y="target", color="black", alpha=0.75
)
ax.plot(data, target - noise, color="black", label="Ground Truth")
ax.plot(data, y_brr, color="red", label="Polynomial BayesianRidge")
ax.plot(data, y_ard, color="navy", label="Polynomial ARD")
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
_ = ax.set_title("Polynomial fit of non-linear feature")
