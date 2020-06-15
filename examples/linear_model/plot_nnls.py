"""
==========================
Non-negative least squares
==========================

In this example, we fit a linear model with positive constraints on the
regression parameters and compare the estimated parameters to a classic linear
regression.
"""
print(__doc__)
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score


# Generate some sparse data
np.random.seed(42)

n_samples, n_features = 200, 50
X = np.random.randn(n_samples, n_features)
coef = 3 * np.random.randn(n_features)
# Threshold coefficients to render them non-negative
coef[coef < 0] = 0
y = np.dot(X, coef)

# Add some noise
y += 5 * np.random.normal(size=(n_samples, ))

# Split the data in train set and test set
X_train, y_train = X[:n_samples // 2], y[:n_samples // 2]
X_test, y_test = X[n_samples // 2:], y[n_samples // 2:]

# %%
from sklearn.linear_model import LinearRegression

# %%
# Fit the Non-Negative least squares.
reg_nnls = LinearRegression(positive=True)
y_pred_nnls = reg_nnls.fit(X_train, y_train).predict(X_test)
r2_score_nnls = r2_score(y_test, y_pred_nnls)
print("NNLS R2 score", r2_score_nnls)

# %%
# Fit an OLS.
reg_ols = LinearRegression()
y_pred_ols = reg_ols.fit(X_train, y_train).predict(X_test)
r2_score_ols = r2_score(y_test, y_pred_ols)
print("OLS R2 score", r2_score_ols)


# %%
# Comparing the regression parameters between OLS and NNLS, we can observe
# they are highly correlated, but the non-negative constraint shrink some to
# 0. The Non-Negative Least square inherently yield sparse results.

fig, ax = plt.subplots()
ax.plot(reg_ols.coef_, reg_nnls.coef_, linewidth=0, marker=".")
ax.set_xlabel("OLS regression parameters", fontweight="bold")
ax.set_ylabel("NNLS regression parameters", fontweight="bold")
