# Author: Soumya
"""
==============================================
Robust regression using HuberRegressor
==============================================

This example compares HuberRegressor and LinearRegression on a dataset
with strong outliers. HuberRegressor is more robust to outliers and
fits better on the inlier data.
"""

import numpy as np
import matplotlib.pyplot as plt

from sklearn.linear_model import HuberRegressor, LinearRegression
from sklearn.model_selection import train_test_split

# Generate random data
rng = np.random.RandomState(42)
X = 2 * rng.rand(100, 1) - 1
y = 3 * X.squeeze() + 0.5 * rng.randn(100)

# Add some strong outliers
X[:10] = 2 * rng.rand(10, 1) - 1
y[:10] += 10 * (0.5 - rng.rand(10))  # large noise

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

# Fit both models
huber = HuberRegressor().fit(X_train, y_train)
ols = LinearRegression().fit(X_train, y_train)

# Plot results
plt.scatter(X_train, y_train, color="gray", alpha=0.6, label="Train data")
plt.plot(X_test, huber.predict(X_test), color="red", label="Huber Regressor")
plt.plot(X_test, ols.predict(X_test), color="blue", linestyle="--", label="Linear Regression")
plt.legend()
plt.title("Comparison of Huber and Linear Regression")
plt.xlabel("X")
plt.ylabel("y")
plt.show(block=True)

