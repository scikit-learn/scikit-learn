"""
=====================================
Simple Ordinary Least Squares Example
=====================================

This example shows how to use the simplest ordinary least squares (OLS) model
called :class:`~sklearn.linear_model.LinearRegression` in scikit-learn.

For this purpose, we use a single feature from the diabetes dataset and try to
predict the diabetes progression using this linear model. We therefore load the
diabetes dataset and split it into training and test sets.

Then, we fit the model on the training set and evaluate its performance on the test
set and finally visualize the results on the test set.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Data Loading and Preparation
# ----------------------------
#
# Load the diabetes dataset. For simplicity, we only keep a single feature in the data.
# Then, we split the data and target into training and test sets.
from sklearn.datasets import load_diabetes
from sklearn.model_selection import train_test_split

X, y = load_diabetes(return_X_y=True)
X = X[:, [2]]  # Use only one feature
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=20, shuffle=False)

# %%
# Linear regression model
# -----------------------
#
# We create a linear regression model and fit it on the training data.
from sklearn.linear_model import LinearRegression

regressor = LinearRegression().fit(X_train, y_train)

# %%
# Model evaluation
# ----------------
#
# We evaluate the model's performance on the test set using the mean squared error
# and the coefficient of determination.
from sklearn.metrics import mean_squared_error, r2_score

y_pred = regressor.predict(X_test)

print(f"Mean squared error: {mean_squared_error(y_test, y_pred):.2f}")
print(f"Coefficient of determination: {r2_score(y_test, y_pred):.2f}")

# %%
# Plotting the results
# --------------------
#
# Finally, we visualize the results on the test set.
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

ax.scatter(X_test, y_test, label="Data points")
ax.plot(X_test, y_pred, linewidth=3, color="tab:orange", label="Model predictions")
ax.set(xlabel="Feature", ylabel="Target", title="Linear Regression")
ax.legend()

plt.show()

# %%
# Conclusion
# ----------
#
# This example shows how to use the simplest linear model called
# :class:`~sklearn.linear_model.LinearRegression` in scikit-learn. For this purpose, we
# use a single feature from the diabetes dataset and try to predict the diabetes
# progression using this linear model. We therefore load the diabetes dataset and split
# it into training and test sets.
#
