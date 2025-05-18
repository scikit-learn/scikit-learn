"""
=============================================
Regression Error Characteristic curve display
=============================================

We illustrate how to display the regression error characteristic curve, and how to use it
to compare various regression models.
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
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=80, random_state=42, shuffle=True
)

# %%
# Linear regression model
# -----------------------
#
# We create a linear regression model and fit it on the training after standard scaling
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

lr_estimator = make_pipeline(StandardScaler(), LinearRegression())
lr_estimator.fit(X_train, y_train)

# %%
# Display REC curve using the test set. See that our linear regressor's REC curve dominates
# the curve of the default constant predictor - the median. This is an indicator that our model is doing something
# "reasonable".
from sklearn.metrics import RecCurveDisplay

RecCurveDisplay.from_estimator(lr_estimator, X_test, y_test, name="Linear regression")


# %%
# Compare two REC curves of linear regression vs ridge regression for polnyomial features. We can see that the curve
# of the ridge regressor with polynomial features dominates the one of the linear regressor. Meaning, for any error
# tolerance, the Poly-Ridge model has more samples below this tolerance.
import matplotlib.pyplot as plt
from sklearn.linear_model import RidgeCV
from sklearn.metrics import root_mean_squared_error, mean_absolute_error
from sklearn.preprocessing import PolynomialFeatures


ridge_poly_estimator = make_pipeline(
    StandardScaler(), PolynomialFeatures(2, include_bias=False), RidgeCV()
)
ridge_poly_estimator.fit(X_train, y_train)

pred_lr = lr_estimator.predict(X_test)
pred_ridge_poly = ridge_poly_estimator.predict(X_test)

lr_metrics = f"RMSE = {root_mean_squared_error(pred_lr, y_test):.2f}, MAE = {mean_absolute_error(pred_lr, y_test):.2f}"
ridge_poly_metrics = f"RMSE = {root_mean_squared_error(pred_ridge_poly, y_test):.2f}, MAE = {mean_absolute_error(pred_ridge_poly, y_test):.2f}"

fig, ax = plt.subplots()
RecCurveDisplay.from_predictions(
    y_test,
    pred_lr,
    ax=ax,
    name=f"Linear regression ({lr_metrics})",
    plot_const_predictor=False,
)
RecCurveDisplay.from_predictions(
    y_test,
    pred_ridge_poly,
    ax=ax,
    name=f"Ridge Poly ({ridge_poly_metrics})",
    plot_const_predictor=False,
)
fig.show()

# %%
