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
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=20, shuffle=False)

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

RecCurveDisplay.from_estimator(lr_estimator, X_test, y_test)


# %%
# Display REC curve of the model only, without the constant predictor. The linear regressor's REC curve also
# dominates the one of another constant predictor - the mean.
RecCurveDisplay.from_estimator(lr_estimator, X_test, y_test, constant_predictor="mean")


# %% Compare two REC curves of linear regression vs linear SVR. We can see different performance profiles
# of both estimators. The LinearSVR appears to have more samples with errors below 60, whereas
# linear regressor appears to have more samples than the SVR with errors below 20. This is despite both having
# almost the same summary metrics.
import matplotlib.pyplot as plt
from sklearn.svm import LinearSVR
from sklearn.metrics import root_mean_squared_error, mean_absolute_error


svr_estimator = make_pipeline(StandardScaler(), LinearSVR())
svr_estimator.fit(X_train, y_train)

pred_lr = lr_estimator.predict(X_test)
pred_svr = svr_estimator.predict(X_test)

lr_metrics = f"RMSE = {root_mean_squared_error(pred_lr, y_test):.2f}, MAE = {mean_absolute_error(pred_lr, y_test):.2f}"
svr_metrics = f"RMSE = {root_mean_squared_error(pred_svr, y_test):.2f}, MAE = {mean_absolute_error(pred_svr, y_test):.2f}"

fig, ax = plt.subplots()
RecCurveDisplay.from_predictions(
    pred_lr,
    y_test,
    ax=ax,
    name=f"Linear regression ({lr_metrics})",
    plot_const_predictor=False,
)
RecCurveDisplay.from_predictions(
    pred_svr,
    y_test,
    ax=ax,
    name=f"Linear SVR ({svr_metrics})",
    plot_const_predictor=False,
)
fig.show()
