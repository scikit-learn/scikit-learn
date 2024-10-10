"""
=========================
Train error vs Test error
=========================

Illustration of how the performance of an estimator on unseen data (test data)
is not the same as the performance on training data. As the regularization
increases the performance on train decreases while the performance on test
is optimal within a range of values of the regularization parameter.
The example with an Elastic-Net regression model and the performance is
measured using the explained variance a.k.a. R^2.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Generate sample data
# --------------------
import numpy as np

from sklearn import linear_model
from sklearn.datasets import make_regression
from sklearn.model_selection import train_test_split

n_samples_train, n_samples_test, n_features = 75, 150, 500
X, y, coef = make_regression(
    n_samples=n_samples_train + n_samples_test,
    n_features=n_features,
    n_informative=50,
    shuffle=False,
    noise=1.0,
    coef=True,
)
X_train, X_test, y_train, y_test = train_test_split(
    X, y, train_size=n_samples_train, test_size=n_samples_test, shuffle=False
)
# %%
# Compute train and test errors
# -----------------------------
alphas = np.logspace(-5, 1, 60)
enet = linear_model.ElasticNet(l1_ratio=0.7, max_iter=10000)
train_errors = list()
test_errors = list()
for alpha in alphas:
    enet.set_params(alpha=alpha)
    enet.fit(X_train, y_train)
    train_errors.append(enet.score(X_train, y_train))
    test_errors.append(enet.score(X_test, y_test))

i_alpha_optim = np.argmax(test_errors)
alpha_optim = alphas[i_alpha_optim]
print("Optimal regularization parameter : %s" % alpha_optim)

# Estimate the coef_ on full data with optimal regularization parameter
enet.set_params(alpha=alpha_optim)
coef_ = enet.fit(X, y).coef_

# %%
# Plotting Validation Curves
# -------------------------------------------------------------
# In this plot, you can see the training and validation scores
# of the ElasticNet model for different values of regularization
# parameter alpha. As can be inferred from the plot, for very low values
# of alpha (close to zero), the regularization is weak, meaning the model
# fits the training data very closely, leading to high training scores but lower
# validation scores. This is a case of overfitting, where the model captures
# noise in the training data rather than the underlying pattern.
#
# Using the ``ValidationCurveDisplay`` class helps by automating the plotting of
# trainingand validation scores across a range of alpha values, eliminating the
# need for manual iteration and plotting, and providing a clear, consistent
# visualization of model performance.

import matplotlib.pyplot as plt

from sklearn.model_selection import ValidationCurveDisplay

alphas = np.logspace(-5, 1, 60)

disp = ValidationCurveDisplay.from_estimator(
    enet,
    X_train,
    y_train,
    param_name="alpha",
    param_range=alphas,
    scoring="r2",
    n_jobs=2,
    score_type="both",
)

disp.ax_.set_title("Validation Curve for ElasticNet (R^2 Score)")
disp.ax_.set_xlabel(r"alpha (regularization strength)")
disp.ax_.set_ylabel("R^2 Score")
disp.ax_.set_ylim(-1.0, 1.2)
disp.ax_.vlines(
    alpha_optim,
    disp.ax_.get_ylim()[0],
    np.max(test_errors),
    color="k",
    linewidth=3,
    label="Optimum on test",
)
disp.ax_.legend(loc="lower right")

plt.show()

# %%
# Plotting Performance Comparison Curves
# -------------------------------------------------------------
# This plot compares the true coefficients (coef) with the estimated coefficients
# (coef_) from the model. It visually helps assess how well the model has captured the
# underlying patterns in the data.

plt.plot(coef, label="True coef")
plt.plot(coef_, label="Estimated coef")
plt.legend()
plt.title("True vs Estimated Coefficients")
plt.xlabel("Feature Index")
plt.ylabel("Coefficient Value")
plt.subplots_adjust(0.09, 0.04, 0.94, 0.94, 0.26, 0.26)
plt.show()
