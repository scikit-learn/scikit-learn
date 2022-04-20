# -*- coding: utf-8 -*-
"""
======================================================
Effect of transforming the targets in regression model
======================================================

In this example, we give an overview of
:class:`~sklearn.compose.TransformedTargetRegressor`. We use two examples
to illustrate the benefit of transforming the targets before learning a linear
regression model. The first example uses synthetic data while the second
example is based on the Ames housing data set.

"""

# Author: Guillaume Lemaitre <guillaume.lemaitre@inria.fr>
# License: BSD 3 clause

print(__doc__)

# %%
# Synthetic example
###################
#
# A synthetic random regression dataset is generated. The targets ``y`` are
# modified by:
#
#   1. translating all targets such that all entries are
#      non-negative (by adding the absolute value of the lowest ``y``) and
#   2. applying an exponential function to obtain non-linear
#      targets which cannot be fitted using a simple linear model.
#
# Therefore, a logarithmic (`np.log1p`) and an exponential function
# (`np.expm1`) will be used to transform the targets before training a linear
# regression model and using it for prediction.
import numpy as np
from sklearn.datasets import make_regression

X, y = make_regression(n_samples=10000, noise=100, random_state=0)
y = np.expm1((y + abs(y.min())) / 200)
y_trans = np.log1p(y)

# %%
# Below we plot the probability density functions of the target
# before and after applying the logarithmic functions.
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

f, (ax0, ax1) = plt.subplots(1, 2)

ax0.hist(y, bins=100, density=True)
ax0.set_xlim([0, 2000])
ax0.set_ylabel("Probability")
ax0.set_xlabel("Target")
ax0.set_title("Target distribution")

ax1.hist(y_trans, bins=100, density=True)
ax1.set_ylabel("Probability")
ax1.set_xlabel("Target")
ax1.set_title("Transformed target distribution")

f.suptitle("Synthetic data", y=0.06, x=0.53)
f.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])

X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

# %%
# At first, a linear model will be applied on the original targets. Due to the
# non-linearity, the model trained will not be precise during
# prediction. Subsequently, a logarithmic function is used to linearize the
# targets, allowing better prediction even with a similar linear model as
# reported by the median absolute error (MAE).
from sklearn.metrics import mean_absolute_error, r2_score


def compute_score(y_true, y_pred):
    return {
        "R2 score": f"{r2_score(y_true, y_pred):.3f}",
        "MAE": f"{mean_absolute_error(y_true, y_pred):.3f}",
    }


# %%
from sklearn.compose import TransformedTargetRegressor
from sklearn.linear_model import RidgeCV
from sklearn.metrics import PredictionErrorDisplay

f, (ax0, ax1) = plt.subplots(1, 2, sharey=True)

ridge_cv = RidgeCV().fit(X_train, y_train)
y_pred_ridge = ridge_cv.predict(X_test)

ridge_cv_with_trans_target = TransformedTargetRegressor(
    regressor=RidgeCV(), func=np.log1p, inverse_func=np.expm1
).fit(X_train, y_train)
y_pred_ridge_with_trans_target = ridge_cv_with_trans_target.predict(X_test)

PredictionErrorDisplay.from_predictions(
    y_test,
    y_pred_ridge,
    scores=compute_score(y_test, y_pred_ridge),
    ax=ax0,
    scatter_kwargs={"alpha": 0.5},
)
PredictionErrorDisplay.from_predictions(
    y_test,
    y_pred_ridge_with_trans_target,
    scores=compute_score(y_test, y_pred_ridge_with_trans_target),
    ax=ax1,
    scatter_kwargs={"alpha": 0.5},
)

ax0.set_title("Ridge regression \n without target transformation")
ax1.set_title("Ridge regression \n with target transformation")
f.suptitle("Synthetic data", y=0.035)
f.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])

# %%
# Real-world data set
#####################
#
# In a similar manner, the Ames housing data set is used to show the impact
# of transforming the targets before learning a model. In this example, the
# target to be predicted is the selling price of each house.

# %%
from sklearn.datasets import fetch_openml
from sklearn.preprocessing import quantile_transform

ames = fetch_openml(name="house_prices", as_frame=True)
# Keep only numeric columns
X = ames.data.select_dtypes(np.number)
# Remove columns with NaN or Inf values
X = X.drop(columns=["LotFrontage", "GarageYrBlt", "MasVnrArea"])
# Let's the price being in k$
y = ames.target / 1000
y_trans = quantile_transform(
    y.to_frame(), n_quantiles=900, output_distribution="normal", copy=True
).squeeze()
# %%
# A :class:`~sklearn.preprocessing.QuantileTransformer` is used to normalize
# the target distribution before applying a
# :class:`~sklearn.linear_model.RidgeCV` model.
f, (ax0, ax1) = plt.subplots(1, 2)

ax0.hist(y, bins=100, density=True)
ax0.set_ylabel("Probability")
ax0.set_xlabel("Target")
ax0.set_title("Target distribution")

ax1.hist(y_trans, bins=100, density=True)
ax1.set_ylabel("Probability")
ax1.set_xlabel("Target")
ax1.set_title("Transformed target distribution")

f.suptitle("Ames housing data: selling price", y=0)
f.tight_layout()

# %%
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=1)

# %%
# The effect of the transformer is weaker than on the synthetic data. However,
# the transformation results in an increase in :math:`R^2` and large decrease
# of the MAE. The residual plot (predicted target - true target vs predicted
# target) without target transformation takes on a curved, 'reverse smile'
# shape due to residual values that vary depending on the value of predicted
# target. With target transformation, the shape is more linear indicating
# better model fit.
from sklearn.preprocessing import QuantileTransformer

f, (ax0, ax1) = plt.subplots(2, 2, sharey="row", figsize=(6.5, 8))

ridge_cv = RidgeCV().fit(X_train, y_train)
y_pred_ridge = ridge_cv.predict(X_test)

ridge_cv_with_trans_target = TransformedTargetRegressor(
    regressor=RidgeCV(),
    transformer=QuantileTransformer(n_quantiles=900, output_distribution="normal"),
).fit(X_train, y_train)
y_pred_ridge_with_trans_target = ridge_cv_with_trans_target.predict(X_test)

# plot the actual vs predicted values
PredictionErrorDisplay.from_predictions(
    y_test,
    y_pred_ridge,
    scores=compute_score(y_test, y_pred_ridge),
    ax=ax0[0],
    scatter_kwargs={"alpha": 0.5},
)
PredictionErrorDisplay.from_predictions(
    y_test,
    y_pred_ridge_with_trans_target,
    scores=compute_score(y_test, y_pred_ridge_with_trans_target),
    ax=ax0[1],
    scatter_kwargs={"alpha": 0.5},
)
ax0[0].set_title("Ridge regression \n without target transformation")
ax0[1].set_title("Ridge regression \n with target transformation")

# plot the residuals vs the predicted values
PredictionErrorDisplay.from_predictions(
    y_test,
    y_pred_ridge,
    scores=compute_score(y_test, y_pred_ridge),
    kind="residuals",
    ax=ax1[0],
    scatter_kwargs={"alpha": 0.5},
)
PredictionErrorDisplay.from_predictions(
    y_test,
    y_pred_ridge_with_trans_target,
    scores=compute_score(y_test, y_pred_ridge_with_trans_target),
    kind="residuals",
    ax=ax1[1],
    scatter_kwargs={"alpha": 0.5},
)
ax1[0].set_title("Ridge regression \n without target transformation")
ax1[1].set_title("Ridge regression \n with target transformation")

f.suptitle("Ames housing data: selling price", y=0)
plt.tight_layout()
plt.show()
