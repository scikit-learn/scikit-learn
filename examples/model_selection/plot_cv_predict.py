"""
====================================
Plotting Cross-Validated Predictions
====================================

This example shows how to use
:func:`~sklearn.model_selection.cross_val_predict` together with
:class:`~sklearn.metrics.PredictionErrorDisplay` to visualize prediction
errors.
"""

# %%
# We will load the diabetes dataset and create an instance of a linear
# regression model.
from sklearn.datasets import load_diabetes
from sklearn.linear_model import LinearRegression

X, y = load_diabetes(return_X_y=True)
lr = LinearRegression()

# %%
# :func:`~sklearn.model_selection.cross_val_predict` returns an array of the
# same size as `y` where each entry is a prediction obtained by cross
# validation.
from sklearn.model_selection import cross_val_predict

y_pred = cross_val_predict(lr, X, y, cv=10)

# %%
# Since `cv=10`, it means that we trained 10 models and each of the model was
# used to predict on one of the 10 folds. We can now use the
# :class:`~sklearn.metrics.PredictionErrorDisplay` to visualize the
# prediction errors.
#
# On the left axis, we plot the true values of `y` vs. the predicted values
# given by the models. On the right axis, we plot the residuals (i.e. the
# difference between the true values and the predicted values) vs. the
# predicted values.
import matplotlib.pyplot as plt
from sklearn.metrics import PredictionErrorDisplay

fig, axs = plt.subplots(ncols=2, figsize=(10, 4))
PredictionErrorDisplay.from_predictions(
    y,
    y_pred=y_pred,
    kind="predictions",
    subsample=100,
    ax=axs[0],
    random_state=0,
)
axs[0].set_title("Predicted vs. True Values")
PredictionErrorDisplay.from_predictions(
    y, y_pred=y_pred, kind="residuals", subsample=100, ax=axs[1], random_state=0
)
axs[1].set_title("Residuals vs. Predicted Values")
fig.suptitle("Plotting cross-validated predictions")
plt.tight_layout()
plt.show()

# %%
# It is important to note that
# :func:`~sklearn.model_selection.cross_val_predict` is a tool that should be
# used for visualization but not for evaluating the performance of the
# predictive models by computing a metric on the aggregated predictions.
#
# A proper validation should use either
# :func:`~sklearn.model_selection.cross_val_score` or
# :func:`~sklearn.model_selection.cross_validate` instead.
