"""
===================================
Early stopping in Gradient Boosting
===================================

Gradient Boosting is an ensemble technique that combines multiple weak
learners, typically regression trees, to create a robust and powerful
predictive model. It does so in an iterative fashion, where each new stage
(tree) corrects the errors of the previous ones.

Early stopping is a feature in Gradient Boosting that allows us to find
the optimal number of iterations required to build a model that generalizes
well to unseen data. The concept is simple: we set aside a portion of our
dataset as a validation set (specified using `validation_fraction`) to assess
the model's performance during training. As the model is iteratively built
with additional stages (trees), its performance on the validation set is
continuously monitored.

Early stopping becomes effective when the model's performance on the
validation set plateaus or worsens over a certain number of consecutive stages
(specified by `n_iter_no_change`) without a significant improvement. This
signals that the model has reached a point where further iterations may lead
to overfitting, and it's time to stop training.

In our example with the `GradientBoostingRegressor` model on the California
Housing Prices dataset, we have demonstrated the practical benefits of early
stopping:

- **Preventing Overfitting:** We showed how the validation error stabilizes
or starts to increase after a certain point, indicating that the model
generalizes better to unseen data. This is achieved by stopping the training
process before overfitting occurs.

- **Improving Training Efficiency:** We compared training times between models
with and without early stopping. The model with early stopping achieved
comparable accuracy while requiring significantly fewer estimators, resulting
in faster training.

The number of estimators (trees) in the final model, when early stopping is
applied, can be accessed using the `n_estimators_` attribute. Overall, early
stopping is a valuable tool to strike a balance between model performance and
efficiency in gradient boosting.

"""
# %%
# Data Preparation
# ----------------
# Loads and prepares the California Housing Prices dataset for training and evaluation.
# It subsets the dataset, splits it into training and validation sets.

# Authors: Vighnesh Birodkar <vighneshbirodkar@nyu.edu>
#          Raghav RV <rvraghav93@gmail.com>
#          Kushan Sharma <kushansharma1@gmail.com>
# License: BSD 3 clause

import time

import matplotlib.pyplot as plt

from sklearn.datasets import fetch_california_housing
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split

data = fetch_california_housing()
X, y = data.data[:600], data.target[:600]

X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=42)

# %%
# Model Training and Comparison
# -----------------------------
# Two :class:`~sklearn.ensemble.GradientBoostingRegressor` models are trained:
# one without early stopping and another with early stopping. The purpose is to
# compare their performance. It also calculates the training time and the
# `n_estimators_` used by both models.

start_time = time.time()
gbm_no_early_stopping = GradientBoostingRegressor(
    n_estimators=1000, max_depth=5, learning_rate=0.1, random_state=42
)
gbm_no_early_stopping.fit(X_train, y_train)
training_time_without = time.time() - start_time
estimators_without = gbm_no_early_stopping.n_estimators_

start_time = time.time()
gbm_with_early_stopping = GradientBoostingRegressor(
    n_estimators=1000,
    max_depth=5,
    learning_rate=0.1,
    validation_fraction=0.1,
    n_iter_no_change=10,
    random_state=42,
)
gbm_with_early_stopping.fit(X_train, y_train)
training_time_with = time.time() - start_time
estimators_with = gbm_with_early_stopping.n_estimators_

# %%
# Error Calculation
# -----------------
# The code calculates the :func:`~sklearn.metrics.mean_squared_error` for both
# training and validation datasets for the models trained in the previous section.
# It computes the errors for each boosting iteration. The purpose is to assess the
# performance and convergence of the models.

train_errors_without = []
val_errors_without = []

train_errors_with = []
val_errors_with = []

for i, (train_pred, val_pred) in enumerate(
    zip(
        gbm_no_early_stopping.staged_predict(X_train),
        gbm_no_early_stopping.staged_predict(X_val),
    )
):
    train_errors_without.append(mean_squared_error(y_train, train_pred))
    val_errors_without.append(mean_squared_error(y_val, val_pred))

for i, (train_pred, val_pred) in enumerate(
    zip(
        gbm_with_early_stopping.staged_predict(X_train),
        gbm_with_early_stopping.staged_predict(X_val),
    )
):
    train_errors_with.append(mean_squared_error(y_train, train_pred))
    val_errors_with.append(mean_squared_error(y_val, val_pred))

# %%
# Visualize Comparision
# ---------------------
# It includes three subplots:
# 1. Plotting training errors of both models over boosting iterations.
# 2. Plotting validation errors of both models over boosting iterations.
# 3. Creating a bar chart to compare the training times and the estimator used
# of the models with and without early stopping.

plt.figure(figsize=(15, 6))

plt.subplot(1, 3, 1)
plt.plot(train_errors_without, label="Training Error (No Early Stopping)")
plt.plot(train_errors_with, label="Training Error (With Early Stopping)")
plt.xlabel("Boosting Iterations")
plt.ylabel("MSE (Training)")
plt.legend()

plt.subplot(1, 3, 2)
plt.plot(val_errors_without, label="Validation Error (No Early Stopping)")
plt.plot(val_errors_with, label="Validation Error (With Early Stopping)")
plt.xlabel("Boosting Iterations")
plt.ylabel("MSE (Validation)")
plt.legend()

plt.subplot(1, 3, 3)
training_times = [training_time_without, training_time_with]
labels = ["No Early Stopping", "With Early Stopping"]
bars = plt.bar(labels, training_times)
plt.ylabel("Training Time (s)")

for bar, n_estimators in zip(bars, [estimators_without, estimators_with]):
    height = bar.get_height()
    plt.text(
        bar.get_x() + bar.get_width() / 2,
        height + 0.01,
        f"Estimators: {n_estimators}",
        ha="center",
        va="bottom",
    )

plt.tight_layout()
plt.show()
