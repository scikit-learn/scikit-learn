"""
===================================
Early stopping of Gradient Boosting
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

# Authors: Vighnesh Birodkar <vighneshbirodkar@nyu.edu>
#          Raghav RV <rvraghav93@gmail.com>
#          Kushan Sharma <kushansharma1@gmail.com>
# License: BSD 3 clause

import matplotlib.pyplot as plt
from sklearn.datasets import fetch_california_housing
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
import time

# Load the California Housing Prices dataset
data = fetch_california_housing()
X, y = data.data[:600], data.target[:600]

# Split data into training and validation sets
X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=42)

# Initialize GradientBoostingRegressor without early stopping
start_time = time.time()
gbm_no_early_stopping = GradientBoostingRegressor(
    n_estimators=1000, max_depth=5, learning_rate=0.1, random_state=42
)
gbm_no_early_stopping.fit(X_train, y_train)
training_time_without = time.time() - start_time
estimators_without = gbm_no_early_stopping.n_estimators_

# Initialize GradientBoostingRegressor with early stopping
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

# Create arrays to store training and validation errors
train_errors_without = []
val_errors_without = []

train_errors_with = []
val_errors_with = []

# Evaluate errors for each iteration
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

# Create plots
plt.figure(figsize=(15, 6))

# Plot training error vs. iterations
plt.subplot(1, 3, 1)
plt.plot(train_errors_without, label="Training Error (No Early Stopping)")
plt.plot(train_errors_with, label="Training Error (With Early Stopping)")
plt.xlabel("Boosting Iterations")
plt.ylabel("MSE (Training)")
plt.legend()

# Plot validation error vs. iterations
plt.subplot(1, 3, 2)
plt.plot(val_errors_without, label="Validation Error (No Early Stopping)")
plt.plot(val_errors_with, label="Validation Error (With Early Stopping)")
plt.xlabel("Boosting Iterations")
plt.ylabel("MSE (Validation)")
plt.legend()

# Plot training time comparison
plt.subplot(1, 3, 3)
training_times = [training_time_without, training_time_with]
labels = ["No Early Stopping", "With Early Stopping"]
bars = plt.bar(labels, training_times)
plt.ylabel("Training Time (s)")

# Add number of estimators on top of the bars
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
