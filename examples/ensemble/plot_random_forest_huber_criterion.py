"""
============================================================
Effect of Huber Loss Criterion on Random Forest Regression
============================================================

In this example we examine the effect of the Huber loss criterion for
robust regression using the RandomForestRegressor on a synthetic dataset
with outliers.

This Python script uses the RandomForestRegressor from the sklearn library to
create two models, one using the squared error loss function and the other
using the Huber loss function.

The script first generates a synthetic dataset using the `make_regression` function. It
then introduces some artificial outliers to the dataset.

The data is then split into training and testing sets using the `train_test_split`
function.

Two RandomForestRegressor models are trained on the training data - one with the
squared error loss function (default) and the other with the Huber loss function.

The models are then used to predict the output for the test data. The mean squared error
of the predictions is calculated using the `mean_squared_error` function.

Finally, the script plots the test data and the predictions of both models using
matplotlib. The plots show the test data in red, the predictions of the model with
squared error loss in green, and the predictions of the model with Huber loss
in purple. The mean squared error of each model is also displayed in the title of
the respective plot.

The purpose of this script is to compare the performance of the RandomForestRegressor
with different loss functions on a dataset with outliers.

"""


import matplotlib.pyplot as plt
import numpy as np

from sklearn.datasets import make_regression
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split

# Generate a synthetic dataset
X, y = make_regression(
    n_samples=200,
    n_features=1,
    tail_strength=0.9,
    effective_rank=1,
    n_informative=1,
    noise=3,
    bias=20,
    random_state=1,
)
# add some artificial outliers
np.random.seed(1)
for i in range(20):
    factor = np.random.randint(2, 4)
    if np.random.rand() > 0.5:
        X[i] += factor * X.std()
    else:
        X[i] -= factor * X.std()

# Random Forest parameters
N_ESIMATORS = 200
MAX_DEPTH = 8
DELTA = 0.25

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# Train a model with squared loss (default)
model_squared_loss = RandomForestRegressor(
    max_depth=MAX_DEPTH, n_estimators=N_ESIMATORS, random_state=42
)
model_squared_loss.fit(X_train, y_train)
y_pred_squared_loss = model_squared_loss.predict(X_test)

# Train a model with Huber loss
model_huber_loss = RandomForestRegressor(
    criterion="huber",
    delta=DELTA,
    max_depth=MAX_DEPTH,
    n_estimators=N_ESIMATORS,
    random_state=42,
)
model_huber_loss.fit(X_train, y_train)
y_pred_huber_loss = model_huber_loss.predict(X_test)

# Calculate the mean squared error
mse_squared_loss = mean_squared_error(y_test, y_pred_squared_loss)
mse_huber_loss = mean_squared_error(y_test, y_pred_huber_loss)

# Plotting
plt.figure(figsize=(12, 6))
plt.suptitle(
    f"RandomForestRegresor(max_depth={MAX_DEPTH}, n_estimators={N_ESIMATORS})\n"
)

plt.subplot(1, 2, 1)
plt.scatter(X_test, y_test, color="red", label="Test data")
plt.scatter(X_test, y_pred_squared_loss, color="green", label="Prediction")
plt.title('criterion="squared_error"\nMSE: {:.2f}'.format(mse_squared_loss))
plt.legend()

plt.subplot(1, 2, 2)
plt.scatter(X_test, y_test, color="red", label="Test data")
plt.scatter(X_test, y_pred_huber_loss, color="purple", label="Prediction")
plt.title('criterion="huber", delta={:.2f}\nMSE: {:.2f}'.format(DELTA, mse_huber_loss))
plt.legend()

plt.subplots_adjust(top=0.85, bottom=0.15)

plt.show()
