"""
=========================================================
Single estimator VS. Bagging: bias-variance decomposition
=========================================================

TODO

"""

import numpy as np
from matplotlib import pyplot as plt

from sklearn.ensemble import BaggingRegressor
from sklearn.tree import DecisionTreeRegressor

# Settings
n_repeat = 50       # Number of repetitions
n_train = 100       # Size of the training set
n_test = 5000       # Size of the test set
eps = 0.05          # Standard deviation of the noise

estimator = DecisionTreeRegressor(min_samples_split=1)
bagging = BaggingRegressor(estimator, n_estimators=50)

np.random.seed(0)

# Generate data
def f(x):
    return np.exp(-x ** 2) + 1.5 * np.exp(-(x - 2) ** 2)

def generate(n_samples, eps):
    X = np.random.rand(n_samples) * 10 - 5
    X = np.sort(X)
    y = f(X) + np.random.normal(0.0, eps, n_samples)
    X = X.reshape((n_samples, 1))

    return X, y

X_train = []
y_train = []

for i in range(n_repeat):
    X, y = generate(n_samples=n_train, eps=eps)
    X_train.append(X)
    y_train.append(y)

X_test, y_test = generate(n_samples=n_test, eps=eps)

# Build estimators and compute predictions
y_predict_estimator = np.zeros((n_test, n_repeat))
y_predict_bagging = np.zeros((n_test, n_repeat))

for i in xrange(n_repeat):
    estimator.fit(X_train[i], y_train[i])
    y_predict_estimator[:, i] = estimator.predict(X_test)
    bagging.fit(X_train[i], y_train[i])
    y_predict_bagging[:, i] = bagging.predict(X_test)

# Bias^2 + Variance decomposition of the mean squared error
y_error_estimator = np.zeros(n_test)
y_error_bagging = np.zeros(n_test)

for i in range(n_repeat):
    y_error_estimator += (f(X_test)[:, 0] - y_predict_estimator[:, i]) ** 2
    y_error_bagging += (f(X_test)[:, 0] - y_predict_bagging[:, i]) ** 2

y_error_estimator /= n_repeat
y_error_bagging /= n_repeat

y_bias_estimator = (f(X_test)[:, 0] - np.mean(y_predict_estimator, axis=1)) ** 2
y_var_estimator = np.var(y_predict_estimator, axis=1)
y_bias_bagging = (f(X_test)[:, 0] - np.mean(y_predict_bagging, axis=1)) ** 2
y_var_bagging = np.var(y_predict_bagging, axis=1)

print "Single estimator: %f (mse) = %f (bias^2) + %f (var)" % (np.mean(y_error_estimator), np.mean(y_bias_estimator), np.mean(y_var_estimator))
print "Bagging: %f (mse) = %f (bias^2) + %f (var)" % (np.mean(y_error_bagging), np.mean(y_bias_bagging), np.mean(y_var_bagging))

# Plot
plt.subplot(2, 1, 1)

for i in range(n_repeat):
    plt.plot(X_test, y_predict_estimator[:, i], "r", alpha=0.05)

plt.plot(X_test, f(X_test), "b", label="f(x)")
plt.plot(X_test, np.mean(y_predict_estimator, axis=1), "r", label="mean prediction (single estimator)")
plt.plot(X_test, np.mean(y_predict_bagging, axis=1), "-.r", label="mean prediction (bagging)")
plt.legend(loc="best", prop={"size": 9})

plt.subplot(2, 1, 2)

plt.plot(X_test, y_error_estimator, "r", label="mse (single estimator)")
plt.plot(X_test, y_bias_estimator, "b", label="bias^2 (single estimator")
plt.plot(X_test, y_var_estimator, "g", label="var (single estimator)")

plt.plot(X_test, y_error_bagging, "-.r", label="mse (bagging)")
plt.plot(X_test, y_bias_bagging, "-.b", label="bias^2 (bagging)")
plt.plot(X_test, y_var_bagging, "-.g", label="var (bagging)")
plt.legend(loc="best", prop={"size": 9})

plt.show()
