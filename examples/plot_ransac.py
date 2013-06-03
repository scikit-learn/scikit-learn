"""
===========================================
Robust linear model estimation using RANSAC
===========================================

In this example we see how to robustly fit a linear model to faulty data using
the RANSAC algorithm.

"""
import numpy as np
from matplotlib import pyplot as plt

from sklearn import linear_model
from sklearn.utils import ransac


# Set random seed for both equal data noise and equal random sample selection
np.random.seed(seed=1)

# Generate coordinates of line
X = np.arange(-200, 200)
y = 0.2 * X + 20
data = np.column_stack([X, y])

# Add faulty data
faulty = np.array(30 * [(180, -100)], dtype=np.double)
faulty += 5 * np.random.normal(size=faulty.shape)
data[:faulty.shape[0]] = faulty

# Add gaussian noise to coordinates
noise = np.random.normal(size=data.shape)
data += 0.5 * noise
data[::2] += 5 * noise[::2]
data[::4] += 20 * noise[::4]

X = data[:, 0][:, np.newaxis]
y = data[:, 1]

# Fit line using all data
model = linear_model.LinearRegression()
model.fit(X, y)

# Robustly fit linear model with RANSAC algorithm
model_robust = linear_model.LinearRegression()
n, inlier_mask = ransac(X, y, model_robust, min_n_samples=2,
                        residual_threshold=2)
outlier_mask = ~inlier_mask

# Generate coordinates of estimated models
line_X = np.arange(-250, 250)
line_y = model.predict(line_X[:, np.newaxis])
line_y_robust = model_robust.predict(line_X[:, np.newaxis])

plt.plot(data[inlier_mask, 0], data[inlier_mask, 1], '.g',
         label='Inlier data')
plt.plot(data[outlier_mask, 0], data[outlier_mask, 1], '.r',
         label='Outlier data')
plt.plot(line_X, line_y, '-k', label='Linear model from all data')
plt.plot(line_X, line_y_robust, '-b', label='Robustly fitted linear model')
plt.legend(loc='lower left')
plt.show()
