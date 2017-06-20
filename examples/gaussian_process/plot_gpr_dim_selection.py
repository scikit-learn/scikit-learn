"""
=========================================================
Illustration of SelectDimensionKernel
=========================================================

This example models two-dimensional data where each feature is generated
independently and in a different fashion.

1. With a product of a RBF kernel and ExpSineSquared, where each is applied
   to one feature using the SelectDimensionKernel.

2. With one RBF kernel applied on both features.

We will show that the SelectDimensionKernel achieves better performance
on test data since it offers better flexibility in modelling.
"""

# License: BSD 3 clause

import numpy as np
from matplotlib import pyplot as plt

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ExpSineSquared
from sklearn.gaussian_process.kernels import RBF
from sklearn.gaussian_process.kernels import SelectDimensionKernel
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import train_test_split

print(__doc__)

# Fixing seed
np.random.seed(0)

# Define a kernel with sum of two RBF kernels applied on individual dimentsion.
kernel = (SelectDimensionKernel(RBF(length_scale=[1.0]), [0]) *
          SelectDimensionKernel(ExpSineSquared(length_scale=[1.0]), [1]))

# Training data
min_count = 50
max_count = 200
rbf_sin_err = []
rbf_only_err = []
param_range = np.arange(min_count, max_count, 10)
for count in param_range:
    X = np.zeros((count, 2))
    X[:, 0] = np.linspace(0, 40, count)
    X[:, 1] = np.linspace(0, 40, count)

    # Construct a mesh grid and compute synthetic data
    y = np.sin(2 * X[:, 1]) + (X[:, 0] - 20) / 2 + np.random.normal(X.shape[0],
                                                                    scale=2)

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.85)

    gp_mix_kernel = GaussianProcessRegressor(kernel=kernel)
    gp_rbf = GaussianProcessRegressor(kernel=RBF([1.0, 1.0]))

    gp_rbf.fit(X_train, y_train)
    y_pred_rbf = gp_rbf.predict(X_test)

    gp_mix_kernel.fit(X_train, y_train)
    y_pred_rbf_sin = gp_mix_kernel.predict(X_test)

    rbf_sin_err.append(mean_absolute_error(y_test , y_pred_rbf_sin))
    rbf_only_err.append(mean_absolute_error(y_test, y_pred_rbf))

plt.hold('on')
plt.plot((param_range * 0.15).astype(int),
         rbf_only_err, lw=3, label="RBF only")
plt.plot((param_range * 0.15).astype(int),
         rbf_sin_err, lw=3, label="RBF[0]*ExpSineSquared[1]")

plt.xlabel("Training Set size")
plt.ylabel("Mean Absolute Error on Test Set")
plt.legend()
plt.show()
