"""
=========================================================
Illustration of SelectDimensionKernel
=========================================================

A simple two-dimensional regression example computed using two
different kernels:
1. With a product of an RBF kernel and ExpSineSquared, where each is applied
   to one feature.
2. With one RBF kernel applied on both features.
"""

# Authors: Behzad Tabibian <me@btabibian.com>
#
# License: BSD 3 clause


from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ExpSineSquared
from sklearn.gaussian_process.kernels import SelectDimensionKernel


import numpy as np
from matplotlib import pyplot as plt

print(__doc__)
# Fixing seed
np.random.seed(0)

# Define a kernel with sum of two RBF kernels applied on individual dimentsion.


# Training data
min_count = 50
max_count = 200
rbf_sin_err = []
rbf_only_err = []
param_range = np.arange(min_count, max_count, 10)
for count in param_range:
    x = np.zeros((count, 2))
    x[:, 0] = np.linspace(0, 40, count)
    x[:, 1] = np.linspace(0, 40, count)

    # Construct a mesh grid and compute synthetic data
    y = np.sin(2 * x[:, 1]) + (x[:, 0] - 20) / 2 + np.random.normal(x.shape[0],
                                                                    scale=2)

    kernel = SelectDimensionKernel(RBF(), [0]) * \
        SelectDimensionKernel(ExpSineSquared(), [1])

    gp_mix_kernel = GaussianProcessRegressor(kernel=kernel)
    gp_rbf = GaussianProcessRegressor(kernel=RBF([1.0, 1.0]))

    shuffle_ind = np.arange(count)
    np.random.shuffle(shuffle_ind)
    training_ratio = 0.15
    x_training = x[shuffle_ind[:int(training_ratio * count)]]
    y_training = y[shuffle_ind[:int(training_ratio * count)]]

    x_test = x[shuffle_ind[int(training_ratio * count):]]
    y_test = y[shuffle_ind[int(training_ratio * count):]]

    gp_rbf.fit(x_training, y_training)
    y_pred_rbf = gp_rbf.predict(x_test)

    gp_mix_kernel.fit(x_training, y_training)
    y_pred_rbf_sin = gp_mix_kernel.predict(x_test)

    rbf_sin_err.append(np.abs(y_test - y_pred_rbf_sin).mean())
    rbf_only_err.append(np.abs(y_test - y_pred_rbf).mean())

plt.hold('on')
plt.plot((param_range * training_ratio).astype(int),
         rbf_only_err, lw=3, label="RBF only")
plt.plot((param_range * training_ratio).astype(int),
         rbf_sin_err, lw=3, label="RBF[0]*ExpSineSquared[1]")

plt.xlabel("Training Set size")
plt.ylabel("Mean Absolute Error on Test Set")
plt.legend()
plt.savefig("./plot_gpr_dim_selection.png")
plt.show()
