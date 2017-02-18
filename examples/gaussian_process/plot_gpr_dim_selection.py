"""
=========================================================
Illustration of SelectDimensionKernel
=========================================================
A simple two-dimensional regression example computed in using two different kernels:
1. With a product of an RBF kernel and ExpSineSquared, where each is applied
   to one feature.
2. With one RBF kernel applied on both features.
"""

# Authors: Behzad Tabibian <me@btabibian.com>
#
# License: BSD 3 clause


from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import SelectDimensionKernel, RBF, ExpSineSquared

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

print(__doc__)
# Fixing seed
np.random.seed(0)

# Define a kernel with sum of two RBF kernels applied on individual dimentsion.
kernel = SelectDimensionKernel(RBF(), [0]) * \
         SelectDimensionKernel(ExpSineSquared(), [1])

# create GaussianProcessRegressor object.
gp_mix_kernel = GaussianProcessRegressor(kernel=kernel)
gp_rbf = GaussianProcessRegressor(kernel=RBF())

# Training data
x_ = np.zeros((8, 2))
x_[:, 0] = np.linspace(0, 5, 8)
x_[:, 1] = np.linspace(0, 5, 8)

# Construct a mesh grid and compute synthetic data
x_0, x_1 = np.meshgrid(x_[:, 0], x_[:, 1])
x = np.vstack((x_0.flatten(), x_1.flatten())).T
y = np.sin(x[:, 1]) + (x[:, 0]-2.5)/2 + np.random.normal(x.shape[0])

gp_rbf.fit(x[::2], y[::2])
y_pred_rbf = gp_rbf.predict(x[1::2])

gp_mix_kernel.fit(x[::2], y[::2])
y_pred_rbf_sin = gp_mix_kernel.predict(x[1::2])

mae = np.abs(y[1::2]-y_pred_rbf_sin).mean()
print("Absolute mean error for mixure of priodic kernel and RBF  %.4f" % mae)

mae = np.abs(y[1::2]-y_pred_rbf).mean()
print("Absolute mean error for RBF kernel %.4f" % mae)
