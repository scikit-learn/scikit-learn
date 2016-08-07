"""
=========================================================
Illustration of SelectDimensionKernel
=========================================================
A simple two-dimensional regression example computed in two different ways:
1. With a product of two RBF kernels on each feature.
2. With one anisotropic RBF kernels applied on both feature.

The figures illustrate the property of SelectDimensionKernel when applied on
different features of input data.
"""
print(__doc__)


# Authors: Behzad Tabibian <me@btabibian.com>
#
# License: BSD 3 clause


import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import SelectDimensionKernel, RBF

# Fixing seed
np.random.seed(1)

# Define a kernel with sum of two RBF kernels applied on individual dimentsion.
kernel = SelectDimensionKernel(RBF(length_scale=0.1), np.array([0])) * \
         SelectDimensionKernel(RBF(length_scale=0.4), np.array([1]))

# create GaussianProcessRegressor object.
gp = GaussianProcessRegressor(kernel=kernel)

# Training data
x_ = np.zeros((3, 2))
x_[:, 0] = np.linspace(0, 5, 3)
x_[:, 1] = np.linspace(0, 5, 3)

# Construct a mesh grid and compute synthetic data
x_0, x_1 = np.meshgrid(x_[:, 0], x_[:, 1])
x = np.vstack((x_0.flatten(), x_1.flatten())).T
y = (x[:, 0] - 2.5) + ((x[:, 1]-2.5)/2)
gp.fit(x, y)
y_out = gp.predict(x)

x_pred = np.zeros((15, 2))
x_pred[:, 0] = np.random.uniform(0, 5, 15)
x_pred[:, 1] = np.random.uniform(0, 5, 15)

y_pred = (x_pred[:, 0] - 2.5) + ((x_pred[:, 1]-2.5)/2)
y_out = gp.predict(x_pred)


fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 2, 1, projection='3d')

ax.plot_surface(x_0, x_1, y.reshape((x_.shape[0], x_.shape[0])),
                color='g', alpha=0.2)
ax.plot_wireframe(x_0, x_1, y.reshape((x_.shape[0], x_.shape[0])),
                  color='g', alpha=0.2)
ax.scatter(x_pred[:, 0], x_pred[:, 1], y_out, c='r', marker='o')
ax.scatter(x[:, 0], x[:, 1], y, c='b', marker='o')
ax.view_init(20, 60)

ax.set_xlabel('$X_0$')
ax.set_ylabel('$X_1$')
ax.set_zlabel('$Y$')

mae = np.abs(y_out-y_pred).mean()
ax.set_title('SelectDimensionKernel on two RBF, MAE: %.3f' % mae)
print("Absolute mean error %.4f" % mae)

# Create new RBF kernel applied on both features.
kernel = RBF(length_scale=[0.1, 0.4])
gp = GaussianProcessRegressor(kernel=kernel)

gp.fit(x, y)
y_out = gp.predict(x)


y_pred = np.sin((x_pred[:, 0] - 2.5)) + np.cos(((x_pred[:, 1]-2.5)/2))
y_pred = (x_pred[:, 0] - 2.5) + ((x_pred[:, 1]-2.5)/2)
y_out = gp.predict(x_pred)


ax = fig.add_subplot(1, 2, 2, projection='3d')

ax.plot_surface(x_0, x_1, y.reshape((x_.shape[0], x_.shape[0])),
                color='g', alpha=0.2)
ax.plot_wireframe(x_0, x_1, y.reshape((x_.shape[0], x_.shape[0])),
                  color='g', alpha=0.2)
ax.scatter(x_pred[:, 0], x_pred[:, 1], y_out, c='r', marker='o')
ax.scatter(x[:, 0], x[:, 1], y, c='b', marker='o')

ax.view_init(20, 60)

ax.set_xlabel('$X_0$')
ax.set_ylabel('$X_1$')
ax.set_zlabel('$Y$')

ax.set_title('Anisotropic RBF, MAE: %.3f' % np.abs(y_out-y_pred).mean())

mae = np.abs(y_out-y_pred).mean()
print("Absolute mean error %.4f" % mae)

plt.show()
