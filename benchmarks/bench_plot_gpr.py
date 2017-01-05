"""
Timing Benchmark Test for GaussianProcessRegressor.

This Benchmark Test is used to evaluate the effect of different values
of n_jobs parameter on the performance of the regression algorithm.
"""
import numpy as np
import matplotlib.pyplot as plt

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels \
    import RBF, ConstantKernel as C, WhiteKernel

from sklearn.datasets import make_regression

import time


fixed_kernel = RBF(length_scale=1.0, length_scale_bounds="fixed")
kernels = [RBF(length_scale=1.0), fixed_kernel,
           RBF(length_scale=1.0, length_scale_bounds=(1e-3, 1e3)),
           C(1.0, (1e-2, 1e2)) *
           RBF(length_scale=1.0, length_scale_bounds=(1e-3, 1e3)),
           C(1.0, (1e-2, 1e2)) *
           RBF(length_scale=1.0, length_scale_bounds=(1e-3, 1e3)) +
           C(1e-5, (1e-5, 1e2)),
           C(0.1, (1e-2, 1e2)) *
           RBF(length_scale=1.0, length_scale_bounds=(1e-3, 1e3)) +
           C(1e-5, (1e-5, 1e2))]


X, y = make_regression(n_samples=500, n_features=5, n_informative=5,
                       noise=1.0)
x1 = ['1', '2', '3', '4', '5', '6']
y1 = []
y2 = []
y3 = []

for kernel in kernels:
    t1 = time.clock()
    GaussianProcessRegressor(kernel=kernel, n_jobs=1,
                             n_restarts_optimizer=5,
                             random_state=42).fit(X, y)
    t2 = time.clock()
    GaussianProcessRegressor(kernel=kernel, n_jobs=2,
                             n_restarts_optimizer=5,
                             random_state=42).fit(X, y)
    t3 = time.clock()
    GaussianProcessRegressor(kernel=kernel, n_jobs=-1,
                             n_restarts_optimizer=5,
                             random_state=42).fit(X, y)
    t4 = time.clock()

    mx = max(t4 - t3, t3 - t2, t2 - t1)
    y1.append((t2 - t1)/mx)
    y2.append((t3 - t2)/mx)
    y3.append((t4 - t3)/mx)

plt.xlabel('kernels')
plt.ylabel('Time Elapsed / Maximum Time Elapsed')
plt.title('Timing Benchmark')
plt.plot(x1, y1, 'ro-', label='n_jobs=1')
plt.plot(x1, y2, 'bs-', label='n_jobs=2')
plt.plot(x1, y3, 'g^-', label='n_jobs=-1')
plt.axis([1, 6, 0, 2])
plt.legend()
plt.show()
