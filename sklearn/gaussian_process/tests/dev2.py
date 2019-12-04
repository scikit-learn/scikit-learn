from sklearn.gaussian_process import GaussianProcessRegressor as GP
from sklearn.gaussian_process.kernels import RBF
from matplotlib import pyplot as plt
import numpy as np
import GPy

"""
Illustrating a problem with test_y_normalization

"""

# Set amplitude of signal that we'll use as training data
A = 10

# Pick a kernel
RBF_params = {'length_scale' : 1.0}
kernel = RBF(**RBF_params)

# Test function
def f(x):
    return A * (x * np.sin(x))

# Test case (from test_gpr.py code)
X = np.atleast_2d([1., 3., 5., 6., 7., 8.]).T
y = f(X).ravel()
X_star = np.vstack(np.linspace(1, 8, 100))

# GP1 (y not normalised)
gp1 = GP(kernel=kernel, normalize_y=False)
gp1.fit(X, y)
y_pred1, y_pred_std1 = gp1.predict(X_star, return_std=True)

# GP2 (y normalised)
gp2 = GP(kernel=kernel, normalize_y=True)
gp2.fit(X, y)
y_pred2, y_pred_std2 = gp2.predict(X_star, return_std=True)

# GPy analysis for verification
kernel_gpy = GPy.kern.RBF(input_dim=1, lengthscale=1.)
gpy = GPy.models.GPRegression(X, np.vstack(y), kernel_gpy)
gpy.optimize()
y_pred_gpy, y_var_gpy = gpy.predict(X_star)
y_pred_std_gpy = np.sqrt(y_var_gpy)

# Plot some results
plt.figure()
plt.plot(X_star, y_pred1, 'k', label='sklearn, normalise_y=False')
plt.plot(X_star, y_pred1 + 3 * y_pred_std1, 'k')
plt.plot(X_star, y_pred1 - 3 * y_pred_std1, 'k')
plt.plot(X_star, y_pred2, 'r', label='sklearn, normalise_y=True')
plt.plot(X_star, y_pred2 + 3 * y_pred_std2, 'r')
plt.plot(X_star, y_pred2 - 3 * y_pred_std2, 'r')
plt.plot(X_star, y_pred_gpy, 'b', label='GPy')
plt.plot(X_star, y_pred_gpy + 3 * y_pred_std_gpy, 'b')
plt.plot(X_star, y_pred_gpy - 3 * y_pred_std_gpy, 'b')
plt.plot(X, y, 'k o')
plt.legend()
plt.show()
