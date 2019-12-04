from sklearn.gaussian_process import GaussianProcessRegressor as GP
from sklearn.gaussian_process.kernels import RBF
from matplotlib import pyplot as plt
import numpy as np
import GPy

"""
Verifying the proposed solution against GPy. 

"""

# Set amplitude of signal that we'll use as training data
A = 10

# Pick a kernel
RBF_params = {'length_scale' : 1.0}
kernel = RBF(**RBF_params)

# Test function
def f(x):
    return A * np.sin(x)

# Test case (from test_gpr.py code)
X = np.atleast_2d([1., 3., 5., 6., 7., 8.]).T
y = f(X).ravel()
X_star = np.vstack(np.linspace(1, 8, 100))

# GP2 (y normalised)
gp = GP(kernel=kernel, normalize_y=True)
gp.fit(X, y)
y_pred, y_pred_std = gp.predict(X_star, return_std=True)

# GPy analysis for verification
kernel_gpy = GPy.kern.RBF(input_dim=1, lengthscale=1.)
gpy = GPy.models.GPRegression(X, np.vstack(y), kernel_gpy)
gpy.optimize()
y_pred_gpy, y_var_gpy = gpy.predict(X_star)
y_pred_std_gpy = np.sqrt(y_var_gpy)

# Plot some results
plt.figure()
plt.plot(X_star, y_pred, 'k', label='sklearn (proposed solution)')
plt.plot(X_star, y_pred + 3 * y_pred_std, 'k')
plt.plot(X_star, y_pred - 3 * y_pred_std, 'k')
plt.plot(X_star, y_pred_gpy, 'b', label='GPy')
plt.plot(X_star, y_pred_gpy + 3 * y_pred_std_gpy, 'b')
plt.plot(X_star, y_pred_gpy - 3 * y_pred_std_gpy, 'b')
plt.plot(X, y, 'k o')
plt.legend()
plt.show()
