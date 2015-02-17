
import numpy as np
import pylab
from matplotlib.colors import LogNorm

from sklearn.gaussian_process import GaussianProcessRegression
from sklearn.gaussian_process.kernels import RBF, WhiteKernel


np.random.seed(0)
X = np.random.uniform(0, 5, 20)[:, None]
y = 0.5*np.sin(3*X[:, 0]) + np.random.normal(0, 0.5, X.shape[0])

rbf_kernel = (1e-10, 1.0, 100) * RBF(param_space=(1e-10, 1.0, None))
white_kernel = WhiteKernel(param_space=(1e-10, 1e-5, 1e+1))

gp = GaussianProcessRegression(kernel=rbf_kernel + white_kernel,
                               y_err=0.0).fit(X, y)

pylab.figure(0)
X_ = np.linspace(0, 5, 100)
y_mean, y_cov = gp.predict(X_[:, None], return_cov=True)
pylab.plot(X_, y_mean, 'k', lw=3, zorder=9)
pylab.fill_between(X_, y_mean - np.sqrt(np.diag(y_cov)),
                   y_mean + np.sqrt(np.diag(y_cov)),
                   alpha=0.5, color='k')

pylab.scatter(X[:, 0], y, c='r', s=50, zorder=10)

# Plot LML landscape
pylab.figure(1)
theta0 = np.logspace(-2, 3, 49)
theta1 = np.logspace(-2, 0, 50)
Theta0, Theta1 = np.meshgrid(theta0, theta1)
LML = [[gp.log_marginal_likelihood([0.36, Theta0[i, j], Theta1[i, j]])
        for i in range(Theta0.shape[0])] for j in range(Theta0.shape[1])]
LML = np.array(LML).T

vmin, vmax = (-LML).min(), (-LML).max()
vmax = 50
pylab.contour(Theta0, Theta1, -LML,
              levels=np.logspace(np.log10(vmin), np.log10(vmax), 50),
              norm=LogNorm(vmin=vmin, vmax=vmax))
pylab.xscale("log")
pylab.yscale("log")
pylab.colorbar()

pylab.show()