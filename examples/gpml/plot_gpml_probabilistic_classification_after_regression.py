#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
===============================================
Gaussian Processes for Machine Learning example
===============================================

A two-dimensional regression exercise with a post-processing
allowing for probabilistic classification thanks to the
Gaussian property of the prediction.
"""
# Author: Vincent Dubourg <vincent.dubourg@gmail.com
# License: BSD style

import numpy as np
from scipy import stats
from scikits.learn.gpml import GaussianProcessModel
from matplotlib import pyplot as pl
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
class FormatFaker(object):
    def __init__(self, str): self.str = str
    def __mod__(self, stuff): return self.str

# Standard normal distribution functions
Grv = stats.distributions.norm()
phi = lambda x: Grv.pdf(x)
PHI = lambda x: Grv.cdf(x)
PHIinv = lambda x: Grv.ppf(x)

# A few constants
beta0 = 8
b, kappa, e = 5, .5, .1

# The function to predict (classification will then consist in predicting wheter g(x) <= 0 or not)
g = lambda x: b - x[:,1] - kappa * ( x[:,0] - e )**2.

# Design of experiments
X = np.array([[-4.61611719, -6.00099547],
              [ 4.10469096,  5.32782448],
              [ 0.00000000, -0.50000000],
              [-6.17289014, -4.6984743 ],
              [ 1.3109306 , -6.93271427],
              [-5.03823144,  3.10584743],
              [-2.87600388,  6.74310541],
              [ 5.21301203,  4.26386883]])

# Observations
Y = g(X)

# Instanciate and fit Gaussian Process Model
aGaussianProcessModel = GaussianProcessModel(theta0=5e-1)

# Don't perform MLE or you'll get a perfect prediction for this simple example!
aGaussianProcessModel.fit(X, Y)

# Evaluate real function, the prediction and its MSE on a grid
res = 50
x1, x2 = np.meshgrid(np.linspace(-beta0, beta0, res), np.linspace(-beta0, beta0, res))
xx = np.vstack([x1.reshape(x1.size), x2.reshape(x2.size)]).T

YY = g(xx)
yy, MSE = aGaussianProcessModel.predict(xx, eval_MSE=True)
sigma = np.sqrt(MSE)
yy = yy.reshape((res,res))
YY = YY.reshape((res,res))
sigma = sigma.reshape((res,res))
k = PHIinv(.975)

# Plot the probabilistic classification iso-values using the Gaussian property of the prediction
fig = pl.figure(1)
ax = fig.add_subplot(111)
ax.axes.set_aspect('equal')
cax = pl.imshow(np.flipud(PHI(-yy/sigma)), cmap=cm.gray_r, alpha=.8, extent=(-beta0,beta0,-beta0,beta0))
norm = pl.matplotlib.colors.Normalize(vmin=0., vmax=0.9)
cb = pl.colorbar(cax, ticks=[0.,0.2,0.4,0.6,0.8,1.], norm=norm)
cb.set_label('${\\rm \mathbb{P}}\left[\widehat{G}(\mathbf{x}) \leq 0\\right]$')
pl.plot(X[Y <= 0, 0], X[Y <= 0, 1], 'r.', markersize=12)
pl.plot(X[Y > 0, 0], X[Y > 0, 1], 'b.', markersize=12)
pl.xticks([])
pl.yticks([])
ax.set_xticklabels([])
ax.set_yticklabels([])
pl.xlabel('$x_1$')
pl.ylabel('$x_2$')
print u'Click to place label: $g(\mathbf{x})=0$'
cs = pl.contour(x1, x2, YY, [0.], colors='k', linestyles='dashdot')
pl.clabel(cs,fmt=FormatFaker(u'$g(\mathbf{x})=0$'),fontsize=11,manual=True)
print u'Click to place label: ${\\rm \mathbb{P}}\left[{\widehat{G}}(\mathbf{x}) \leq 0\\right]= 2.5\%$'
cs = pl.contour(x1, x2, PHI(-yy/sigma), [0.025], colors='b', linestyles='solid')
pl.clabel(cs,fmt=FormatFaker(u'${\\rm \mathbb{P}}\left[{\widehat{G}}(\mathbf{x}) \leq 0\\right]= 2.5\%$'),fontsize=11,manual=True)
print u'Click to place label: $\mu_{\widehat{G}}(\mathbf{x})=0$'
cs = pl.contour(x1, x2, PHI(-yy/sigma), [0.5], colors='k', linestyles='dashed')
pl.clabel(cs,fmt=FormatFaker(u'$\mu_{\widehat{G}}(\mathbf{x})=0$'),fontsize=11,manual=True)
print u'Click to place label: ${\\rm \mathbb{P}}\left[{\widehat{G}}(\mathbf{x}) \leq 0\\right]= 97.5\%$'
cs = pl.contour(x1, x2, PHI(-yy/sigma), [0.975], colors='r', linestyles='solid')
pl.clabel(cs,fmt=FormatFaker(u'${\\rm \mathbb{P}}\left[{\widehat{G}}(\mathbf{x}) \leq 0\\right]= 97.5\%$'),fontsize=11,manual=True)

# Plot the prediction and the bounds of the 95% confidence interval
fig = pl.figure(2)
ax = Axes3D(fig)
ax.axes.set_aspect('equal')
ax.plot_surface(x1, x2, yy, linewidth = 0.5, rstride = 1, cstride = 1, color = 'k', alpha = .8)
ax.plot_surface(x1, x2, yy - k*sigma, linewidth = 0.5, rstride = 1, cstride = 1, color = 'b', alpha = .8)
ax.plot_surface(x1, x2, yy + k*sigma, linewidth = 0.5, rstride = 1, cstride = 1, color = 'r', alpha = .8)
ax.scatter3D(X[Y <= 0, 0], X[Y <= 0, 1], Y[Y <= 0], 'r.', s = 20)
ax.scatter3D(X[Y > 0, 0], X[Y > 0, 1], Y[Y > 0], 'b.', s = 20)
ax.set_xlabel(u'$x_1$')
ax.set_ylabel(u'$x_2$')
ax.set_zlabel(u'$\widehat{G}(x_1,\,x_2)$')

pl.show()