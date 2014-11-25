#!/usr/bin/python
# -*- coding: utf-8 -*-

r"""
============================================================================
Support Vector Regression: comparing different variants of Matern kernel
============================================================================

Support Vector Regression with four different variants of the Matern kernel
is compared on a (discontinuous) step-function:
 * The Matern kernel for coef0==1.5, learning a once differentiable function
 * The Matern kernel for coef0==2.5, learning a twice differentiable function
 * The Matern kernel for coef0==3.5, learning a three-rimes differentiable 
   function
 * The absolute-exponential kernel which corresponds to a Matern kernel
   with coef0==0.5
 * The squared-exponential (RBF) kernel which corresponds to a Matern kernel
   for the limit of coef0 becoming infinitely large

See Rasmussen and Williams 2006, pp84 for details regarding the different
variants of the Matern kernel.

The example shows that smaller values of coef0 can better approximate the
discontinuous step-function.
"""
print(__doc__)

# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# Licence: BSD 3 clause

from functools import partial

import numpy as np
from sklearn.svm import NuSVR
from sklearn.metrics.pairwise import matern_kernel

import matplotlib.pyplot as plt

np.random.seed(0)

# Train SVR with RBF and Matern kernels and plot resulting
# predictions
x = np.random.uniform(0, 10, 50)
y = (x < 5)

svr_rbf = NuSVR(nu=0.25, C=1e2, kernel="rbf", gamma=0.25)
svr_matern0_5 = NuSVR(nu=0.25, C=1e2,
                      kernel=partial(matern_kernel, coef0=0.5, gamma=0.25))
svr_matern1_5 = NuSVR(nu=0.25, C=1e2,
                      kernel=partial(matern_kernel, coef0=1.5, gamma=0.25))
svr_matern2_5 = NuSVR(nu=0.25, C=1e2,
                      kernel=partial(matern_kernel, coef0=2.5, gamma=0.25))
svr_matern3_5 = NuSVR(nu=0.25, C=1e2,
                      kernel=partial(matern_kernel, coef0=3.5, gamma=0.25))

svr_rbf.fit(x[:, None], y)
svr_matern0_5.fit(x[:, None], y)
svr_matern1_5.fit(x[:, None], y)
svr_matern2_5.fit(x[:, None], y)
svr_matern3_5.fit(x[:, None], y)

xp = np.linspace(0, 10, 100)
plt.scatter(x, y, c='k', s=25, zorder=10)
plt.plot(xp, xp < 5, label="True", c='k')
plt.plot(xp, svr_rbf.predict(xp[:, None]), label="RBF", c='g')
plt.plot(xp, svr_matern0_5.predict(xp[:, None]), label="Matern(0.5)", c='m')
plt.plot(xp, svr_matern1_5.predict(xp[:, None]), label="Matern(1.5)", c='r')
plt.plot(xp, svr_matern2_5.predict(xp[:, None]), label="Matern(2.5)", c='c')
plt.plot(xp, svr_matern3_5.predict(xp[:, None]), label="Matern(3.5)", c='b')
plt.legend(loc='best', title="kernel")
plt.xlabel("input")
plt.ylabel("target")
plt.show()
