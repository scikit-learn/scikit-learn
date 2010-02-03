"""
Expectation maximation example
"""

import pylab
import numpy as np
from scikits.learn.machine.em.densities2 import gauss_ell

#=========================================
# Test plotting a simple diag 2d variance:
#=========================================
va  = np.array([5, 3])
mu  = np.array([2, 3])

# Generate a multivariate gaussian of mean mu and covariance va
X       = np.random.randn(2, 1e3)
Yc      = np.dot(np.diag(np.sqrt(va)), X)
Yc      = Yc.transpose() + mu

# Plotting
Xe, Ye  = gauss_ell(mu, va, npoints = 100)
pylab.figure()
pylab.plot(Yc[:, 0], Yc[:, 1], '.')
pylab.plot(Xe, Ye, 'r')

#=========================================
# Test plotting a simple full 2d variance:
#=========================================
va  = np.array([[0.2, 0.1],[0.1, 0.5]])
mu  = np.array([0, 3])

# Generate a multivariate gaussian of mean mu and covariance va
X       = np.random.randn(1e3, 2)
Yc      = np.dot(np.linalg.cholesky(va), X.transpose())
Yc      = Yc.transpose() + mu

# Plotting
Xe, Ye  = gauss_ell(mu, va, npoints = 100, level=0.95)
pylab.figure()
pylab.plot(Yc[:, 0], Yc[:, 1], '.')
pylab.plot(Xe, Ye, 'r')
pylab.show()
