#!/usr/bin/env python
"""
=================================
Lasso with Least Angle Regression
=================================

"""

# Author: Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

from datetime import datetime
import itertools
import numpy as np
import pylab as pl

from scikits.learn import glm
from scikits.learn import datasets

diabetes = datasets.load_diabetes()
X = diabetes.data
y = diabetes.target
# someting's wrong with our dataset
X[:, 6] = -X[:, 6]


# m, n = 200, 200
# np.random.seed(0)
# X = np.random.randn(m, n)
# y = np.random.randn(m)


# _xmean = X.mean(0)
# _ymean = y.mean(0)
# X = X - _xmean
# y = y - _ymean
# _norms = np.apply_along_axis (np.linalg.norm, 0, X)
# nonzeros = np.flatnonzero(_norms)
# X[:, nonzeros] /= _norms[nonzeros]

################################################################################
# Demo path functions
################################################################################

G = np.dot(X.T, X)
print "Computing regularization path using the LARS ..."
start = datetime.now()
alphas, active, path = glm.lars_path(X, y, Gram=G, method='lasso', max_iter=12)
print "This took ", datetime.now() - start

alphas = np.sum(np.abs(path.T), axis=1)
alphas /= alphas[-1]

# # Display results
color_iter = itertools.cycle(['r', 'g', 'b', 'c'])

for coef_, color in zip(path, color_iter):
    pl.plot(alphas, coef_.T, color)

ymin, ymax = pl.ylim()
pl.vlines(alphas, ymin, ymax, linestyle='dashed')
pl.xlabel('-Log(lambda)') # XXX : wrong label
pl.ylabel('Coefficients')
pl.title('LASSO Path')
pl.axis('tight')
pl.show()

