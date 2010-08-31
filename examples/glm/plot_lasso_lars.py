#!/usr/bin/env python
"""
======================
Least Angle Regression
======================

"""

# Author: Fabian Pedregosa <fabian.pedregosa@inria.fr>
# License: BSD Style.

from datetime import datetime
import itertools
import numpy as np
import pylab as pl

from scikits.learn import glm
from scikits.learn import datasets

diabetes = datasets.load_diabetes()
X = diabetes.data
Y = diabetes.target
X[:, 6] = -X[:, 6]

################################################################################
# Demo path functions
################################################################################

print "Computing regularization path using the LARS ..."
start = datetime.now()
alphas, active, path = glm.lars(X, Y, method='lasso')
print "This took ", datetime.now() - start

alphas = np.sum(np.abs(path.T), axis=1)
alphas /= alphas[-1]

# # Display results
color_iter = itertools.cycle (['r', 'g', 'b', 'c'])

for coef_, color in zip(path, color_iter):
    pl.plot(alphas, coef_.T, color)

ymin, ymax = pl.ylim()
pl.vlines(alphas, ymin, ymax, linestyle='dashed')
pl.xlabel('-Log(lambda)')
pl.ylabel('weights')
pl.title('Least Angle Regression (LAR) Paths')
pl.axis('tight')
pl.show()

