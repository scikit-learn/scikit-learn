#!/usr/bin/env python
"""
============================
Least Angle Regression (LAR)
============================

Compute LAR path on diabetes dataset.

See: http://en.wikipedia.org/wiki/Least-angle_regression

"""
print __doc__

# Author: Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

from datetime import datetime
import numpy as np
import pylab as pl

from scikits.learn import glm
from scikits.learn import datasets

diabetes = datasets.load_diabetes()
X = diabetes.data
y = diabetes.target

X[:,6] *= -1 # To reproduce wikipedia LAR page

################################################################################
# Compute path functions

print "Computing regularization path using the LARS ..."
start = datetime.now()
_, _, coefs_ = glm.lars_path(X, y, max_features=10, method="lar")
print "This took ", datetime.now() - start

###############################################################################
# Display path
xx = np.sum(np.abs(coefs_), axis=0)
xx /= xx[-1]
pl.plot(xx, coefs_.T)
ymin, ymax = pl.ylim()
pl.vlines(xx, ymin, ymax, linestyle='dashed')
pl.xlabel('|coef| / max|coef|')
pl.ylabel('Coefficients')
pl.title('Least Angle Regression (LAR) Path')
pl.axis('tight')
pl.show()

