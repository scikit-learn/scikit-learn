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
alphas_, _, coefs_ = glm.lars_path(X, y, max_iter=9, method="lar")
print "This took ", datetime.now() - start

###############################################################################
# Display path
pl.plot(-np.log10(alphas_), coefs_.T)
ymin, ymax = pl.ylim()
pl.vlines(-np.log10(alphas_), ymin, ymax, linestyle='dashed')
pl.xlabel('-Log(lambda)') # XXX : wrong label
pl.ylabel('Coefficients')
pl.title('Least Angle Regression (LAR) Path')
pl.axis('tight')
pl.show()

