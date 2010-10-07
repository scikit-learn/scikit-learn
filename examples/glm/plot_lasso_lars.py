#!/usr/bin/env python
"""
=================================
Lasso with Least Angle Regression
=================================

Computes Lasso Path with the LARS algorithm

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
X[:,6] *= -1 # To reproduce wikipedia LASSO page

################################################################################
# Demo path functions

print "Computing regularization path using the LARS ..."
start = datetime.now()
alphas_, _, coefs_ = glm.lars_path(X, y, method='lasso')
print "This took ", datetime.now() - start

xx = np.sum(np.abs(coefs_.T), axis=1)
xx /= xx[-1]
pl.plot(xx, coefs_.T)
ymin, ymax = pl.ylim()
pl.vlines(xx, ymin, ymax, linestyle='dashed')
pl.xlabel('|coef| / max|coef|')
pl.ylabel('Coefficients')
pl.title('LASSO Path')
pl.axis('tight')
pl.show()

