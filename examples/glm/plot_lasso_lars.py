#!/usr/bin/env python
"""
=====================
Lasso path using LARS
=====================

Computes Lasso Path along the regularization parameter using the LARS
algorithm on the diabetest dataset.

"""
print __doc__

# Author: Fabian Pedregosa <fabian.pedregosa@inria.fr>
#         Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

import numpy as np
import pylab as pl

from scikits.learn import glm
from scikits.learn import datasets

diabetes = datasets.load_diabetes()
X = diabetes.data
y = diabetes.target
X[:,6] *= -1 # To reproduce wikipedia LASSO page

print "Computing regularization path using the LARS ..."
_, _, coefs_ = glm.lars_path(X, y, method='lasso', verbose=True)

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

