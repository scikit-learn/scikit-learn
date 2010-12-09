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

from scikits.learn import linear_model
from scikits.learn import datasets

diabetes = datasets.load_diabetes()
X = diabetes.data
y = diabetes.target

print "Computing regularization path using the LARS ..."
alphas, _, coefs = linear_model.lars_path(X, y, method='lasso', verbose=True)

# trim the last alpha that is significantly far closer to zero on the log
# plot
m_log_alphas = -np.log(alphas)[:-1]
pl.plot(m_log_alphas, coefs.T[:-1, :])

ymin, ymax = pl.ylim()
pl.vlines(m_log_alphas, ymin, ymax, linestyle='dashed')
pl.xlabel('-log(lambda)')
pl.ylabel('Coefficients')
pl.title('LASSO Path')
pl.axis('tight')
pl.show()

