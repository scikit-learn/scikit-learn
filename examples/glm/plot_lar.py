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


################################################################################
# Demo path functions
################################################################################

print "Computing regularization path using the LARS ..."
start = datetime.now()
clf = glm.LeastAngleRegression().fit(X, Y, normalize=True)
print "This took ", datetime.now() - start

alphas = -np.log10(clf.alphas_)

# # Display results
color_iter = itertools.cycle (['r', 'g', 'b', 'c'])

for coef_, color in zip(clf.coef_path_, color_iter):
    pl.plot(alphas, coef_.T, color)

ymin, ymax = pl.ylim()
pl.vlines(alphas, ymin, ymax, linestyle='dashed')
pl.xlabel('-Log(lambda)')
pl.ylabel('weights')
pl.title('Least Angle Regression (LAR) Paths')
pl.axis('tight')
pl.show()

