#!/usr/bin/env python
"""
======================
Least Angle Regression
======================

"""

# Author: Fabian Pedregosa <fabian.pedregosa@inria.fr>
# License: BSD Style.

from datetime import datetime
from itertools import cycle
import numpy as np
import pylab as pl

from scikits.learn import glm

n_samples, n_features = 500, 10

np.random.seed(0)
Y = np.random.randn(n_samples)
X = np.random.randn(n_samples, n_features)

from scikits.learn import datasets

diabetes = datasets.load_diabetes()
X = diabetes.data
Y = diabetes.target

################################################################################
# Fit models
################################################################################

################################################################################
# Demo path functions
################################################################################

eps = 1e-2 # the smaller it is the longer is the path

print "Computing regularization path using the LARS ..."
start = datetime.now()
clf = glm.LeastAngleRegression().fit(X, Y, normalize=True)
print "This took ", datetime.now() - start

alphas = -np.log10(clf.alphas_)

# # Display results
color_iter = cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k' , 'b', 'c', 'c'])

for coef_, color in zip(clf.coef_path_, color_iter):
    pl.plot(alphas, coef_.T, color)

ymin, ymax = pl.ylim()
pl.vlines(alphas, ymin, ymax, linestyle='dashed')
pl.xlabel('-Log(lambda)')
pl.ylabel('weights')
pl.title('Least Angle Regression (LAR) Paths')
pl.axis('tight')
pl.show()

