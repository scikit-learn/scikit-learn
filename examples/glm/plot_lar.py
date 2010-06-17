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

print "Computing regularization path using the LARS ..."
start = datetime.now()
clf = glm.LeastAngleRegression().fit(X, Y, n_features=9, normalize=False)
print "This took ", datetime.now() - start

# Display results
alphas = -np.log10(clf.alphas_[:-1])
pl.plot(alphas, clf.coef_.T)

ymin, ymax = pl.ylim()
pl.vlines(alphas, ymin, ymax, linestyle='dashed')
pl.xlabel('-Log(lambda)')
pl.ylabel('weights')
pl.title('Least Angle Regression (LAR) Paths')
pl.axis('tight')
pl.show()

