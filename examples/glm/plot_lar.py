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

n_samples, n_features = 50, 10

np.random.seed(0)
y = np.random.randn(n_samples)
X = np.random.randn(n_samples, n_features)

################################################################################
# Fit models
################################################################################

################################################################################
# Demo path functions
################################################################################

eps = 1e-2 # the smaller it is the longer is the path

print "Computing regularization path using the LARS ..."
start = datetime.now()
clf = glm.LeastAngleRegression().fit(X, y, n_features=7)
print "This took ", datetime.now() - start

# alphas = np.append(clf.alphas_, np.zeros(7))
alphas = -np.log10(np.abs(clf.alphas_))

# Display results
color_iter = cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])

for coef_, color in zip(clf.coef_, color_iter):
    pl.plot(alphas, coef_.T, color)


pl.xlabel('-Log(lambda)')
pl.ylabel('weights')
pl.title('Least Angle Regression (LAR) Paths')
pl.axis('tight')
pl.show()

