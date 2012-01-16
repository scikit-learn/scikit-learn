"""
===================
Logistic Regression
===================

Comparison of the sparsity (percentage of zero coefficients) of solutions when
L1 and L2 penalty are used for different values of C. We can see that large
values of C give more freedom to the model.  Conversely, smaller values of C
constrain the model more. In the L1 penalty case, this leads to sparser
solutions.
"""

print __doc__

# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Mathieu Blondel <mathieu@mblondel.org>
# License: BSD Style.

import numpy as np

from sklearn.linear_model import LogisticRegression
from sklearn import datasets

# FIXME: the iris dataset has only 4 features!
iris = datasets.load_iris()
X = iris.data
y = iris.target

# Set regularization parameter
for C in (0.1, 1, 10):
    clf_l1_LR = LogisticRegression(C=C, penalty='l1')
    clf_l2_LR = LogisticRegression(C=C, penalty='l2')
    clf_l1_LR.fit(X, y)
    clf_l2_LR.fit(X, y)

    coef_l1_LR = clf_l1_LR.coef_[:]
    coef_l2_LR = clf_l2_LR.coef_[:]

    # coef_l1_LR contains zeros due to the
    # L1 sparsity inducing norm

    sparsity_l1_LR = np.mean(coef_l1_LR == 0) * 100
    sparsity_l2_LR = np.mean(coef_l2_LR == 0) * 100

    print "C=%f" % C
    print "Sparsity with L1 penalty: %f" % sparsity_l1_LR
    print "Sparsity with L2 penalty: %f" % sparsity_l2_LR
