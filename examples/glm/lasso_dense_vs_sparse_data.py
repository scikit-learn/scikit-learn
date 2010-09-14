"""
==============================
Lasso on dense and sparse data
==============================

We show that glm.Lasso and glm.sparse.Lasso
provide the same results.

XXX : At the end of the day it should also lead to a speed improvement

"""

from time import time
import numpy as np
from scipy import sparse
from scipy import linalg

from scikits.learn.glm.sparse import Lasso as SparseLasso
from scikits.learn.glm import Lasso as DenseLasso


###############################################################################
# The two Lasso implementation on Dense data
print "--- Dense matrices"

n_samples, n_features = 100, 10000
np.random.seed(0)
y = np.random.randn(n_samples)
X = np.random.randn(n_samples, n_features)

alpha = 1
sparse_lasso = SparseLasso(alpha=alpha, fit_intercept=False)
dense_lasso = DenseLasso(alpha=alpha, fit_intercept=False)

t0 = time()
sparse_lasso.fit(X, y, maxit=1000)
print "Sparse Lasso done in %fs" % (time() - t0)

t0 = time()
dense_lasso.fit(X, y, maxit=1000)
print "Dense Lasso done in %fs" % (time() - t0)

print "Distance between coefficients : %s" % linalg.norm(sparse_lasso.coef_
                                                        - dense_lasso.coef_)

###############################################################################
# The two Lasso implementation on Sparse data
print "--- Sparse matrices"

Xs = sparse.coo_matrix(X)
mask = Xs.data > 2 # Sparsify data matrix
col = Xs.col[mask]
row = Xs.row[mask]
Xs = Xs.tocsc()

print "Matrix density : %s %%" % (mask.sum() / float(X.size) * 100)

alpha = 0.1
sparse_lasso = SparseLasso(alpha=alpha, fit_intercept=False)
dense_lasso = DenseLasso(alpha=alpha, fit_intercept=False)

t0 = time()
sparse_lasso.fit(Xs, y, maxit=1000, tol=0.0)
print "Sparse Lasso done in %fs" % (time() - t0)

t0 = time()
dense_lasso.fit(Xs.todense(), y, maxit=1000, tol=0.0)
print "Dense Lasso done in %fs" % (time() - t0)

print "Distance between coefficients : %s" % linalg.norm(sparse_lasso.coef_
                                                        - dense_lasso.coef_)
