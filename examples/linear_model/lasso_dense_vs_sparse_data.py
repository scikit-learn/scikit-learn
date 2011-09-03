"""
==============================
Lasso on dense and sparse data
==============================

We show that linear_model.Lasso and linear_model.sparse.Lasso
provide the same results and that in the case of
sparse data linear_model.sparse.Lasso improves the speed.

"""
print __doc__

from time import time
import numpy as np
from scipy import sparse
from scipy import linalg

from sklearn.linear_model.sparse import Lasso as SparseLasso
from sklearn.linear_model import Lasso as DenseLasso


###############################################################################
# The two Lasso implementations on Dense data
print "--- Dense matrices"

n_samples, n_features = 200, 10000
np.random.seed(0)
y = np.random.randn(n_samples)
X = np.random.randn(n_samples, n_features)

alpha = 1
sparse_lasso = SparseLasso(alpha=alpha, fit_intercept=False)
dense_lasso = DenseLasso(alpha=alpha, fit_intercept=False)

t0 = time()
sparse_lasso.fit(X, y, max_iter=1000)
print "Sparse Lasso done in %fs" % (time() - t0)

t0 = time()
dense_lasso.fit(X, y, max_iter=1000)
print "Dense Lasso done in %fs" % (time() - t0)

print "Distance between coefficients : %s" % linalg.norm(sparse_lasso.coef_
                                                        - dense_lasso.coef_)

###############################################################################
# The two Lasso implementations on Sparse data
print "--- Sparse matrices"

Xs = X.copy()
Xs[Xs < 2.5] = 0.0
Xs = sparse.coo_matrix(Xs)
Xs = Xs.tocsc()

print "Matrix density : %s %%" % (Xs.nnz / float(X.size) * 100)

alpha = 0.1
sparse_lasso = SparseLasso(alpha=alpha, fit_intercept=False)
dense_lasso = DenseLasso(alpha=alpha, fit_intercept=False)

t0 = time()
sparse_lasso.fit(Xs, y, max_iter=1000)
print "Sparse Lasso done in %fs" % (time() - t0)

t0 = time()
dense_lasso.fit(Xs.todense(), y, max_iter=1000)
print "Dense Lasso done in %fs" % (time() - t0)

print "Distance between coefficients : %s" % linalg.norm(sparse_lasso.coef_
                                                        - dense_lasso.coef_)
