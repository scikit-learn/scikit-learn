
import numpy as np
import scipy.sparse as sp
from numpy.testing import assert_almost_equal

from scikits.learn import datasets
from scikits.learn.linear_model.ridge import Ridge
from scikits.learn.linear_model.sparse.ridge import Ridge as SparseRidge

diabetes = datasets.load_diabetes()

def test_diabetes():
    X, y = diabetes.data, diabetes.target
    X_sp = sp.csr_matrix(diabetes.data)

    #for fi in (True, False):
    for fi in (False,):
        ridge = Ridge(fit_intercept=fi)
        ridge_sp = SparseRidge(fit_intercept=fi)
        ridge.fit(X, y)
        ridge_sp.fit(X_sp, y)
        assert_almost_equal(ridge.score(X, y), ridge_sp.score(X_sp, y))

