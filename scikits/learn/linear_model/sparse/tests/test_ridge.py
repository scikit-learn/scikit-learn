
import numpy as np
import scipy.sparse as sp
from numpy.testing import assert_almost_equal, assert_array_almost_equal

from scikits.learn import datasets
from scikits.learn.linear_model.ridge import Ridge
from scikits.learn.linear_model.sparse.ridge import Ridge as SpRidge
from scikits.learn.linear_model.sparse.ridge import RidgeClassifier
from scikits.learn.linear_model.sparse.ridge import RidgeClassifierLOO

diabetes = datasets.load_diabetes()
iris = datasets.load_iris()

X_diabetes = sp.csr_matrix(diabetes.data)
X_iris = sp.csr_matrix(iris.data)
y_diabetes = diabetes.target
y_iris = iris.target

def test_ridge_diabetes():
    #for fi in (True, False):
    for fi in (False,):
        ridge = Ridge(fit_intercept=fi)
        ridge.fit(diabetes.data, y_diabetes)
        for solver in ("cg", "default"):
            ridge_sp = SpRidge(fit_intercept=fi, solver=solver)
            ridge_sp.fit(X_diabetes, y_diabetes)
            assert_almost_equal(ridge.score(diabetes.data, y_diabetes),
                                ridge_sp.score(X_diabetes, y_diabetes))

def test_multi_ridge_diabetes():
    # simulate several responses
    Y = np.vstack((y_diabetes,y_diabetes)).T

    ridge = SpRidge(fit_intercept=False)
    ridge.fit(X_diabetes, Y)
    Y_pred = ridge.predict(X_diabetes)
    ridge.fit(X_diabetes, y_diabetes)
    y_pred = ridge.predict(X_diabetes)
    assert_array_almost_equal(np.vstack((y_pred,y_pred)).T,
                              Y_pred)

def test_ridge_classifiers():
    for clf in (RidgeClassifier(), RidgeClassifierLOO()):
        clf.fit(X_iris, y_iris)
        y_pred = clf.predict(X_iris)
        assert np.mean(y_iris == y_pred) >= 0.8

