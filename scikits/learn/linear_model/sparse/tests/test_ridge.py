
import numpy as np
import scipy.sparse as sp
from numpy.testing import assert_almost_equal, assert_array_almost_equal

from scikits.learn import datasets
from scikits.learn.linear_model.ridge import Ridge
from scikits.learn.linear_model.sparse.ridge import Ridge as SpRidge
from scikits.learn.linear_model.sparse.ridge import RidgeClassifier

diabetes = datasets.load_diabetes()
iris = datasets.load_iris()

def test_ridge_diabetes():
    X, y = diabetes.data, diabetes.target

    #for fi in (True, False):
    for fi in (False,):
        ridge = Ridge(fit_intercept=fi)
        ridge.fit(X, y)
        for solver in ("cg", "default"):
            ridge_sp = SpRidge(fit_intercept=fi, solver=solver)
            ridge_sp.fit(X, y)
            assert_almost_equal(ridge.score(X, y), ridge_sp.score(X, y))

def test_multi_ridge_diabetes():
    X, y = diabetes.data, diabetes.target

    # simulate several responses
    Y = np.vstack((y,y)).T

    ridge = SpRidge(fit_intercept=False)
    ridge.fit(X, Y)
    Y_pred = ridge.predict(X)
    ridge.fit(X, y)
    y_pred = ridge.predict(X)
    assert_array_almost_equal(np.vstack((y_pred,y_pred)).T,
                              Y_pred)

def test_ridge_classifiers():
    clf = RidgeClassifier()
    clf.fit(iris.data, iris.target)
    y_pred = clf.predict(iris.data)
    assert np.mean(iris.target == y_pred) >= 0.8
