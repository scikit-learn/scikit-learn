"""
Testing Recursive feature elimination
"""

import numpy as np
from numpy.testing import assert_array_almost_equal

from ..rfe import RFE, RFECV
from ...datasets import load_iris
from ...metrics import zero_one
from ...svm import SVC
from ...utils import check_random_state


def test_rfe():
    generator = check_random_state(0)

    iris = load_iris()
    X = np.c_[iris.data, generator.normal(size=(len(iris.data), 6))]
    y = iris.target

    clf = SVC(kernel="linear")
    rfe = RFE(estimator=clf, n_features_to_select=4, step=0.1)
    rfe.fit(X, y)
    X_r = rfe.transform(X)

    assert X_r.shape == iris.data.shape
    assert_array_almost_equal(X_r[:10], iris.data[:10])

    assert_array_almost_equal(rfe.predict(X), clf.predict(iris.data))
    assert rfe.score(X, y) == clf.score(iris.data, iris.target)


def test_rfecv():
    generator = check_random_state(0)

    iris = load_iris()
    X = np.c_[iris.data, generator.normal(size=(len(iris.data), 6))]
    y = iris.target

    # Test using the score function
    rfecv = RFECV(estimator=SVC(kernel="linear"), step=1, cv=3)
    rfecv.fit(X, y)
    X_r = rfecv.transform(X)

    assert X_r.shape == iris.data.shape
    assert_array_almost_equal(X_r[:10], iris.data[:10])

    # Test using a customized loss function
    rfecv = RFECV(estimator=SVC(kernel="linear"), step=1, cv=3,
            loss_func=zero_one)
    rfecv.fit(X, y)
    X_r = rfecv.transform(X)

    assert X_r.shape == iris.data.shape
    assert_array_almost_equal(X_r[:10], iris.data[:10])
