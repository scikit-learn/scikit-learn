import numpy as np
from numpy.testing import assert_array_equal
import pytest
from scipy import sparse

from sklearn.datasets import load_iris
from sklearn.utils import check_array

from sklearn.utils._mocking import CheckingClassifier


@pytest.fixture
def iris():
    return load_iris(return_X_y=True)


def test_checking_classifier(iris):
    # Check that the CheckingClassifier output what we expect
    X, y = iris
    clf = CheckingClassifier()
    clf.fit(X, y)

    assert_array_equal(clf.classes_, np.unique(y))
    assert clf.n_features_in_ == X.shape[1]

    y_pred = clf.predict(X)
    assert_array_equal(y_pred, np.zeros(y_pred.size, dtype=np.int))

    assert clf.score(X) == pytest.approx(0)
    clf.set_params(foo_param=10)
    assert clf.fit(X, y).score(X) == pytest.approx(1)


def test_checking_classifier_sparse(iris):
    # Smoke test to check that we can pass a sparse matrix when check are
    # disabled
    X, y = iris
    X_sparse = sparse.csr_matrix(X)

    clf = CheckingClassifier()
    clf.fit(X_sparse, y)
    y_pred = clf.predict(X_sparse)

    assert_array_equal(y_pred, np.zeros(y_pred.size, dtype=np.int))
    assert clf.score(X_sparse, y) == pytest.approx(0)


def test_checking_classifier_with_params(iris):
    X, y = iris
    X_sparse = sparse.csr_matrix(X)

    def check_X_is_sparse(X):
        if not sparse.issparse(X):
            raise ValueError("X is not sparse")
        return True

    clf = CheckingClassifier(check_X=check_X_is_sparse)
    with pytest.raises(ValueError, match="X is not sparse"):
        clf.fit(X, y)
    clf.fit(X_sparse, y)

    def _check_array(X, **params):
        check_array(X, **params)
        return True

    clf = CheckingClassifier(
        check_X=_check_array, check_X_params={"accept_sparse": False}
    )
    clf.fit(X, y)
    with pytest.raises(TypeError, match="A sparse matrix was passed"):
        clf.fit(X_sparse, y)


def test_checking_classifier_fit_params(iris):
    # check the error raised when the number of samples is not the one expected
    X, y = iris
    clf = CheckingClassifier(expected_fit_params=["sample_weight"])
    sample_weight = np.ones(len(X) // 2)

    with pytest.raises(AssertionError, match="Fit parameter sample_weight"):
        clf.fit(X, y, sample_weight=sample_weight)
