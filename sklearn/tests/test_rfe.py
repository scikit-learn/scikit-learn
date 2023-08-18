from sklearn.datasets import make_friedman1
from sklearn.feature_selection import RFE
from sklearn.svm import SVR
import pytest


def test_rfe():
    expected = [1, 1, 1, 1, 1, 6, 4, 3, 2, 5]
    X, y = make_friedman1(n_samples=50, n_features=10, random_state=0)
    estimator = SVR(kernel="linear")
    selector = RFE(estimator, n_features_to_select=5, step=1)
    selector_fit = selector.fit(X, y)
    assert selector_fit.support_.tolist() == [
        True,
        True,
        True,
        True,
        True,
        False,
        False,
        False,
        False,
        False,
    ]
    assert selector_fit.ranking_.tolist() == [1, 1, 1, 1, 1, 6, 4, 3, 2, 5]
