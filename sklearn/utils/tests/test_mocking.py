import pytest
from sklearn.datasets import make_classification
from sklearn.utils._mocking import CheckingClassifier


X, y = make_classification(n_samples=20, random_state=42)


def test_check_X_on_predict_proba_success():
    def fail(X):
        return True

    clf = CheckingClassifier(check_X=fail).fit(X, y)
    clf.predict_proba(X)  # does not raise


def test_check_X_on_predict_proba_fail():
    def fail(X):
        return False

    clf = CheckingClassifier(check_X=fail)
    with pytest.raises(AssertionError):
        clf.predict_proba(X)
